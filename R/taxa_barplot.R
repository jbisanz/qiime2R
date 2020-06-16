#' taxa_barplot
#'
#' @description Uses ggplot2 to create a stacked barplot, for example on phylum level abundances. The most abundant features (defaults to 10, based on rowMeans) will be plotted unless user specified.
#' @param features Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param metadata A Table of metadata where sample names are the row names OR there is already a column with the name "SampleID"
#' @param category A metadata category to block samples by which is a column in the metadata (faceting via ggplot2)
#' @param normalize How should table be normalized? (default="percent", options=c("percent","proportion","none"))
#' @param ntoplot A number of features to plot.
#' @return Barplot
#' @usage taxa_barplot(features, metadata, "treatment")
#' @export

taxa_barplot<-function(features, metadata, category, normalize, ntoplot){
  
  q2r_palette<-c(
    "blue4",
    "olivedrab",
    "firebrick",
    "gold",
    "darkorchid",
    "steelblue2",
    "chartreuse1",
    "aquamarine",
    "yellow3",
    "coral",
    "grey"
  )
  
  if(missing(ntoplot) & nrow(features)>10){ntoplot=10} else if (missing(ntoplot)){ntoplot=nrow(features)}
  if(missing(normalize)){normalize<-"percent"}
  if(normalize=="percent"){features<-make_percent(features)} else if(normalize=="proportion"){features<-make_proportion(features)}
  
  if(missing(metadata)){metadata<-data.frame(SampleID=colnames(features))}
  if(!"SampleID" %in% colnames(metadata)){metadata <- metadata %>% rownames_to_column("SampleID")}
  if(!missing(category)){
    if(!category %in% colnames(metadata)){message(stop(category, " not found as column in metdata"))}
  }
  
  plotfeats<-names(sort(rowMeans(features), decreasing = TRUE)[1:ntoplot]) # extract the top N most abundant features on average
  
  suppressMessages(
  suppressWarnings(
  fplot<-
    features %>%
    as.data.frame() %>%
    rownames_to_column("Taxon") %>%
    gather(-Taxon, key="SampleID", value="Abundance") %>%
    mutate(Taxon=if_else(Taxon %in% plotfeats, Taxon, "Remainder")) %>%
    group_by(Taxon, SampleID) %>%
    summarize(Abundance=sum(Abundance)) %>%
    ungroup() %>%
    mutate(Taxon=factor(Taxon, levels=rev(c(plotfeats, "Remainder")))) %>%
    left_join(metadata)
  ))

  
  bplot<-
    ggplot(fplot, aes(x=SampleID, y=Abundance, fill=Taxon)) +
    geom_bar(stat="identity") +
    theme_q2r() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    coord_cartesian(expand=FALSE) +
    xlab("Sample") +
    ylab("Abundance")
  
  if(ntoplot<=10){bplot<-bplot+scale_fill_manual(values=rev(q2r_palette), name="")}
  
  if(!missing(category)){bplot<-bplot + facet_grid(~get(category), scales="free_x", space="free")}
  
  return(bplot)
}
