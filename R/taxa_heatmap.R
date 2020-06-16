#' taxa_heatmap
#'
#' @description Uses ggplot2 to create a heatmap, for example on phylum level abundances. The most abundant features (defaults to 10, based on rowMeans) will be plotted unless user specified.
#' @param features Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param metadata A Table of metadata where sample names are the row names OR there is already a column with the name "SampleID"
#' @param category A metadata category to block samples by which is a column in the metadata (faceting via ggplot2)
#' @param normalize How should table be normalized?  default="log10(percent)", options=c("log10(percent)","clr","none"). Note if using  log10(percent), a pseudocount 1 read is added before conversion to percent
#' @param ntoplot A number of features to plot.
#' @return Barplot
#' @usage taxa_heatmap(features, metadata, "treatment")
#' @export

taxa_heatmap<-function(features, metadata, category, normalize, ntoplot){
  
  if(missing(ntoplot) & nrow(features)>30){ntoplot=30} else if (missing(ntoplot)){ntoplot=nrow(features)}
  if(missing(normalize)){normalize<-"log10(percent)"}
  if(normalize=="log10(percent)"){features<-log10(make_percent(features+1))} else if(normalize=="clr"){features<-make_clr(features)}
  
  if(missing(metadata)){metadata<-data.frame(SampleID=colnames(features))}
  if(!"SampleID" %in% colnames(metadata)){metadata <- metadata %>% rownames_to_column("SampleID")}
  if(!missing(category)){
    if(!category %in% colnames(metadata)){message(stop(category, " not found as column in metdata"))}
  }
  
  plotfeats<-names(sort(rowMeans(features), decreasing = TRUE)[1:ntoplot]) # extract the top N most abundant features on average
  
  roworder<-hclust(dist(features[plotfeats,]))
  roworder<-roworder$labels[roworder$order]
  
  colorder<-hclust(dist(t(features[plotfeats,])))
  colorder<-colorder$labels[colorder$order]
  
  suppressMessages(
  suppressWarnings(
  fplot<-
    features %>%
    as.data.frame() %>%
    rownames_to_column("Taxon") %>%
    gather(-Taxon, key="SampleID", value="Abundance") %>%
    filter(Taxon %in% plotfeats) %>%
    mutate(Taxon=factor(Taxon, levels=rev(plotfeats))) %>%
    left_join(metadata) %>%
    mutate(Taxon=factor(Taxon, levels=roworder)) %>%
    mutate(SampleID=factor(SampleID, levels=colorder))
  ))

  bplot<-
    ggplot(fplot, aes(x=SampleID, y=Taxon, fill=Abundance)) +
    geom_tile(stat="identity") +
    theme_q2r() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    coord_cartesian(expand=FALSE) +
    xlab("Sample") +
    ylab("Feature") +
    scale_fill_viridis_c()
  
  if(!missing(category)){bplot<-bplot + facet_grid(~get(category), scales="free_x", space="free")}
  
  return(bplot)
}
