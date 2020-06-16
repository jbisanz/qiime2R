#' summarize_taxa
#'
#' Takes a feature table (ex OTUs/SVs) and a parsed taxonomy table to generate a named list of feature tables with abundances summed to each taxonomic level.
#'
#' @param features a table of abundances where rows are features and columns are samples
#' @param taxonomy a table of taxonomies that has been parsed (ex using parse_taxonomy). Should contain the following columns: Kingdom,Phylum,Class,Order,Family,Genus,Species.
#' @return a named list of abundances aggregated at each level.
#' @examples \dontrun{taxasums<-summarize_taxa(svtable, taxonomy); taxasums$Species}
#' @export
#'



summarize_taxa <- function(features, taxonomy) {

taxlevels<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  
if(missing(features)){stop("Feature table not provided")}
if(missing(taxonomy)){stop("taxonomy table not provided")}
if(sum(colnames(taxonomy) %in% taxlevels)!=7){stop("Taxonomy does not contain expected columns containing Kingdom,Phylum,Class,Order,Family,Genus,Species.")}  

output<-list()

for(lvl in taxlevels){
suppressMessages(
  output[[lvl]]<-
    features %>%
    as.data.frame() %>%
    rownames_to_column("FeatureID") %>%
    gather(-FeatureID, key="SampleID", value="Counts") %>%
    left_join(
      taxonomy %>% 
      rownames_to_column("FeatureID") %>%
      unite("Taxon", taxlevels[1:grep(lvl, taxlevels)], sep="; ") %>%
      select(FeatureID, Taxon)
    ) %>%
    group_by(SampleID, Taxon) %>%
    summarize(Counts=sum(Counts)) %>%
    ungroup() %>%
    spread(key=SampleID, value=Counts) %>%
    as.data.frame() %>%
    column_to_rownames("Taxon")
  )
}
return(output)
}
