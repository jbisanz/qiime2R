#' Parse Q2 taxonomy
#'
#' @param taxonomy the imported taxonomy object (the result of read_qza() containing columns Feature.ID, Taxon, and Confidence)
#' @param tax_sep The separator between taxonomic levels. Defaults to one compatible with both GreenGenes and SILVA ("; " OR ";")
#' @param trim_extra Remove leading characters from taxonomic levels: ex: k__ or D_0__. TRUE/FALSE. default=TRUE 
#' 
#' @return a data.frame with feature IDs as row names and the columns: Kingdom, Phylum, Class, Order, Family, Genus, Species
#'
#' @examples \dontrun{taxonomy<-parse_taxonomy(taxonomy_artifact)}
#' @export
#'


parse_taxonomy <- function(taxonomy, tax_sep, trim_extra){
  if(missing(taxonomy)){stop("Taxonomy Table not supplied.")}
  if(missing(trim_extra)){trim_extra=TRUE}
  if(missing(tax_sep)){tax_sep="; |;"}
  if(sum(colnames(taxonomy) %in% c("Feature.ID","Taxon","Confidence", "Consensus"))!=3){stop("Table does not match expected format (colnames(obj) are Feature.ID, Taxon, (Confidence OR Consensus))")}

  taxonomy$Confidence<-NULL
  if(trim_extra){
  taxonomy$Taxon<-gsub("[kpcofgs]__","", taxonomy$Taxon) #remove leading characters from GG
  taxonomy$Taxon<-gsub("D_\\d__","", taxonomy$Taxon) #remove leading characters from SILVA
  }
  taxonomy<-suppressWarnings(taxonomy %>% separate(Taxon, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep=tax_sep))
  taxonomy<-apply(taxonomy, 2, function(x) if_else(x=="", NA_character_, x)) 
  taxonomy<-as.data.frame(taxonomy)
  rownames(taxonomy)<-taxonomy$Feature.ID
  taxonomy$Feature.ID<-NULL
  return(taxonomy)  
}
