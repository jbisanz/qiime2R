#' Parse Q2 taxonomy
#'
#' @param taxonomy a table-like object containing the columns Feature.ID and Taxon. Can be imported using read_qza(file)$data.
#' @param tax_sep The separator between taxonomic levels. Defaults to one compatible with both GreenGenes and SILVA ("; " OR ";")
#' @param trim_extra Remove leading characters from taxonomic levels: ex: k__ or D_0__. TRUE/FALSE. default=TRUE 
#' 
#' Note: Assumes an assignment has been made to all levels. Fills missing assignments with NA.
#' @return a data.frame with feature IDs as row names and the columns: Kingdom, Phylum, Class, Order, Family, Genus, Species
#'
#' @examples \dontrun{taxonomy<-parse_taxonomy(taxonomy)}
#' @export
#'


parse_taxonomy <- function(taxonomy, tax_sep, trim_extra){
  if(missing(taxonomy)){stop("Taxonomy Table not supplied.")}
  if(missing(trim_extra)){trim_extra=TRUE}
  if(missing(tax_sep)){tax_sep="; |;"}
  if(sum(colnames(taxonomy) %in% c("Feature.ID","Taxon"))!=2){stop("Table does not match expected format. ie does not have columns Feature.ID and Taxon.")}

  taxonomy<-taxonomy[,c("Feature.ID","Taxon")]
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


#' Parse Q2 taxonomy custom
#'
#' @param taxonomy a table-like object containing the columns Feature.ID and Taxon. Can be imported using read_qza(file)$data.
#' @param ranks a vector of rank labels expected in your table. Must include all expected ranks or extras will be dropped silently!. default=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")                
#' @param tax_sep The separator between taxonomic levels. Set as NULL if not required. Defaults to one compatible with both GreenGenes and SILVA ("; " OR ";").
#' @param trim_extra A regex string compatible with gsub to remove leading characters from taxonomic levels: ex: k__ or D_0__. default="[kpcofgs]__|D_|d__"
#' 
#' Note: Assumes an assignment has been made to all levels. Fills missing assignments with NA.
#' @return a data.frame with feature IDs as row names and the columns specified by the ranks parameter
#'
#' @examples \dontrun{taxonomy<-parse_taxonomy(taxonomy)}
#' @export
#'                  
parse_taxonomy_custom <- function(taxonomy, ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), tax_sep = "; |;", trim_extra = "[kpcofgs]__|D_|d__"){
  if(missing(taxonomy)){stop("Taxonomy Table not supplied.")}
  if(sum(colnames(taxonomy) %in% c("Feature.ID","Taxon"))!=2){stop("Table does not match expected format. ie does not have columns Feature.ID and Taxon.")}

  taxonomy<-taxonomy[,c("Feature.ID","Taxon")]
  if(!is.null(trim_extra)){
  taxonomy$Taxon<-gsub(trim_extra,"", taxonomy$Taxon) #remove leading characters from database [kpcofgs]__ for GG, D_ and d__ for SILVA
  }

  expected_ranks<-str_count(taxonomy$Taxon, tax_sep) %>% max
  if(length(ranks) < expected_ranks){warning(paste0("Expected ", expected_ranks, " ranks, supplied ", length(ranks), " (", paste0(ranks, collapse=","), "). Check tax_sep and ranks are appropriate."))}
  
  taxonomy<-suppressWarnings(taxonomy %>% separate(Taxon, ranks, sep=tax_sep))
  taxonomy<-apply(taxonomy, 2, function(x) if_else(x=="", NA_character_, x)) 
  taxonomy<-as.data.frame(taxonomy)
  rownames(taxonomy)<-taxonomy$Feature.ID
  taxonomy$Feature.ID<-NULL
  return(taxonomy)  
}
