#' Normalization Functions
#'
#' Takes a numeric table and transforms to a number of options as below. Note: Rows are expected to be features and columns are expected to be sample IDs.
#' * `make_clr()` A centered log2-ratio (columns sum to ~0). See and cite 10.1139/cjm-2015-0821
#' * `make_percent()` Convert to percent (columns sum to 100)
#' * `make_proportion()` Convert to proportion (columns sum to 1)
#' 
#' @param features Table of feature/OTU/SV counts where Samples are columns, and IDs are row names. Expects matrix or data.frame.
#' @param prior Only relevant to CLR: A numeric value to add before log transformation (default 0.5). *Ignored if CZM=TRUE
#' @param czm Only relevant to CLR: Should count zero multiplicative method be used instead of adding simple prior (TRUE/FALSE, defaults to FALSE)
#' @return Table of normalized abundances
#' @export

make_clr<-function(features, prior, czm){

  if(missing(czm)){czm=FALSE}
  if(missing(prior)){prior=0.5}
  
  if(czm){
    features <- t(cmultRepl(t(features),  label=0, method="CZM", output="p-counts", suppress.print=T))
  } else {
    features<-features + prior
  }
  
  features<-apply(features,2, function(column){log2(column)-mean(log2(column))})
  return(features)
}

make_percent<-function(features){
  features<-apply(features,2, function(x){100*(x/sum(x))})
  return(features)
}

make_proportion<-function(features){
  features<-apply(features,2, function(x){(x/sum(x))})
  return(features)
}