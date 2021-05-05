#' Filter low abundance features from a table.
#'
#' Takes a table of features (samples are columns) and remove features that do not appear at in at least X samples OR do not have at least Y reads.
#'
#' @param features Table of feature/OTU/ASV counts where Samples are columns and feature IDs are row names
#' @param minsamples A minimum number of samples (>=) for the feature/OTU/ASV to be observed in. (Default=2)
#' @param minreads A minimum number of reads  (>=) for a feature/OTU/ASV to be kept across all samples. (Default=2)
#' @param verbose Should summary be printed? (T/F)
#' @return filtered feature table
#' @export

filter_features<-function(features, minsamples, minreads, verbose){
  if(missing(verbose)){verbose=TRUE}
  if(missing(minsamples)){minsamples=2}
  if(missing(minreads)){minreads=2}
  
  
  if(verbose==T){message(paste("Filtering features such that they are present in at least", minsamples, "samples with a total of at least", minreads, "reads."))}

  failreads<-rownames(features)[!rowSums(features)>=minreads]
  failsamples<-rownames(features)[(apply(features, 1, function(x) sum(x>0))>=minsamples)==FALSE]
  
  filtered<-features[!rownames(features) %in% c(failreads, failsamples),]
  
  
  if(verbose==T){message("...after filtering there are ", nrow(filtered), " of ", nrow(features)," features retaining ", round(sum(filtered)/sum(features), 4)*100,"% of reads.")}

  return(filtered)
}
