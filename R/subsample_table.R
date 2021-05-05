#' subsample_table
#'
#' @description This function has been ported from MicrobeR. Takes a Feature/OTU/ASV table of counts where samples are columns and feature names are row names and returns a table with even coverage on a per-sample basis. Now this function is really just an alias for rarefy_even_depth from phyloseq to simplify its usage for non-phyloseq objects. NOTE: Sampling with replacement for single rarefraction method!
#'
#' @param features Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param depth Count depth, defaults to min(colSums(OTUTABLE)). Note, samples with less than the desired number of read will be removed!
#' @param seed A randomization SEED, defaults to 182. This is ignored for multiple subsamples and instead the seeds 1:NSAMPS is used.
#' @param with_replace Should subsampling be with or without replacement? (TRUE/FALSE) Default=FALSE to be analogous to qiime2; however, this may take longer to calculate. Please see phyloseq::rarefy_even_depth for more information.
#' @param verbose Should progress and metrics be printed to screen via message()? Default=TRUE
#' @return Subsampled Table
#' @export

subsample_table<-function(features, depth, seed, with_replace, verbose){
  if(missing(verbose)){verbose=T} 
  if(missing(depth)){depth=min(colSums(features))}
  if(missing(with_replace)){with_replace=FALSE}
  if(missing(seed)){seed=182}

  if(verbose==T){message(paste("Subsampling feature table to", depth, ", currently has ", nrow(features), " taxa."))}
  sub<-as.data.frame(rarefy_even_depth(otu_table(features, taxa_are_rows = T), sample.size=depth, rngseed=seed, replace=with_replace, verbose = FALSE))
  if(verbose==T){message(paste("...after subsampling there are", nrow(sub), "taxa with",round(sum(sub)/sum(features),4)*100, "% of reads retained from", ncol(sub),"of",ncol(features),"samples."))}
  return(sub)
}
