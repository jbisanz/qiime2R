#' plot the provenance of a QIIME2 artifact (.qza)
#'
#' extracts embedded data and object metadata into an R session
#'
#' @param artifact the object returned by read_qza
#' @return a text summary for the provenance information
#' @examples \dontrun{print_provenance(artifact)}
#' @export
#'
#'
#'


print_provenance<-function(artifact){
  if(missing(artifact)){stop("Artifact not provided...")}

  return(list.tree(artifact$provenance, maxcomp=1000, attr.print=FALSE))

}
