#' plot the provenance of a QIIME2 artifact (.qza)
#'
#' extracts embedded data and object metadata into an R session
#'
#' @param artifact the object returned by read_qza
#' @return a ggplot object which can be printed to the screen and/or modified
#' @export
#'
#'
#'

#Not currently exported

plot_provenance<-function(artifact){
  if(missing(artifact)){stop("Artifact not provided...")}

  prv<-artifact$provenance

  lapply(prv, function(x){
    ins<-x$action$inputs
  })


  plot<-1
  return(plot)
}
