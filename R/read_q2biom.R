#' read qiime2 biom file (version 2.1)
#'
#' Loads a version 2.1 spec biom file (http://biom-format.org/documentation/format_versions/biom-2.1.html) as expected to be found within a qiime2 artifact.
#'
#' @param file path to the input file, ex: file="~/Downloads/3372d9e0-3f1c-43d8-838b-35c7ad6dac89/data/feature-table.biom"

#' @return a matrix of values
#'
#' @examples \dontrun{metadata<-read_q2biom("feature-table.biom")}
#' @export
#'
#'


read_q2biom <- function(file) {
  if(missing(file)){stop("Path to biom file given")}
  if(!file.exists(file)){stop("File not found")}
  
  hdata<-h5read(file,"/")
  
  ftable<-
    sparseMatrix(
      p=hdata$observation$matrix$indptr,
      j=hdata$observation$matrix$indices,
      x=as.numeric(hdata$observation$matrix$data),
      index1=FALSE,
      dims=c(length(hdata$observation$ids), length(hdata$sample$ids)),
      dimnames=list(hdata$observation$ids,hdata$sample$ids)
    )
  
  return(as.matrix(ftable))
}
