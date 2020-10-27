#' read qiime2 metadata (.tsv)
#'
#' Loads a qiime2 metadata file wherein the 2nd line contains the #q2:types line dictating the type of variable (categorical/numeric)
#'
#' @param file path to the input file, ex: file="~/data/moving_pictures/table.qza"

#' @return a data.frame wherein the first column is SampleID
#'
#' @examples \dontrun{metadata<-read_q2metadata("q2metadata.tsv")}
#' @export
#'
#'


read_q2metadata <- function(file) {
  if(missing(file)){stop("Path to metadata file not found")}
  if(!is_q2metadata(file)){stop("Metadata does not define types (ie second line does not start with #q2:types)")}
  
  defline<-suppressWarnings(readLines(file)[2])
  defline<-strsplit(defline, split="\t")[[1]]
  
  defline[grep("numeric", tolower(defline))]<-"double"
  defline[grep("categorical|q2:types", tolower(defline))]<-"factor"
  defline[defline==""]<-"factor"
  
  coltitles<-strsplit(suppressWarnings(readLines(file)[1]), split='\t')[[1]]
  metadata<-read.table(file, header=F, col.names=coltitles, skip=2, sep='\t', colClasses = defline, check.names = FALSE)
  colnames(metadata)[1]<-"SampleID"
  
  return(metadata)
}
