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
  
  defline<-suppressWarnings(readLines(file)[2])
  if(!grepl("^#q2:types", defline)){stop("Metadata does not define types (ie second line does not start with #q2:types)")}
  
  defline<-strsplit(defline, split="\t")[[1]]
  
  defline[grep("numeric", defline)]<-"as.numeric"
  defline[grep("categorical|q2:types", defline)]<-"as.factor"
  
  metadata<-subset(read.table(file, header=F, comment.char = "#", sep='\t'), !V1 %in% c("#id","id","sampleid","sample id","sample-id","#SampleID","#Sample ID", "sample_name", "SampleID","Sample ID"))
  colnames(metadata)<-strsplit(suppressWarnings(readLines(file)[1]), "\t")[[1]]
  colnames(metadata)[1]<-"SampleID"
  return(metadata)
}
