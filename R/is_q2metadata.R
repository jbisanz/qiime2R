#' checks if metadata is in qiime2 (.tsv)
#'
#' Checks to see if a file is in qiime2 metadata format, ie contains #q2:types line dictating the type of variable (categorical/numeric)
#'
#' @param file path to the input file, ex: file="~/data/moving_pictures/table.qza"

#' @return TRUE/FALSE
#'
#' @examples \dontrun{metadata<-is_q2metadata("q2metadata.tsv")}
#' @export
#'
#'

is_q2metadata <- function(file){
  suppressWarnings(
  if(grepl("^#q2:types", readLines(file)[2])){return(TRUE)}else{return(FALSE)}
  )
}
