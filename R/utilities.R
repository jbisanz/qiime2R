#' Utility Functions
#'
#' Perform misc utility functions:
#' @param corner() show the the top left corner of a data.frame/tibble/matrix.
#' @param min_nonzero() Find the smallest non-zero, non-NA value in a vector.
#' 
#' @export

corner<-function(table, nrow, ncol){
  
  if(missing(table)){stop("Table-like object not provided")}
  if(missing(nrow)){nrow=5}
  if(missing(ncol)){ncol=5}
  
  return(table[1:nrow, 1:ncol])
}

min_nonzero<-function(data){
  if(missing(data) | !is.vector(data) | !is.numeric(data)){stop("Vector not provided, or is not a vector, or contains non-numeric entries.")}
  data<-data[!is.na(data)]
  data<-data[data!=0]
  return(min(data))
}

mean_sd=function(x){data.frame(y=mean(x), ymin=mean(x)-sd(x), ymax=mean(x)+sd(x))}
