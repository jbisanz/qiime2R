#' Display an interactive table
#'
#' Uses DT:datatable to embed an interactive table with export options in an R session or markdown document
#'
#' @param table a data.frame/tibble/matrix
#' @param nrow number of rows to display in table (default=10)
#'
#' @examples \dontrun{interactive_table(metadata)}
#' @export
#'
#'

interactive_table<-function(table, nrow){
  if(missing(nrow))(nrow=10)
  dtable<-datatable(table, extensions='Buttons', filter="top", options=list(pageLength=nrow, dom='Bfrtip', buttons=c('copy','csv','excel', 'pdf') ))
  return(dtable)
}