#' A ggplot2 theme for publication-style figures
#'
#' @examples \dontrun{ggplot(data) + theme_q2r()}
#' @export
#'
#'
#'

theme_q2r<- function () { 
  theme_classic(base_size=8, base_family="Helvetica") +
  theme(panel.border = element_rect(color="black", size=1, fill=NA)) +
  theme(axis.line = element_blank(), strip.background = element_blank())
}