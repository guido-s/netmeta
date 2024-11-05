#' Print and change default network meta-analysis settings in R package
#' \bold{netmeta}
#' 
#' @description
#' Print and change default settings to conduct and print or plot
#' network meta-analyses in R package \bold{netmeta}.
#' 
#' @param ... Arguments to change default settings.
#' @param quietly A logical indicating whether information on settings
#'   should be printed.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link[meta]{settings.meta}}, \code{\link[meta]{gs}}
#' 
#' @export settings.netmeta

settings.netmeta <- function(..., quietly = TRUE) {
  
  #
  # Check argument
  #
  missing.quietly <- missing(quietly)
  chklogical(quietly)
  
  netargs <- gs(".argslist.netmeta")

  
  settings.meta(..., quietly = TRUE)
}
