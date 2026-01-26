#' Generic function for network graphs
#' 
#' @description
#' Generic function for network graphs
#' 
#' @param x An R object.
#' @param \dots Additional arguments.
#' 
#' @details
#' 
#' For more details, look at the following functions to generate
#' network graphs:
#' \itemize{
#' \item \code{\link{netgraph.netmeta}}
#' \item \code{\link{netgraph.netimpact}}
#' \item \code{\link{netgraph.netconnection}}
#' \item \code{\link{netgraph.netcomb}}
#' \item \code{\link{netgraph.discomb}}
#' }
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de
#' }
#' 
#' @keywords hplot
#'
#' @examples
#' # Examples: example(netgraph.netmeta)
#' 
#' @rdname netgraph
#' @export netgraph

netgraph <- function(x, ...)
  UseMethod("netgraph")
