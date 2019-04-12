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
#' }
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de
#' }
#' 
#' @seealso \code{\link{netgraph.netmeta}},
#'   \code{\link{netgraph.netimpact}}
#' 
#' @rdname netgraph
#' @export netgraph


netgraph <- function(x, ...)
  UseMethod("netgraph")

