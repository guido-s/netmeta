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
#' @keywords hplot
#'
#' @examples
#' data(Senn2013)
#' 
#' # Generation of an object of class 'netmeta' with reference
#' # treatment 'plac'
#' #
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'                 data = Senn2013, sm = "MD", reference = "plac")
#' 
#' # Network graph with default settings
#' #
#' netgraph(net1)
#' 
#' @rdname netgraph
#' @export netgraph


netgraph <- function(x, ...)
  UseMethod("netgraph")

