#' Network graph for objects of class netcomb
#' 
#' @description
#' This function generates a graph of the evidence network.
#' 
#' @param x An object of class \code{netcomb}.
#' @param \dots Additional arguments passed on to
#'   \code{\link{netgraph.netmeta}}.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de},
#'   Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netcomb}}, \code{\link{netgraph.netmeta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Linde2016)
#' 
#' # Only consider studies including Face-to-face PST (to reduce
#' # runtime of example)
#' #
#' face <- subset(Linde2016, id %in% c(16, 24, 49, 118))
#' 
#' # Conduct random effects network meta-analysis
#' #
#' net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'                 data = face, ref = "placebo",
#'                 sm = "OR", fixed = FALSE)
#' 
#' # Additive model for treatment components (with placebo as inactive
#' # treatment)
#' #
#' nc1 <- netcomb(net1, inactive = "placebo")
#'
#' netgraph(nc1)
#' 
#' @method netgraph netcomb
#' @export


netgraph.netcomb <- function(x, ...) {
  
  
  chkclass(x, "netcomb")
  
  
  res <- netgraph(x$x, ...)
  
  
  invisible(res)
}
