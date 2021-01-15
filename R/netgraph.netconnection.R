#' Network graph for objects of class netconnection
#' 
#' @description
#' This function generates a graph of the evidence network.
#' 
#' @param x An object of class \code{netconnection}.
#' @param \dots Additional arguments passed on to
#'   \code{\link{netgraph.netmeta}} (see Details).
#' 
#' @details
#' The following arguments are used internally and cannot be specified
#' by the user: \code{thickness}, \code{seq}, \code{iterate}.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de},
#'   Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netconnection}}, \code{\link{netgraph.netmeta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' # Artificial example with two subnetworks
#' #
#' t1 <- c("G", "B", "B", "D", "A", "F")
#' t2 <- c("B", "C", "E", "E", "H", "A")
#' #
#' nc1 <- netconnection(t1, t2)
#' print(nc1, details = TRUE)
#' 
#' netgraph(nc1, plastic = FALSE)
#' 
#' @method netgraph netconnection
#' @export
#' @export netgraph.netconnection


netgraph.netconnection <- function(x, ...) {
  
  
  meta:::chkclass(x, "netconnection")
  
  x$seq <- rownames(x$A.matrix)
  x$d <- 2
  x$trts <- rownames(x$A.matrix)
  ##
  class(x) <- "netmeta"


  res <- netgraph(x, thickness = "equal", seq = x$seq, iterate = FALSE, ...)

  
  invisible(res)
}
