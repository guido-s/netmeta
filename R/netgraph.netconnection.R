#' Network graph for objects of class netconnection
#' 
#' @description
#' This function generates a graph of the evidence network.
#' 
#' @param x An object of class \code{netconnection}.
#' @param seq A character or numerical vector specifying the sequence
#'   of treatments arrangement (anticlockwise if \code{start.layout =
#'   "circle"}).
#' @param col A single color (or vector of colors) for lines
#'   connecting treatments (edges) if argument \code{plastic = FALSE}.
#' @param plastic A logical indicating whether the appearance of the
#'   comparisons should be in '3D look'.
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


netgraph.netconnection <- function(x, seq,
                                   col = x$subnet.comparisons,
                                   plastic = FALSE, ...) {
  
  
  chkclass(x, "netconnection")
  
  if (missing(seq))
    seq <- rownames(x$A.matrix)
  else
    seq <- setseq(seq, rownames(x$A.matrix))
  ##
  x$d <- 2
  x$trts <- rownames(x$A.matrix)
  ##
  class(x) <- "netmeta"


  res <- netgraph(x, plastic = plastic, thickness = "equal",
                  seq = seq, iterate = FALSE,
                  col = col, ...)

  
  invisible(res)
}
