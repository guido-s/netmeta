#' Network graph for objects of class discomb
#' 
#' @description
#' This function generates a graph of the evidence network.
#' 
#' @param x An object of class \code{discomb}.
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
#' @seealso \code{\link{discomb}}, \code{\link{netgraph.netmeta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' # Artificial dataset
#' #
#' t1 <- c("A + B", "A + C", "A"    , "A"    , "D", "D", "E")
#' t2 <- c("C"    , "B"    , "B + C", "A + D", "E", "F", "F")
#' #
#' mean <- c(4.1, 2.05, 0, 0, 0.1, 0.1, 0.05)
#' se.mean <- rep(0.1, 7)
#' #
#' study <- paste("study", c(1:4, 5, 5, 5))
#' #
#' dat <- data.frame(mean, se.mean, t1, t2, study,
#'                   stringsAsFactors = FALSE)
#' #
#' trts <- c("A", "A + B", "A + C", "A + D",
#'   "B", "B + C", "C", "D", "E", "F")
#' #
#' comps <- LETTERS[1:6]
#' 
#' # Use netconnection() to display network information
#' #
#' netconnection(t1, t2, study)
#' 
#' dc1 <- discomb(mean, se.mean, t1, t2, study, seq = trts)
#'
#' netgraph(dc1)
#' 
#' @method netgraph discomb
#' @export


netgraph.discomb <- function(x, ...) {
  
  
  chkclass(x, "discomb")
  
  
  y <- netconnection(x$treat1, x$treat2, x$studlab)
  ##
  y$seq <- rownames(y$A.matrix)
  y$d <- 2
  y$trts <- rownames(y$A.matrix)
  ##
  class(y) <- "netmeta"
  
  
  res <- netgraph(y, thickness = "equal", seq = y$seq, iterate = FALSE, ...)
  
  
  invisible(res)
}
