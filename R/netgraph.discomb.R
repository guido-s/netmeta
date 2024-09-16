#' Network graph for objects of class discomb
#' 
#' @description
#' This function generates a graph of the evidence network.
#' 
#' @param x An object of class \code{discomb}.
#' @param labels An optional vector with treatment labels.
#' @param adj One, two, or three values in [0, 1] (or a vector /
#'   matrix with length / number of rows equal to the number of
#'   treatments) specifying the x (and optionally y and z) adjustment
#'   for treatment labels.
#' @param offset Distance between edges (i.e. treatments) in graph and
#'   treatment labels for 2-D plots (value of 0.0175 corresponds to a
#'   difference of 1.75\% of the range on x- and y-axis).
#' @param rotate A single numeric with value between -180 and 180
#'   specifying the angle to rotate nodes in a circular network.
#' @param points A logical indicating whether points should be printed
#'   at nodes (i.e. treatments) of the network graph.
#' @param cex.points Corresponding point size. Can be a vector with length
#'   equal to the number of treatments.
#' @param \dots Additional arguments passed on to
#'   \code{\link{netgraph.netmeta}}.
#' 
#' @details
#' The arguments \code{seq} and \code{iterate} are used internally and
#' cannot be specified by the user.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de},
#'   Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}
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


netgraph.discomb <- function(x,
                             labels = x$trts,
                             adj = NULL,
                             offset =
                               if (!is.null(adj) && all(unique(adj) == 0.5))
                                 0
                             else
                               0.0175,
                             rotate = 0,
                             points = !missing(cex.points),
                             cex.points = 1, ...) {
  
  
  chkclass(x, "discomb")
  
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  if (!missing(labels))
    labels <- catch("labels", mc, x, sfsp)
  if (!missing(adj))
    adj <- catch("adj", mc, x, sfsp)
  if (!missing(offset))
    offset <- catch("offset", mc, x, sfsp)
  if (!missing(rotate))
    rotate <- catch("rotate", mc, x, sfsp)
  if (!missing(cex.points))
    cex.points <- catch("cex.points", mc, x, sfsp)
  
  y <- netconnection(x$treat1, x$treat2, x$studlab)
  #
  y$trts <- x$trts
  y$D.matrix <- y$D.matrix[x$trts, x$trts]
  y$A.matrix <- y$A.matrix[x$trts, x$trts]
  y$L.matrix <- y$L.matrix[x$trts, x$trts]
  #
  class(y) <- "netmeta"
  
  res <- netgraph(y, seq = y$seq, iterate = FALSE,# plastic = plastic,
                  #
                  labels = labels,
                  adj = adj,
                  offset = offset,
                  rotate = rotate,
                  points = points,
                  cex.points = cex.points,
                  #
                  ...)
  
  invisible(res)
}
