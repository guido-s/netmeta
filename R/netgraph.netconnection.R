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
#'   connecting treatments (edges) if argument \code{plastic = FALSE}
#'   (see Details).
#' @param reference.group Reference treatment (only relevant for
#'   disconnected networks).
#' @param plastic A logical indicating whether the appearance of the
#'   comparisons should be in '3D look'.
#' @param \dots Additional arguments passed on to
#'   \code{\link{netgraph.netmeta}} (see Details).
#' 
#' @details
#' Argument \code{col} can be a single color for all edges, a vector
#' of length equal to the number of edges, or a vector of length equal
#' to the number of subnetworks. Argument \code{reference.group} is
#' only considered in disconnected networks, i.e., if more than one
#' (sub)network exists, and if argument \code{col} provides colors for
#' subnetworks. In this case, the first color provided in argument
#' \code{col} defines the color for the subnetwork with the reference
#' treatment.
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
#' netgraph(nc1, points = TRUE, adj = 0.5, bg.points = "lightgray")
#' netgraph(nc1, points = TRUE, adj = 0.5, bg.points = "lightgray",
#'   plastic = TRUE)
#' 
#' @method netgraph netconnection
#' @export


netgraph.netconnection <- function(x, seq,
                                   col = seq_len(x$n.subnets),
                                   reference.group = NULL,
                                   plastic = FALSE,
                                   ...) {
  
  
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
  
  
  ## Check length of argument 'col'
  ##
  col.subnets <- FALSE
  ##
  if (length(col) == 1) {
    col.subnets <- TRUE
    col <- rep_len(col, x$n.subnets)
  }
  else if (length(col) == x$n.subnets)
    col.subnets <- TRUE
  else
    chklength(col, length(x$subnet.comparisons),
              text = paste0("Length of argument 'col' must be equal to ",
                            "number of edges (",
                            length(x$subnet.comparisons), ")",
                            if (x$n.subnets > 1)
                              paste0(" or subnetworks (",
                                     x$n.subnets, ").")
                            else "."))
  
  
  ## Reorder subnetworks if reference treatment is not part of first
  ## subnetwork
  ##
  if (!is.null(reference.group)) {
    if (!col.subnets)
      warning("Argument 'reference.group' not considered as",
              "a color is defined for each edge.",
              call. = FALSE)
    else {
      reference.group <- setref(reference.group, x$trts)
      ##
      sel.ref <-
        x$treat1 == reference.group | x$treat2 == reference.group
      sel.ref <- unique(x$subnet[sel.ref])
      ##
      if (sel.ref != 1)
        x$subnet.comparisons <-
          as.numeric(relevel(as.factor(x$subnet.comparisons),
                             ref = sel.ref))
    }
  }
  
  
  ## Set color for edges in (sub)networks
  ##
  if (col.subnets) {
    col <-
      as.character(factor(x$subnet.comparisons,
                          levels = seq_len(x$n.subnet),
                          labels = col))
  }
  
  
  res <- netgraph(x, plastic = plastic, col = col,
                  seq = seq, ...)

  
  invisible(res)
}
