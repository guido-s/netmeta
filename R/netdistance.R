#' Calculate distance matrix for an adjacency matrix
#' 
#' @description
#' Calculate distance matrix for an adjacency matrix based on distance
#' algorithm by Müller et al. (1987).
#' 
#' @aliases netdistance netdistance.default netdistance.netmeta
#'   netdistance.netcomb print.netdistance
#' 
#' @param x Either a netmeta or netcomb object or an adjacency matrix.
#' @param lab.Inf A character string to label infinite values.
#' @param \dots Additional arguments (ignored).
#'
#' @author Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netconnection}}
#' 
#' @references
#' Müller WR, Szymanski K, Knop JV, and Trinajstic N (1987):
#' An algorithm for construction of the molecular distance matrix.
#' \emph{Journal of Computational Chemistry},
#' \bold{8}, 170--73
#' 
#' @examples
#' data(smokingcessation)
#' 
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' net1 <- netmeta(p1, common = FALSE)
#' 
#' netdistance(net1)
#' 
#' \dontrun{
#' data(Senn2013)
#' 
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD")
#' 
#' netdistance(net1)
#' netdistance(net1$A.matrix)
#' }
#' 
#' @rdname netdistance
#' @method netdistance default
#' @export

netdistance.default <- function(x) {
  
  # Calculate distance matrix D of adjacency matrix A based on
  # distance algorithm by Mueller et al. (1987) using triangle
  # inequality
  
  chkclass(x, "matrix")
  #
  A <- x
  
  # Starting value for D is sign(A), with 0 replaced by Inf
  #
  n <- nrow(A)
  D <- sign(A)
  #
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (D[i, j] == 0) {
        D[i, j] <- Inf
        D[j, i] <- Inf
      }
    }
  }
  #
  for (d in 1:(n - 1)) {
    for (i in 1:n) {
      for (j in 1:n) {
        if (D[i, j] == d) {
          for (k in 1:n) {
            akj <- D[k, i] + d # = D[k, i] + D[i, j]
            D[k, j] <- min(D[k, j], akj)
          }
        }
      }
    }
  }
  #
  maxdist <- nrow(D)
  D2 <- D
  D2[is.infinite(D2)] <- maxdist
  attr(D, "order") <- hclust(dist(D2))$order
  
  class(D) <- c("netdistance", class(D))
  #
  D
}





#' @rdname netdistance
#' @method netdistance netmeta
#' @export

netdistance.netmeta <- function(x) {
  
  chkclass(x, "netmeta")

  A <- x$A.matrix
  seq <- netconnection(x$treat1, x$treat2)$seq
  A <- A[seq, seq]
  
  res <- netdistance(A)
  attr(res, "order") <- NULL
  #
  res
}





#' @rdname netdistance
#' @method netdistance netcomb
#' @export

netdistance.netcomb <- function(x) {
  
  chkclass(x, "netcomb")
  
  if (inherits(x, "discomb")) {
    A <- x$A.matrix
    seq <- netconnection(x$treat1, x$treat2)$seq
  }
  else {
    A <- x$x$A.matrix
    seq <- netconnection(x$x$treat1, x$x$treat2)$seq
  }
  #
  A <- A[seq, seq]
  
  res <- netdistance(A)
  attr(res, "order") <- NULL
  #
  res
}





#' @rdname netdistance
#' @method netdistance netconnection
#' @export

netdistance.netconnection <- function(x) {
  
  chkclass(x, "netconnection")
    
  netdistance(x$A.matrix)
}





#' @rdname netdistance
#' @method print netdistance
#' @export

print.netdistance <- function(x, lab.Inf = ".", ...) {
  o <- attr(x, "order")
  #
  if (!is.null(o))
    x <- x[o, o]
  #
  x[is.infinite(x)] <- lab.Inf
  #
  prmatrix(x, quote = FALSE, right = TRUE)
  #
  invisible(NULL)
}





#' @rdname netdistance
#' @export


netdistance <- function(x)
  UseMethod("netdistance")
