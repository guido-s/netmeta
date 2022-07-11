#' Calculate distance matrix for an adjacency matrix
#' 
#' @description
#' Calculate distance matrix for an adjacency matrix based on distance
#' algorithm by Müller et al. (1987).
#' 
#' @param x Either a netmeta object or an adjacency matrix.
#'
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}
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
#' data(Senn2013)
#' 
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD")
#' 
#' netdistance(net1)
#' netdistance(net1$A.matrix)
#' 
#' @export netdistance


netdistance <- function(x) {
  
  ## Calculate distance matrix D of adjacency matrix A based on
  ## distance algorithm by Mueller et al. (1987) using triangle
  ## inequality
  
  if (inherits(x, "netmeta"))
    A <- x$A.matrix
  else
    A <- x
  
  
  ## Check whether A is a matrix
  ##
  if (!is.matrix(A))
    stop("Argument 'x' must be a netmeta object or a matrix.")
  
  
  ## Starting value for D is sign(A), with 0 replaced by Inf
  ##
  n <- dim(A)[1] 
  D <- sign(A)
  ##
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (D[i, j] == 0) {
        D[i, j] <- Inf
        D[j, i] <- Inf
      }
    }
  }
  ##
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
  
  D
}
