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
