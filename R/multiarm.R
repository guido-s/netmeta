multiarm <- function(r) {
  ##
  ## Dimension of r and R
  ##
  m <- length(r)                 # Number of edges
  k <- (1 + sqrt(8 * m + 1)) / 2 # Number of vertices
  ##
  ## Construct edge.vertex incidence matrix of
  ## complete graph of dimension k
  ##
  B <- createB(ncol = k)
  ##
  ## Distribute the edge variances on a symmetrical k x k matrix, R
  ##
  R <- diag(diag(t(B) %*% diag(r) %*% B)) - t(B) %*% diag(r) %*% B
  ##
  ## Construct pseudoinverse Lt from given variance (resistance) matrix R
  ## using a theorem equivalent to Theorem 7 by Gutman & Xiao
  ## Lt <- -0.5 * (R - (R %*% J + J %*% R) / k + J %*%R %*% J / k^2)
  ##
  Lt <- -0.5 * t(B) %*% B %*% R %*% t(B) %*% B / k^2
  ##
  ## Compute Laplacian matrix L from Lt
  ## 
  L <- ginv(Lt)
  ##
  ## Compute weight matrix W and variance matrix V from Laplacian L
  ## 
  W <- diag(diag(L)) - L
  W[W < 0 & abs(W) < .Machine$double.eps^0.75] <- 0
  ##
  V <- 1 / W
  ##
  ## Compute original variance vector v from V
  ##
  v <- rep(0, m)
  edge <- 0
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      edge <- edge + 1
      v[edge] <- V[i, j]
    }
  }
  ##
  ## Result
  ##
  res <- list(k = k, r = r, R = R, Lt = Lt, L = L, W = W, V = V, v = v)
  res
}
