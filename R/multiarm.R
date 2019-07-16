multiarm <- function(r) {
  ##
  ## Dimension of r and R
  ##
  m <- length(r) # Number of edges
  ##
  k <- (1 + sqrt(8 * m + 1)) / 2 # Number of vertices
  if (!(abs(k - round(k)) < .Machine$double.eps^0.5))
    stop("Wrong number of comparisons in multi-arm study.", call. = FALSE)
  ##
  ## Construct edge-vertex incidence matrix of complete graph of
  ## dimension k
  ##
  B <- createB(ncol = k)
  ##
  ## Distribute the edge variances on a symmetrical matrix R of
  ## dimension k x k
  ##
  R <- diag(diag(t(B) %*% diag(r, nrow = m) %*% B)) -
    t(B) %*% diag(r, nrow = m) %*% B
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
  W[W < 0 & abs(W) < .Machine$double.eps^0.25] <- 0
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
