samples_netmeta <- function(x, nsim, pooled) {
  
  chkclass(x, "netmeta")
  #
  pooled <- setchar(pooled, c("common", "random"))
  #
  if (pooled == "common") {
    TE.pooled <- x$TE.common
    seTE.pooled <- x$seTE.common
  }
  else {
    TE.pooled <- x$TE.random
    seTE.pooled <- x$seTE.random
  }
  #
  theta <- TE.pooled[, 1]
  
  # Comparison matrix with all pairwise comparisons in rows and treatments
  # in columns
  #
  pairs <- combn(x$trts, 2)
  n.pairs <- ncol(pairs)
  #
  compMatrix <- matrix(0, nrow = n.pairs, ncol = length(x$trts),
    dimnames = list(apply(pairs, 2, paste, collapse = x$sep.trts), x$trts))
  #
  for (i in seq_len(n.pairs))
    compMatrix[i, pairs[, i]] <- 1
  
  # Covariance of treatment effects relative to the first treatment
  #
  V <- seTE.pooled^2
  #
  Sigma <- (outer(V[, 1], V[, 1], "+") - V) / 2
  
  # Generate samples
  #
  samples <- rmvnorm(nsim, theta, Sigma)
  rownames(samples) <- seq_len(nrow(samples))
  colnames(samples) <- x$trts
  
  # Return samples
  #
  res <- list(samples = samples, compMatrix = compMatrix)
  res
}
