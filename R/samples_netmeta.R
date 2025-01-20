samples_netmeta <- function(x, nsim, pooled) {
  
  chkclass(x, "netmeta")
  #
  pooled <- setchar(pooled, c("common", "random"))
  #
  if (pooled == "common") {
    TE.pooled <- x$TE.common
    Cov.pooled <- x$Cov.common
  }
  else {
    TE.pooled <- x$TE.random
    Cov.pooled <- x$Cov.random
  }
  #
  theta <- TE.pooled[, 1]
  
  compMatrix <- matrix(0, nrow = nrow(Cov.pooled), ncol = length(x$trts))
  rownames(compMatrix) <- rownames(Cov.pooled)
  colnames(compMatrix) <- x$trts
  #
  allcomps <- compsplit(rownames(compMatrix), x$sep.trts)
  for (i in seq_len(nrow(Cov.pooled)))
    compMatrix[i, allcomps[[i]]] <- 1
  #
  var.theta <- as.vector(ginv(compMatrix) %*% diag(Cov.pooled))
  
  samples <- rmvnorm(nsim, theta, diag(var.theta))
  rownames(samples) <- seq_len(nrow(samples))
  colnames(samples) <- x$trts
  
  res <- list(samples = samples, compMatrix = compMatrix)
  #
  res
}
