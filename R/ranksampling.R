ranksampling <- function(x, nsim,
                         pooled = "random", small.values = "good",
                         keep.sample = FALSE) {
  chkclass(x, "netmeta")
  pooled <- setchar(pooled, c("common", "random"))
  small.values <- setchar(small.values, c("good", "bad"))
  chklogical(keep.sample)
  ##
  if (pooled == "common") {
    TE.pooled <- x$TE.common
    Cov.pooled <- x$Cov.common
  }
  else {
    TE.pooled <- x$TE.random
    Cov.pooled <- x$Cov.random
  }
  ##
  if (small.values == "good")
    theta <- TE.pooled[, 1]
  else
    theta <- TE.pooled[1, ]
  ##
  compMatrix <- matrix(0, nrow = nrow(Cov.pooled), ncol = length(x$trts))
  rownames(compMatrix) <- rownames(Cov.pooled)
  colnames(compMatrix) <- x$trts
  ##
  allcomps <- compsplit(rownames(compMatrix), x$sep.trts)
  for (i in seq_len(nrow(Cov.pooled)))
    compMatrix[i, allcomps[[i]]] <- 1
  ##
  var.theta <- as.vector(ginv(compMatrix) %*% diag(Cov.pooled))
  ##
  sample <- mvtnorm::rmvnorm(nsim, theta, diag(var.theta))
  rownames(sample) <- seq_len(nrow(sample))
  colnames(sample) <- x$trts
  ##
  ## Ranks
  ##
  rnk <- apply(sample, 1, rank, ties.method = "random")
  ##
  ## Rankogram
  ##
  tab <- apply(rnk, 1, table)
  ##
  if (is.list(tab)) {
    rankogram <- matrix(0, nrow = x$n, ncol = x$n)
    rownames(rankogram) <- names(tab)
    colnames(rankogram) <- seq_len(x$n)
    ##
    for (i in names(tab))
      rankogram[i, names(tab[[i]])] <- tab[[i]][names(tab[[i]])]
  }
  else
    rankogram <- t(as.data.frame(tab))
  ##
  ## Cumulative ranks
  ##
  cumrank <- t(apply(rankogram, 1, cumsum))
  ##
  ## SUCRAs
  ##
  ranking <- apply(cumrank[, -x$n], 1, sum) / (x$n - 1)
  ##
  ## Return results
  ##
  res <- list(ranking = ranking / nsim,
              rankogram = rankogram / nsim,
              cumrank = cumrank / nsim,
              ##
              nsim = nsim,
              pooled = pooled,
              small.values = small.values,
              keep.sample = keep.sample,
              ##
              compMatrix = compMatrix)
  ##
  if (keep.sample)
    res[["sample"]] <- sample
  ##
  res
}
