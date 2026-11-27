rankings <- function(x) {
  
  #
  # Calculate ranks and rankogram
  #
  
  nsim <- nrow(x)
  n.trts <- ncol(x)
  #
  rnk <- apply(x, 1, rank, ties.method = "random")
  #
  rankogram <- matrix(0, nrow = n.trts, ncol = n.trts)
  rownames(rankogram) <- colnames(x)
  colnames(rankogram) <- seq_len(n.trts)
  #
  for (i in seq_len(n.trts))
    rankogram[i, ] <- tabulate(rnk[i, ], nbins = n.trts)
  #
  rankogram <- rankogram / nsim
  
  #
  # Cumulative ranks and SUCRAs
  #
  
  cumrank <- t(apply(rankogram, 1, cumsum))
  sucras <- rowSums(cumrank[, -n.trts, drop = FALSE]) / (n.trts - 1)
  
  #
  # Mean and median ranks
  #
  
  meanranks <- rowSums(col(rankogram) * rankogram)
  #
  getmed <- function(x) which(x >= 0.5)[1]
  medianranks <- apply(cumrank, 1, getmed)
  
  #
  # Return results
  #
  
  res <- list(sucras = sucras,
              rankogram = rankogram,
              cumrank = cumrank,
              #
              meanranks = meanranks,
              medianranks = medianranks,
              #
              nsim = nsim)
  #
  res
}
