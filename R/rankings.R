rankings <- function(x) {
  
  #
  # Calculate ranks and rankogram
  #
  
  nsim <- nrow(x)
  n.trts <- ncol(x)
  #
  rnk <- apply(x, 1, rank, ties.method = "random")
  #
  tab <- apply(rnk, 1, table)
  #
  if (is.list(tab)) {
    rankogram <- matrix(0, nrow = n.trts, ncol = n.trts)
    rownames(rankogram) <- names(tab)
    colnames(rankogram) <- seq_len(n.trts)
    ##
    for (i in names(tab))
      rankogram[i, names(tab[[i]])] <- tab[[i]][names(tab[[i]])]
  }
  else
    rankogram <- t(as.data.frame(tab))
  
  #
  # Cumulative ranks and SUCRAs
  #
  
  cumrank <- t(apply(rankogram, 1, cumsum))
  sucras <- apply(cumrank[, -n.trts], 1, sum) / (n.trts - 1) / nsim
  
  #
  # Return results
  #
  
  res <- list(sucras = sucras,
              rankogram = rankogram / nsim,
              cumrank = cumrank / nsim,
              nsim = nsim)
  #
  res
}
