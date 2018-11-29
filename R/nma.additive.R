nma.additive <- function(TE, weights, studlab,
                         treat1, treat2,
                         level.comb,
                         X, C.matrix, B.matrix,
                         Q, df.Q.additive, df.Q.diff,
                         n, sep.trts) {
  
  
  m <- length(TE)
  ##
  ## Adjusted weights
  ##
  W <- diag(weights,
            nrow = length(weights))
  ##
  ## Laplacian matrix and pseudoinverse of L
  ##
  L <- t(X) %*% W %*% X
  Lplus <- ginv(L) # = Cov matrix of beta (components)
  colnames(Lplus) <- colnames(L)
  rownames(Lplus) <- rownames(L)
  ##
  ## H matrix
  ##
  H <- X %*% Lplus %*% t(X) %*% W
  
  
  ##
  ## beta = effects of components
  ##
  beta <- as.vector(Lplus %*% t(X) %*% W %*% TE)
  se.beta <- sqrt(diag(Lplus))
  names(beta) <- names(se.beta)
  ##
  ## theta = estimates for combination treatments
  ##
  theta <- as.vector(C.matrix %*% beta)
  se.theta <- sqrt(diag(C.matrix %*% Lplus %*% t(C.matrix)))
  names(theta) <- names(se.theta)
  ##
  ## delta = treatment estimates for observed comparisons
  ##
  delta <- as.vector(X %*% beta) # = B.matrix %*% theta = H %*% TE
  se.delta <- unname(sqrt(diag(X %*% Lplus %*% t(X))))
  ##
  ## delta.all = all direct and indirect treatment estimates
  ##
  B.full <- createB(ncol = n)
  X.full <- B.full %*% C.matrix
  colnames(X.full) <- colnames(C.matrix)
  ##
  labels <- colnames(B.matrix)
  ##
  k <- 0
  lab <- vector(mode = "numeric", length = choose(n, 2))
  ##
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      k <- k + 1
      lab[k] <- paste(labels[i], labels[j], sep = sep.trts)
    }
  }
  ##
  rownames(X.full) <- lab
  ##
  delta.full <- as.vector(X.full %*% beta)
  se.delta.full <- sqrt(diag(X.full %*% Lplus %*% t(X.full)))
  names(delta.full) <- names(se.delta.full)
  ##
  delta.all <- se.delta.all <- matrix(0, ncol = n, nrow = n)
  ##
  k <- 0
  ##
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      k <- k + 1
      delta.all[i, j] <-  delta.full[k]
      delta.all[j, i] <- -delta.full[k]
      se.delta.all[i, j] <- se.delta.all[j, i] <- se.delta.full[k]
    }
  }
  ##
  colnames(delta.all) <- rownames(delta.all) <-
    colnames(se.delta.all) <- rownames(se.delta.all) <- labels
  
  
  comparisons <- c(list(studlab = studlab, treat1 = treat1, treat2 = treat2),
                   meta::ci(delta, se.delta, level = level.comb))
  ##
  all.comparisons <- meta::ci(delta.all, se.delta.all, level = level.comb)
  ##
  components <- meta::ci(beta, se.beta, level = level.comb)
  ##
  combinations <- meta::ci(theta, se.theta, level = level.comb)
  ##
  ## Test of total heterogeneity / inconsistency:
  ##
  Q.additive <- as.vector(t(delta - TE) %*% W %*% (delta - TE))
  ##
  if (is.na(df.Q.additive) | df.Q.additive == 0)
    pval.Q.additive <- NA
  else
    pval.Q.additive <- 1 - pchisq(Q.additive, df.Q.additive)
  ##
  ## Difference to standard network meta-analysis model
  ##
  Q.diff <- Q.additive - Q
  if (!is.na(Q.diff) && abs(Q.diff) < .Machine$double.eps^0.75)
    Q.diff <- 0
  ##
  if (is.na(df.Q.diff) | df.Q.diff == 0)
    pval.Q.diff <- NA
  else
    pval.Q.diff <- 1 - pchisq(Q.diff, df.Q.diff)
  ##
  ## Heterogeneity variance
  ##
  I <- diag(m)
  E <- matrix(0, nrow = m, ncol = m)
  for (i in 1:m)
    for (j in 1:m)
      E[i, j] <- as.numeric(studlab[i] ==  studlab[j])
  ##
  if (df.Q.additive == 0) {
    tau2 <- NA
    tau <- NA
    I2 <- NA
  }
  else {
    tau2 <- max(0, (Q.additive - df.Q.additive) /
                   sum(diag((I - H) %*%
                            (B.matrix %*% t(B.matrix) * E / 2) %*% W)))
    tau <- sqrt(tau2)
    I2 <- meta:::isquared(Q.additive, df.Q.additive, 0.95)$TE
  }
  
  
  res <- list(comparisons = comparisons,
              all.comparisons = all.comparisons,
              components = components,
              combinations = combinations,
              ##
              Q.additive = Q.additive,
              df.Q.additive = df.Q.additive,
              pval.Q.additive = pval.Q.additive,
              ##
              Q.diff = Q.diff,
              df.Q.diff = df.Q.diff,
              pval.Q.diff = pval.Q.diff,
              ##
              tau = tau,
              I2 = I2)
  
  
  res
}
