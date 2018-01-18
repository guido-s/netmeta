nma.additive <- function(TE, weights, studlab,
                         treat1, treat2, level,
                         X, C.matrix,
                         Q, df.Q.comp, df.Q.diff) {
  
  
  ##
  ## Adjusted weights
  ##
  W <- diag(weights)
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
  beta <- Lplus %*% t(X) %*% W %*% TE
  se.beta <- sqrt(diag(Lplus))
  ##
  ## theta = estimates for combination treatments
  ##
  theta <- C.matrix %*% beta
  se.theta <- sqrt(diag(C.matrix %*% Lplus %*% t(C.matrix)))
  ##
  ## delta = estimates for observed comparisons
  ##
  delta <- X %*% beta # = B.matrix %*% theta = H %*% TE
  se.delta <- sqrt(diag(X %*% Lplus %*% t(X)))
  
  
  comparisons <- c(list(studlab = studlab, treat1 = treat1, treat2 = treat2),
                   meta::ci(delta, se.delta, level = level))
  ##
  components <- meta::ci(beta, se.beta, level = level)
  ##
  combinations <- meta::ci(theta, se.theta, level = level)
  ##
  ## Model fit
  ##
  Q.comp <- as.vector(t(delta - TE) %*% W %*% (delta - TE))
  pval.Q.comp <- 1 - pchisq(Q.comp, df.Q.comp)
  ##
  ## Difference to standard network meta-analysis model
  ##
  Q.diff <- Q.comp - Q
  ##
  pval.Q.diff <- 1 - pchisq(Q.diff, df.Q.diff)

  
  res <- list(comparisons = comparisons,
              components = components,
              combinations = combinations,
              ##
              Q.comp = Q.comp,
              df.Q.comp = df.Q.comp,
              pval.Q.comp = pval.Q.comp,
              ##
              Q.diff = Q.diff,
              df.Q.diff = df.Q.diff,
              pval.Q.diff = pval.Q.diff)
  
  
  res
}