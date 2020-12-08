##
##
## (Restricted) maximum log-likelihood function for NMA model
##
##

loglik_reml <- function(x,
                        TE, seTE, treat1.pos, treat2.pos,
                        X.matrix, n, method = "REML") {
  
  
  tau2 <- exp(x)
  w.pooled <- 1 / (seTE^2 + tau2)
  
  
  W <- diag(w.pooled, nrow = length(w.pooled)) # Diagonal matrix of weights
  
  
  ## L is the weighted Laplacian (Kirchhoff) matrix (n x n)
  ##
  L.matrix <- t(X.matrix) %*% W %*% X.matrix
  
  ## Lplus is its Moore-Penrose pseudoinverse
  ##
  Lplus <- ginv(L.matrix)
  
  G <- X.matrix %*% Lplus %*% t(X.matrix)
  H <- G %*% W
  
  ## Resulting effects at numbered edges
  ##
  v <- as.vector(H %*% TE)
  
  ## Resulting effects, all edges, as a n x n matrix:
  ##
  all <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:length(TE)) {
    all[treat1.pos[i], treat2.pos[i]] <- v[i]
  }
  ##
  for (i in 1:n) {
    for (j in 1:n) {
      for (k in 1:n) {
        if (!is.na(all[i, k]) & !is.na(all[j, k])) {
          all[i, j] <- all[i, k] - all[j, k]
          all[j, i] <- all[j, k] - all[i, k]
        }
        if (!is.na(all[i, j]) & !is.na(all[k, j])) {
          all[i, k] <- all[i, j] - all[k, j]
          all[k, i] <- all[k, j] - all[i, j]
        }
        if (!is.na(all[i, k]) & !is.na(all[i, j])) {
          all[j, k] <- all[i, k] - all[i, j]
          all[k, j] <- all[i, j] - all[i, k]
        }
      }
    }
  }
  ##
  ## Summary estimate of RE model
  ##
  theta_hat <- all[, 1]
  
  
  ## ML log-likelihood function (-2ML)
  ##
  res <- - log(det(abs(W))) + t(TE - X.matrix %*% c(theta_hat)) %*%
    W %*% (TE - X.matrix %*% c(theta_hat))
  ##
  ## REML log-likelihood function (-2RL)
  ##
  if (method == "REML")
    res <- res + log(det(abs(L.matrix)))
  
  res
}
