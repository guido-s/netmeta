nma.ruecker <- function(TE, W, seTE,
                        treat1, treat2,
                        treat1.pos, treat2.pos,
                        narms, studlab,
                        sm = "",
                        level, level.ma,
                        seTE.orig, tau.direct = 0, sep.trts = ":",
                        method.tau = "DL",
                        func.inverse) {
  
  
  w.pooled <- 1 / seTE^2
  
  
  m <- length(TE)                        # Number of pairwise comparisons (edges)
  n <- length(unique(c(treat1, treat2))) # Number of treatments (vertices)
  df1 <- 2 * sum(1 / narms)              # Sum of degrees of freedom per study
  
  
  ##
  ## B is the edge-vertex incidence matrix (m x n)
  ##
  B <- createB(treat1.pos, treat2.pos, ncol = n)
  ##
  ## B.full is the full edge-vertex incidence matrix (m x n)
  ##
  B.full <- createB(ncol = n)
  ##
  ## M is the unweighted Laplacian, D its diagonal,
  ## A is the adjacency matrix
  ##
  M <- t(B) %*% B    # unweighted Laplacian matrix
  D <- diag(diag(M)) # diagonal matrix
  A <- D - M         # adjacency matrix (n x n)
  ##
  ## L is the weighted Laplacian (Kirchhoff) matrix (n x n)
  ## Lplus is its Moore-Penrose pseudoinverse
  ##
  L <- t(B) %*% W %*% B
  Lplus <- do.call(func.inverse, list(X = L))
  Lplus[is.zero(Lplus)] <- 0
  ##
  ## R resistance distance (variance) matrix (n x n)
  ##
  R <- matrix(0, nrow = n, ncol = n) 
  for (i in 1:n) {
    for (j in 1:n) {
      R[i, j] <- Lplus[i, i] + Lplus[j, j] - 2 * Lplus[i, j]
    }
  }
  ##
  ## V is the vector of effective variances
  ##
  V <- vector(length = m, mode = "numeric")
  for (i in 1:m) {
    V[i] <- R[treat1.pos[i], treat2.pos[i]]
  }
  ##
  ## G is the matrix B %*% Lplus %*% t(B)
  ## H is the projection matrix (also called "hat matrix")
  ##
  ## Interpretation:
  ## (i)    diag(G) = V                 The effective variances
  ## (ii)   diag(H) = V %*% W = V * w   The leverages
  ## (iii)  sum(diag(H)) = n - 1        Rank of projection
  ## (iv)   mean(diag(H)) = (n - 1) / m Mean leverage = average efficiency
  ##
  G <- B %*% Lplus %*% t(B)
  H <- G %*% W
  ##
  ## Variance-covariance matrix for all comparisons
  ##
  Cov <- B.full %*% Lplus %*% t(B.full)
  ##
  ## Resulting effects and variances at numbered edges
  ##
  v <- as.vector(H %*% TE)
  ci.v <- ci(v, sqrt(V), level = level)
  ##
  ## Resulting effects, all edges, as a n x n matrix:
  ##
  all <- matrix(NA, nrow = n, ncol = n)
  ##
  for (i in 1:m) {
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
  ## Test of total heterogeneity / inconsistency:
  ##
  Q <- as.vector(t(TE - v) %*% W %*% (TE - v))
  df <- df1 - (n - 1)
  if (df == 0)
    pval.Q <-  NA
  else
    pval.Q <- pchisq(Q, df, lower.tail = FALSE)
  ##
  ## Heterogeneity variance
  ##
  I <- diag(m)
  E <- matrix(0, nrow = m, ncol = m)
  for (i in 1:m)
    for (j in 1:m)
      E[i, j] <- as.numeric(studlab[i] == studlab[j])
  ##
  if (df == 0) {
    tau2 <- NA
    tau <- NA
    I2 <- lower.I2 <- upper.I2 <- NA
  }
  else {
    tau2 <-
      max(0, (Q - df) / sum(diag((I - H) %*% (B %*% t(B) * E / 2) %*% W)))
    tau <- sqrt(tau2)
    ci.I2 <- isquared(Q, df, level.ma)
    I2 <- ci.I2$TE
    lower.I2 <- ci.I2$lower
    upper.I2 <- ci.I2$upper
  }
  ##
  ## Decomposition of total Q into parts from pairwise meta-analyses
  ## and residual inconsistency
  ##
  Q.matrix <- matrix(0, nrow = n, ncol = n)
  n.pairwise <- 0
  data <- data.frame(treat1.pos, treat2.pos,
                     TE, w.pooled,
                     stringsAsFactors = FALSE)
  ##
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (A[i, j] > 1) {
        n.pairwise <- n.pairwise + 1
        sub <- data[(data$treat1.pos == i & data$treat2.pos == j) |
                    (data$treat1.pos == j & data$treat2.pos == i), ]
        l <- nrow(sub)
        Q.matrix[i, j] <-
          suppressWarnings(metagen(TE = sub$TE, seTE = 1 / sqrt(sub$w.pooled),
                                   method.tau = "DL", method.tau.ci = "")$Q)
      }
    }
  }
  q <- t1 <- t2 <- dfs <- vector(length = n.pairwise, mode = "numeric")
  p <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (A[i, j] > 0) {
        p <- p + 1
        t1[p] <- i
        t2[p] <- j
        q[p] <- Q.matrix[i, j]
        dfs[p] <- A[i, j] - 1
      }
    }
  }
  ##
  Q.heterogeneity <- sum(Q.matrix)
  Q.inconsistency <- Q - Q.heterogeneity
  ##
  names.treat <- sort(unique(c(treat1, treat2)))
  ##
  Q.decomp <- data.frame(treat1 = names.treat[t1],
                         treat2 = names.treat[t2],
                         Q = q,
                         df = dfs,
                         pval.Q = pchisq(q, dfs, lower.tail = FALSE),
                         stringsAsFactors = FALSE)
  ##
  Q.decomp$pval.Q[Q.decomp$df == 0] <- NA
  
  
  TE.pooled <- all
  seTE.pooled <- sqrt(R)
  ##
  ci.pooled <- meta::ci(all, sqrt(R), level = level.ma)
  ##
  lower.pooled <- ci.pooled$lower
  upper.pooled <- ci.pooled$upper
  statistic.pooled <- ci.pooled$statistic
  pval.pooled <- ci.pooled$p
  ##
  rownames(TE.pooled) <- colnames(TE.pooled) <- names.treat
  rownames(seTE.pooled) <- colnames(seTE.pooled) <- names.treat
  rownames(lower.pooled) <- colnames(lower.pooled) <- names.treat
  rownames(upper.pooled) <- colnames(upper.pooled) <- names.treat
  rownames(statistic.pooled) <- colnames(statistic.pooled) <- names.treat
  rownames(pval.pooled) <- colnames(pval.pooled) <- names.treat
  
  
  A.matrix <- A
  L.matrix <- L
  Lplus.matrix <- Lplus
  ##
  rownames(A.matrix) <- colnames(A.matrix) <- names.treat
  rownames(L.matrix) <- colnames(L.matrix) <- names.treat
  rownames(Lplus.matrix) <- colnames(Lplus.matrix) <- names.treat
  rownames(Q.matrix) <- colnames(Q.matrix) <- names.treat
  
  G.matrix <- G
  H.matrix <- H
  ##
  rownames(G.matrix) <- colnames(G.matrix) <- studlab
  rownames(H.matrix) <- colnames(H.matrix) <- studlab
  
  names.Cov <- rep("", nrow(B.full))
  ##
  k <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      k <- k + 1
      names.Cov[k] <- paste(names.treat[i], names.treat[j], sep = sep.trts)
    }
  }
  ##
  rownames(Cov) <- colnames(Cov) <- names.Cov
  
  
  ## Calculate direct treatment estimates
  ##
  Bo <- B
  Bo[Bo == 1] <- -2
  rownames(B) <- studlab
  colnames(B) <- names.treat
  ##
  B.matrix <- B[do.call(order, data.frame(Bo)), , drop = FALSE]
  ##
  TE.direct <- seTE.direct <-
    Q.direct <- tau2.direct <- I2.direct <-
      matrix(NA, ncol = length(names.treat), nrow = length(names.treat))
  ##
  rownames(TE.direct) <- colnames(TE.direct) <-
    rownames(seTE.direct) <- colnames(seTE.direct) <-
    rownames(Q.direct) <- colnames(Q.direct) <-
    rownames(tau2.direct) <- colnames(tau2.direct) <-
    rownames(I2.direct) <- colnames(I2.direct) <-
    names.treat
  ##
  C.matrix <- B.matrix[!duplicated(B.matrix), , drop = FALSE]
  ##
  for (i in 1:dim(C.matrix)[1]) {
    sel1 <- C.matrix[i, ] ==  1
    sel2 <- C.matrix[i, ] == -1
    sel.treat1 <- colnames(C.matrix)[sel1]
    sel.treat2 <- colnames(C.matrix)[sel2]
    ##
    selstud <- treat1 == sel.treat1 & treat2 == sel.treat2
    ##
    m.i.tau.preset <-
      suppressWarnings(metagen(TE, seTE.orig, subset = selstud,
                               tau.preset = tau.direct,
                               method.tau = "DL", method.tau.ci = ""))
    ##
    if (is.na(tau.direct) | tau.direct == 0) {
      TE.i   <- m.i.tau.preset$TE.common
      seTE.i <- m.i.tau.preset$seTE.common
    }
    else {
      TE.i   <- m.i.tau.preset$TE.random
      seTE.i <- m.i.tau.preset$seTE.random
    }
    ##
    TE.direct[sel.treat1, sel.treat2]   <- TE.i
    seTE.direct[sel.treat1, sel.treat2] <- seTE.i
    ##
    TE.direct[sel.treat2, sel.treat1]   <- -TE.i
    seTE.direct[sel.treat2, sel.treat1] <- seTE.i
    ##
    m.i <-
      suppressWarnings(metagen(TE, seTE.orig, subset = selstud,
                               method.tau = method.tau, method.tau.ci = ""))
    ##
    Q.direct[sel.treat1, sel.treat2] <-
      Q.direct[sel.treat2, sel.treat1] <- m.i$Q
    tau2.direct[sel.treat1, sel.treat2] <-
      tau2.direct[sel.treat2, sel.treat1] <- m.i$tau2
    I2.direct[sel.treat1, sel.treat2] <-
      I2.direct[sel.treat2, sel.treat1] <- m.i$I2
    
  }
  ##
  ci.direct <- meta::ci(TE.direct, seTE.direct, level = level.ma)
  ##
  lower.direct <- ci.direct$lower
  upper.direct <- ci.direct$upper
  statistic.direct <- ci.direct$statistic
  pval.direct <- ci.direct$p
  
  
  ##
  ## Contribution of individual studies to Q
  ##
  Q.pooled <- w.pooled * (TE - v)^2
  
  
  res <- list(studlab = studlab,
              treat1 = treat1, treat2 = treat2,
              TE = TE, seTE = seTE,
              seTE.orig = seTE.orig,
              TE.nma = ci.v$TE,
              seTE.nma = ci.v$seTE,
              lower.nma = ci.v$lower,
              upper.nma = ci.v$upper,
              statistic.nma = ci.v$statistic,
              pval.nma = ci.v$p,
              leverage = diag(H),
              w.pooled = w.pooled,
              Q.pooled = Q.pooled,
              treat1.pos = treat1.pos,
              treat2.pos = treat2.pos,
              ##
              TE.pooled = TE.pooled,
              seTE.pooled = seTE.pooled,
              lower.pooled = lower.pooled,
              upper.pooled = upper.pooled,
              statistic.pooled = statistic.pooled,
              pval.pooled = pval.pooled,
              ##
              k = length(unique(studlab)),
              m = length(TE),
              n = dim(TE.pooled)[[1]],
              Q = Q,
              df = df,
              pval.Q = pval.Q,
              I2 = I2,
              lower.I2 = lower.I2,
              upper.I2 = upper.I2,
              tau = tau,
              Q.heterogeneity = Q.heterogeneity,
              Q.inconsistency = Q.inconsistency,
              ##
              sm = sm,
              level = level,
              level.ma = level.ma,
              ##
              A.matrix = A.matrix,
              B.matrix = B,
              L.matrix = L.matrix,
              Lplus.matrix = Lplus.matrix,
              Q.matrix = Q.matrix,
              G.matrix = G.matrix,
              H.matrix = H.matrix,
              ##
              Cov = Cov,
              ##
              TE.direct = TE.direct,
              seTE.direct = seTE.direct,
              lower.direct = lower.direct,
              upper.direct = upper.direct,
              statistic.direct = statistic.direct,
              pval.direct = pval.direct,
              ##
              Q.direct = Q.direct,
              tau2.direct = tau2.direct,
              I2.direct = I2.direct,
              ##
              Q.decomp = Q.decomp
              )
  
  res$version <- packageDescription("netmeta")$Version
  
  res
}


nma_ruecker <- function(TE, W, seTE,
                        treat1, treat2,
                        treat1.pos, treat2.pos,
                        narms, studlab,
                        sm = "",
                        level, level.ma,
                        seTE.orig, tau.direct = 0, sep.trts = ":",
                        method.tau = "DL",
                        func.inverse) {
  
  
  w.pooled <- 1 / seTE^2
  
  m <- length(TE)                        # Number of pairwise comparisons (edges)
  n <- length(unique(c(treat1, treat2))) # Number of treatments (vertices)
  df1 <- 2 * sum(1 / narms)              # Sum of degrees of freedom per study
  
  # Drop Matrix attributes
  W <- as.matrix(W)
  class(W) <- "matrix"
  
  ##
  ## B is the edge-vertex incidence matrix (m x n)
  ##
  B <- createB(treat1.pos, treat2.pos, ncol = n)
  ##
  ## B.full is the full edge-vertex incidence matrix (m x n)
  ##
  B.full <- createB(ncol = n)
  ##
  ## M is the unweighted Laplacian, D its diagonal,
  ## A is the adjacency matrix
  ##
  M <- t(B) %*% B    # unweighted Laplacian matrix
  D <- diag(diag(M)) # diagonal matrix
  A <- D - M         # adjacency matrix (n x n)
  ##
  ## L is the weighted Laplacian (Kirchhoff) matrix (n x n)
  ## Lplus is its Moore-Penrose pseudoinverse
  ##
  L <- t(B) %*% W %*% B
  Lplus <- do.call(func.inverse, list(X = L))
  Lplus[is.zero(Lplus)] <- 0
  ##
  ## R resistance distance (variance) matrix (n x n)
  ##
  R <- matrix(0, nrow = n, ncol = n) 
  for (i in 1:n) {
    for (j in 1:n) {
      R[i, j] <- Lplus[i, i] + Lplus[j, j] - 2 * Lplus[i, j]
    }
  }
  ##
  ## V is the vector of effective variances
  ##
  V <- vector(length = m, mode = "numeric")
  for (i in 1:m) {
    V[i] <- R[treat1.pos[i], treat2.pos[i]]
  }
  ##
  ## G is the matrix B %*% Lplus %*% t(B)
  ## H is the projection matrix (also called "hat matrix")
  ##
  ## Interpretation:
  ## (i)    diag(G) = V                 The effective variances
  ## (ii)   diag(H) = V %*% W = V * w   The leverages
  ## (iii)  sum(diag(H)) = n - 1        Rank of projection
  ## (iv)   mean(diag(H)) = (n - 1) / m Mean leverage = average efficiency
  ##
  G <- B %*% Lplus %*% t(B)
  H <- G %*% W
  ##
  ## Variance-covariance matrix for all comparisons
  ##
  Cov <- B.full %*% Lplus %*% t(B.full)
  ##
  ## Resulting effects and variances at numbered edges
  ##
  v <- as.vector(H %*% TE)
  ci.v <- ci(v, sqrt(V), level = level)
  ##
  ## Resulting effects, all edges, as a n x n matrix:
  ##
  all <- matrix(NA, nrow = n, ncol = n)
  ##
  for (i in 1:m) {
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
  ## Test of total heterogeneity / inconsistency:
  ##
  Q <- as.vector(t(TE - v) %*% W %*% (TE - v))
  df <- df1 - (n - 1)
  if (df == 0)
    pval.Q <-  NA
  else
    pval.Q <- pchisq(Q, df, lower.tail = FALSE)
  ##
  ## Heterogeneity variance
  ##
  I <- diag(m)
  E <- matrix(0, nrow = m, ncol = m)
  for (i in 1:m)
    for (j in 1:m)
      E[i, j] <- as.numeric(studlab[i] == studlab[j])
  ##
  if (df == 0) {
    tau2 <- NA
    tau <- NA
    I2 <- lower.I2 <- upper.I2 <- NA
  }
  else {
    tau2 <-
      max(0, (Q - df) / sum(diag((I - H) %*% (B %*% t(B) * E / 2) %*% W)))
    tau <- sqrt(tau2)
    ci.I2 <- isquared(Q, df, level.ma)
    I2 <- ci.I2$TE
    lower.I2 <- ci.I2$lower
    upper.I2 <- ci.I2$upper
  }
  ##
  ## Decomposition of total Q into parts from pairwise meta-analyses
  ## and residual inconsistency
  ##
  Q.matrix <- matrix(0, nrow = n, ncol = n)
  n.pairwise <- 0
  data <- data.frame(treat1.pos, treat2.pos,
                     TE, w.pooled,
                     stringsAsFactors = FALSE)
  ##
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (A[i, j] > 1) {
        n.pairwise <- n.pairwise + 1
        sub <- data[(data$treat1.pos == i & data$treat2.pos == j) |
                      (data$treat1.pos == j & data$treat2.pos == i), ]
        l <- nrow(sub)
        Q.matrix[i, j] <-
          suppressWarnings(metagen(TE = sub$TE, seTE = 1 / sqrt(sub$w.pooled),
                                   method.tau = "DL", method.tau.ci = "")$Q)
      }
    }
  }
  q <- t1 <- t2 <- dfs <- vector(length = n.pairwise, mode = "numeric")
  p <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (A[i, j] > 0) {
        p <- p + 1
        t1[p] <- i
        t2[p] <- j
        q[p] <- Q.matrix[i, j]
        dfs[p] <- A[i, j] - 1
      }
    }
  }
  ##
  Q.heterogeneity <- sum(Q.matrix)
  Q.inconsistency <- Q - Q.heterogeneity
  ##
  names.treat <- sort(unique(c(treat1, treat2)))
  ##
  Q.decomp <- data.frame(treat1 = names.treat[t1],
                         treat2 = names.treat[t2],
                         Q = q,
                         df = dfs,
                         pval.Q = pchisq(q, dfs, lower.tail = FALSE),
                         stringsAsFactors = FALSE)
  ##
  Q.decomp$pval.Q[Q.decomp$df == 0] <- NA
  
  
  TE.pooled <- all
  seTE.pooled <- sqrt(R)
  ##
  ci.pooled <- meta::ci(all, sqrt(R), level = level.ma)
  ##
  lower.pooled <- ci.pooled$lower
  upper.pooled <- ci.pooled$upper
  statistic.pooled <- ci.pooled$statistic
  pval.pooled <- ci.pooled$p
  ##
  rownames(TE.pooled) <- colnames(TE.pooled) <- names.treat
  rownames(seTE.pooled) <- colnames(seTE.pooled) <- names.treat
  rownames(lower.pooled) <- colnames(lower.pooled) <- names.treat
  rownames(upper.pooled) <- colnames(upper.pooled) <- names.treat
  rownames(statistic.pooled) <- colnames(statistic.pooled) <- names.treat
  rownames(pval.pooled) <- colnames(pval.pooled) <- names.treat
  
  
  A.matrix <- A
  L.matrix <- L
  Lplus.matrix <- Lplus
  ##
  rownames(A.matrix) <- colnames(A.matrix) <- names.treat
  rownames(L.matrix) <- colnames(L.matrix) <- names.treat
  rownames(Lplus.matrix) <- colnames(Lplus.matrix) <- names.treat
  rownames(Q.matrix) <- colnames(Q.matrix) <- names.treat
  
  G.matrix <- G
  H.matrix <- H
  ##
  rownames(G.matrix) <- colnames(G.matrix) <- studlab
  rownames(H.matrix) <- colnames(H.matrix) <- studlab
  
  names.Cov <- rep("", nrow(B.full))
  ##
  k <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      k <- k + 1
      names.Cov[k] <- paste(names.treat[i], names.treat[j], sep = sep.trts)
    }
  }
  ##
  rownames(Cov) <- colnames(Cov) <- names.Cov
  
  
  ## Calculate direct treatment estimates
  ##
  Bo <- B
  Bo[Bo == 1] <- -2
  rownames(B) <- studlab
  colnames(B) <- names.treat
  ##
  B.matrix <- B[do.call(order, data.frame(Bo)), , drop = FALSE]
  ##
  TE.direct <- seTE.direct <-
    Q.direct <- tau2.direct <- I2.direct <-
    matrix(NA, ncol = length(names.treat), nrow = length(names.treat))
  ##
  rownames(TE.direct) <- colnames(TE.direct) <-
    rownames(seTE.direct) <- colnames(seTE.direct) <-
    rownames(Q.direct) <- colnames(Q.direct) <-
    rownames(tau2.direct) <- colnames(tau2.direct) <-
    rownames(I2.direct) <- colnames(I2.direct) <-
    names.treat
  ##
  C.matrix <- B.matrix[!duplicated(B.matrix), , drop = FALSE]
  ##
  for (i in 1:dim(C.matrix)[1]) {
    sel1 <- C.matrix[i, ] ==  1
    sel2 <- C.matrix[i, ] == -1
    sel.treat1 <- colnames(C.matrix)[sel1]
    sel.treat2 <- colnames(C.matrix)[sel2]
    ##
    selstud <- treat1 == sel.treat1 & treat2 == sel.treat2
    ##
    m.i.tau.preset <-
      suppressWarnings(metagen(TE, seTE.orig, subset = selstud,
                               tau.preset = tau.direct,
                               method.tau = "DL", method.tau.ci = ""))
    ##
    if (is.na(tau.direct) | tau.direct == 0) {
      TE.i   <- m.i.tau.preset$TE.common
      seTE.i <- m.i.tau.preset$seTE.common
    }
    else {
      TE.i   <- m.i.tau.preset$TE.random
      seTE.i <- m.i.tau.preset$seTE.random
    }
    ##
    TE.direct[sel.treat1, sel.treat2]   <- TE.i
    seTE.direct[sel.treat1, sel.treat2] <- seTE.i
    ##
    TE.direct[sel.treat2, sel.treat1]   <- -TE.i
    seTE.direct[sel.treat2, sel.treat1] <- seTE.i
    ##
    m.i <-
      suppressWarnings(metagen(TE, seTE.orig, subset = selstud,
                               method.tau = method.tau, method.tau.ci = ""))
    ##
    Q.direct[sel.treat1, sel.treat2] <-
      Q.direct[sel.treat2, sel.treat1] <- m.i$Q
    tau2.direct[sel.treat1, sel.treat2] <-
      tau2.direct[sel.treat2, sel.treat1] <- m.i$tau2
    I2.direct[sel.treat1, sel.treat2] <-
      I2.direct[sel.treat2, sel.treat1] <- m.i$I2
    
  }
  ##
  ci.direct <- meta::ci(TE.direct, seTE.direct, level = level.ma)
  ##
  lower.direct <- ci.direct$lower
  upper.direct <- ci.direct$upper
  statistic.direct <- ci.direct$statistic
  pval.direct <- ci.direct$p
  
  
  ##
  ## Contribution of individual studies to Q
  ##
  Q.pooled <- w.pooled * (TE - v)^2
  
  
  res <- list(studlab = studlab,
              treat1 = treat1, treat2 = treat2,
              TE = TE, seTE = seTE,
              seTE.orig = seTE.orig,
              TE.nma = ci.v$TE,
              seTE.nma = ci.v$seTE,
              lower.nma = ci.v$lower,
              upper.nma = ci.v$upper,
              statistic.nma = ci.v$statistic,
              pval.nma = ci.v$p,
              leverage = diag(H),
              w.pooled = w.pooled,
              Q.pooled = Q.pooled,
              treat1.pos = treat1.pos,
              treat2.pos = treat2.pos,
              ##
              TE.pooled = TE.pooled,
              seTE.pooled = seTE.pooled,
              lower.pooled = lower.pooled,
              upper.pooled = upper.pooled,
              statistic.pooled = statistic.pooled,
              pval.pooled = pval.pooled,
              ##
              k = length(unique(studlab)),
              m = length(TE),
              n = dim(TE.pooled)[[1]],
              Q = Q,
              df = df,
              pval.Q = pval.Q,
              I2 = I2,
              lower.I2 = lower.I2,
              upper.I2 = upper.I2,
              tau = tau,
              Q.heterogeneity = Q.heterogeneity,
              Q.inconsistency = Q.inconsistency,
              ##
              sm = sm,
              level = level,
              level.ma = level.ma,
              ##
              A.matrix = A.matrix,
              B.matrix = B,
              L.matrix = L.matrix,
              Lplus.matrix = Lplus.matrix,
              Q.matrix = Q.matrix,
              G.matrix = G.matrix,
              H.matrix = H.matrix,
              ##
              Cov = Cov,
              ##
              TE.direct = TE.direct,
              seTE.direct = seTE.direct,
              lower.direct = lower.direct,
              upper.direct = upper.direct,
              statistic.direct = statistic.direct,
              pval.direct = pval.direct,
              ##
              Q.direct = Q.direct,
              tau2.direct = tau2.direct,
              I2.direct = I2.direct,
              ##
              Q.decomp = Q.decomp
  )
  
  res$version <- packageDescription("netmeta")$Version
  
  res
}
