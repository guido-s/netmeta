prepare <- function(TE, seTE, treat1, treat2, studlab, tau = 0,
                    func.inverse) {
  
  if (is.na(tau))
    tau <- 0
  
  weights <- 1 / (seTE^2 + tau^2)
  
  data <- data.frame(studlab,
                     treat1, treat2,
                     treat1.pos = NA, treat2.pos = NA,
                     TE, seTE, weights,
                     narms = NA, stringsAsFactors = FALSE)
  ##
  ## Ordering dataset
  ##
  o <- order(data$studlab, data$treat1, data$treat2)
  data <- data[o, ]
  ##
  ## Adapt numbers to treatment IDs
  ##
  names.treat <- sort(unique(c(data$treat1, data$treat2)))
  data$treat1.pos <- match(data$treat1, names.treat)
  data$treat2.pos <- match(data$treat2, names.treat)
  
  newdata <- data[1, ][-1, ]
  ##
  sl <- unique(data$studlab)
  ##
  ## Determining number of arms and adjusting weights of
  ## multi-arm studies
  ##
  for (s in sl) {
    subgraph <- data[data$studlab == s, ]
    subgraph$narms <- (1 + sqrt(8 * dim(subgraph)[1] + 1)) / 2
    ## Reciprocal new weights
    if (dim(subgraph)[1] > 1)
      subgraph$weights <-
        1 / multiarm(1 / subgraph$weights, s, func.inverse)$v
    ##
    newdata <- rbind(newdata, subgraph)
  }
  res <- newdata
  ##
  res$order <- o
  ##
  res
}


prepare2 <- function(TE, seTE, treat1, treat2, studlab, tau = 0,
                     correlated, func.inverse) {
  
  if (is.na(tau))
    tau <- 0
  
  data <- data.frame(studlab,
                     treat1, treat2,
                     treat1.pos = NA, treat2.pos = NA,
                     TE, seTE, weights = 1 / (seTE^2 + tau^2), correlated,
                     narms = NA, stringsAsFactors = FALSE)
  #
  # Ordering dataset
  #
  o <- order(data$studlab, data$treat1, data$treat2)
  data <- data[o, ]
  #
  # Adapt numbers to treatment IDs
  #
  names.treat <- sort(unique(c(data$treat1, data$treat2)))
  data$treat1.pos <- match(data$treat1, names.treat)
  data$treat2.pos <- match(data$treat2, names.treat)
  #
  data$order <- o
  
  sl <- unique(data$studlab)
  #
  # List with weight matrices
  #
  W.list <- vector("list", length(sl))
  names(W.list) <- sl
  #
  # Determining number of arms and adjusting weights of
  # multi-arm studies
  #
  for (s in sl) {
    sel.s <- data$studlab == s
    correlated.s <- unique(data$correlated[sel.s]) 
    #
    if (length(correlated.s) != 1)
      stop("Different values for argument 'correlated' for study '", s, "'.",
           call. = FALSE)
    # Only treatment arms from multi-arm studies can be correlated
    if (correlated.s & sum(sel.s) == 1)
      correlated.s <- FALSE
    #
    res.s <- covar_study(1 / data$weights[sel.s], s, correlated.s, func.inverse)
    #
    W.list[[s]] <- res.s$W
    data$narms[sel.s] <- res.s$n
    data$weights[sel.s] <- diag(res.s$W)
  }
  #
  res <- list(W = bdiag(W.list), data = data)
  #
  res
}


covar_study <- function(v, studlab, correlated, func.inverse) {
  m <- length(v)
  n <- (1 + sqrt(8 * m + 1)) / 2
  #
  if (correlated) {
    B <- createB(ncol = n)
    V <- diag(diag(t(B) %*% diag(v, nrow = m) %*% B)) - t(B) %*%
      diag(v, nrow = m) %*% B
    #
    Cov <- matrix(0, nrow = m, ncol = m)
    edges <- matrix(nrow = m, ncol = 2)
    #
    r <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        r <- r + 1
        edges[r, ] <- c(i, j)
      }
    }
    #
    for (p in 1:(m - 1)) {
      i <- edges[p, 1]
      j <- edges[p, 2]
      #
      for (q in (p+1):m) {
        k <- edges[q, 1]
        l <- edges[q, 2] 
        #
        Cov[p, q] <- 0.5 * (V[i, l] - V[i, k] + V[j, k] - V[j, l])
        Cov[q, p] <- 0.5 * (V[i, l] - V[i, k] + V[j, k] - V[j, l])
      }
    }
    #
    for (p in 1:m) {
      i <- edges[p, 1]
      j <- edges[p, 2]
      #
      Cov[p, p] <- V[i, j]
    }
    #
    if (qr(Cov)$rank == n - 1)
      W <- ginv(as.matrix(Cov))
    else {
      if (length(v) > 1) {
        Cov <- diag(v)
        W <- diag(1 / v)
      }
      else {
        Cov <- matrix(v)
        W <- 1 / Cov
      }
    }
  }
  else {
    if (length(v) > 1) {
      v <- multiarm(v, studlab, func.inverse)$v
      Cov <- diag(v)
      W <- diag(1 / v)
    }
    else {
      Cov <- matrix(v)
      W <- 1 / Cov
    }
  }
  #
  res <- list(v = v, n = n, m = m, Cov = as.matrix(Cov), W = W)
}
