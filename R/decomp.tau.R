decomp.tau <- function(x, tau.preset = 0, warn = TRUE) {
  
  
  nmak <- nma.krahn(x, tau.preset = tau.preset)
  if (is.null(nmak)) {
    if (warn)
      warning("Only a single design in network meta-analysis.", call. = FALSE)
    return(invisible(NULL))
  }
  ##
  design <- nmak$design
  studies <- nmak$studies
  V <- nmak$V
  V.studies <- nmak$V.studies
  X.obs <- nmak$X.obs
  multicomp <- nmak$multicomp
  
  
  multiarm <- any(studies$narms > 2)
  
  
  designs <- unique(as.character(design$design))
  TE.net.minus <- matrix(NA, nrow = nrow(X.obs), ncol = length(designs))
  colnames(TE.net.minus) <- designs
  ##
  df.net.minus <- numeric(length(designs))
  names(df.net.minus) <- designs
  ##
  for (i in designs) {
    narms <- design$narms[design$design == i][1]
    Xb <- cbind(X.obs,matrix(0, nrow(X.obs), narms - 1))
    Xb[(rownames(X.obs) == i), ] <- 0
    Xb[(rownames(X.obs) == i), ncol(X.obs) + (1:(narms - 1))] <- diag(narms - 1)
    TE.net.minus[, i] <- Xb %*% ginv(t(Xb) %*% solve(V) %*% Xb) %*% t(Xb) %*% solve(V) %*% design$TE.dir
    df.net.minus[i] <- qr(Xb)$rank
  }
  ##
  rownames(TE.net.minus) <- design$design
  
  
  freq.d <- rep_len(NA, nrow(design))
  Q.het.design <- rep_len(NA, nrow(design))
  ##
  for (i in unique(sort(as.character(studies$design)))) {
    Vs <- V.studies[rownames(V.studies) == i, colnames(V.studies) == i]
    ##
    TE.i <- studies$TE[as.character(studies$design) == i]
    TE.dir.i <- studies$TE.dir[as.character(studies$design) == i]
    ##
    Q.het.design[which(design$design == i)[1]] <- t(TE.i - TE.dir.i) %*% solve(Vs) %*% (TE.i - TE.dir.i)
    freq.d[which(design$design == i)[1]] <- length(TE.i)
  }
  ##
  NAfd <- is.na(freq.d)
  design$narms[NAfd] <- NA
  design$freq[NAfd] <- NA
  
  
  df.het.design <- freq.d - (design$narms - 1)
  pval.het.design <- pchisq(Q.het.design, df = df.het.design, lower.tail = FALSE)
  pval.het.design[df.het.design == 0] <- NA
  ##
  Q.het <- t(studies$TE - studies$TE.dir) %*% solve(V.studies) %*% (studies$TE - studies$TE.dir)
  ## Formula (7) in Krahn et al. (2013)
  Q.net <- t(studies$TE - studies$TE.net) %*% solve(V.studies) %*% (studies$TE - studies$TE.net)
  ## Formula (8) in Krahn et al. (2013)
  Q.inc <- Q.net - Q.het
  ##
  df.het <- sum(df.het.design, na.rm = TRUE)
  df.net <- sum(freq.d, na.rm = TRUE) - ncol(X.obs)
  df.inc <- df.net - df.het
  ##
  Qid <- t(design$TE.dir - design$TE.net) %*% solve(V) * (design$TE.dir - design$TE.net)
  nam <- colnames(Qid)
  Q.inc.design <- as.vector(Qid)
  names(Q.inc.design) <- nam
  
  
  residuals <- apply(TE.net.minus, 2, function(x) design$TE.dir - x)
  residuals <- residuals[, rownames(residuals)]
  diag(residuals) <- rep_len(0, nrow(residuals))
  ##
  if (multiarm)
    for (i in 1:length(multicomp))
      residuals[rownames(residuals) == multicomp[i], colnames(residuals) == multicomp[i]] <- 0
  
  
  Q.inc.detach <- apply(residuals, 2, function(x) t(x) %*% solve(V) %*% x)
  Q.inc.detach[NAfd] <- NA
  ##
  df.inc.detach <- rep_len(NA, length(NAfd))
  df.inc.detach[!NAfd]  <- nrow(X.obs) - df.net.minus 
  ##
  pval.inc.detach <- pchisq(Q.inc.detach, df = df.inc.detach, lower.tail = FALSE)
  pval.inc.detach[df.inc.detach == 0] <- NA
  
  
  Q.decomp <- data.frame(Q = c(Q.net, Q.het, Q.inc),
                         df = c(df.net, df.het, df.inc),
                         pval = pchisq(c(Q.net, Q.het, Q.inc),
                                       df = c(df.net, df.het, df.inc),
                                       lower.tail = FALSE))
  ##
  Q.decomp$pval[c(df.net, df.het, df.inc) == 0] <- NA
  ##
  rownames(Q.decomp) <- c("Total",
                          "Within designs",
                          "Between designs")
  
  
  Q.het.design <- data.frame(design = design$design,
                             Q = round(Q.het.design, 15),
                             df = df.het.design,
                             pval = pval.het.design)
  ##
  Q.inc.detach <- data.frame(design = design$design,
                             Q = Q.inc.detach,
                             df = df.inc.detach,
                             pval = pval.inc.detach)
  ##
  if (multiarm) {
    Q.het.design <- Q.het.design[!(duplicated(Q.het.design$design)), ]
    Q.inc.detach <- Q.inc.detach[!(duplicated(Q.inc.detach$design)), ]
  }
  ##
  if (any(is.na(Q.inc.detach$Q)))
    Q.inc.detach <- Q.inc.detach[!(is.na(Q.inc.detach$Q)), ]
  
  res <- list(
    Q.decomp = Q.decomp,
    Q.het.design = Q.het.design,
    Q.inc.detach = Q.inc.detach,
    Q.inc.design = Q.inc.design,
    residuals.inc.detach = residuals
  )
  
  res
}
