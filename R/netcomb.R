netcomb <- function(x,
                    inactive = NULL,
                    sep.components = "+",
                    C.matrix,
                    comb.fixed = x$comb.fixed,
                    comb.random = x$comb.random | !is.null(tau.preset),
                    tau.preset = NULL) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  meta:::chkclass(x, "netmeta")
  ##
  meta:::chkchar(sep.components, nchar = 1)
  ##
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  ##
  if (!is.null(tau.preset))
    meta:::chknumeric(tau.preset, min = 0, single = TRUE)
  
  
  ##
  ##
  ## (2) Get data
  ##
  ##
  n <- x$n # number of treatments / combinations
  m <- x$m # number of comparisons
  ##
  trts <- x$trts
  ##
  TE <- x$TE
  seTE <- x$seTE
  treat1 <- x$treat1
  treat2 <- x$treat2
  studlab <- x$studlab
  ##
  B.matrix <- x$B.matrix
  
  
  ##
  ##
  ## (3) Create C.matrix
  ##
  ##
  if (missing(C.matrix)) {
    ##
    ## Create C-matrix from netmeta object
    ##
    C.matrix <- as.matrix(createC(x, sep.components, inactive))
    C.matrix <- C.matrix[trts, , drop = FALSE]
  }
  else {
    ##
    ## Check if appropriate C.matrix is provided
    ##
    if (!(is.matrix(C.matrix) | is.data.frame(C.matrix)))
      stop("Argument 'C.matrix' must be a matrix or data frame.",
           call. = FALSE)
    ##
    if (nrow(C.matrix) != n)
      stop("Argument 'C.matrix' has wrong number of rows",
           " (must be equal to number of treatments).",
           call. = FALSE)
    ##
    ## Row names must be a permutation of treatments
    ##
    wrong.labels <- FALSE
    ##
    if (is.null(rownames(C.matrix)))
      wrong.labels <- TRUE
    else {
      if (length(unique(trts)) == length(unique(tolower(trts))))
        idx <- charmatch(tolower(rownames(C.matrix)),
                         tolower(trts), nomatch = NA)
      else
        idx <- charmatch(rownames(C.matrix), trts, nomatch = NA)
      ##
      if (any(is.na(idx)) || any(idx == 0))
        wrong.labels <- TRUE
    }
    ##
    if (wrong.labels)
      stop(paste("Row names of argument 'C.matrix' must be a ",
                 "permutation of treatment names:\n  ",
                 paste(paste("'", trts, "'", sep = ""),
                       collapse = " - "), sep = ""),
           call. = FALSE)
    ##
    C.matrix <- C.matrix[trts, , drop = FALSE]
    ##
    if (is.data.frame(C.matrix))
      C.matrix <- as.matrix(C.matrix)
  }
  ##
  c <- ncol(C.matrix) # number of components
  
  
  p0 <- prepare(TE, seTE, treat1, treat2, studlab)
  ##
  o <- order(p0$order)
  ##
  B.matrix <- createB(p0$treat1.pos[o], p0$treat2.pos[o])
  ##
  colnames(B.matrix) <- trts
  rownames(B.matrix) <- studlab
  
  
  ##
  ## Design matrix based on treatment components
  ##
  X <- B.matrix %*% C.matrix
  ##
  colnames(X) <- colnames(C.matrix)
  rownames(X) <- studlab
  
  
  ##
  ## Fixed effects models
  ##
  df.Q.additive <- x$df.Q + x$n - 1 - qr(X)$rank
  ##
  net <- netmeta(TE, seTE, treat1, treat2, studlab)
  ##
  Q <- net$Q
  df.Q <- net$df.Q
  pval.Q <- net$pval.Q
  ##
  df.Q.diff <- x$n - 1 - qr(X)$rank
  
  
  res.f <- nma.additive(p0$TE[o], p0$weights[o], p0$studlab[o],
                        p0$treat1[o], p0$treat2[o], x$level.comb,
                        X, C.matrix, B.matrix,
                        Q, df.Q.additive, df.Q.diff)
  
  
  ##
  ## Calculate heterogeneity statistics (additive model)
  ##
  Q.additive <- res.f$Q.additive
  ##
  if (!is.null(tau.preset))
    tau <- tau.preset
  else
    tau <- res.f$tau
  ##
  I2 <- res.f$I2
  ##
  tau2.calc <- if (is.na(tau)) 0 else tau^2
  
  
  ##
  ## Random effects models
  ##
  p1 <- prepare(TE, seTE, treat1, treat2, studlab, tau)
  res.r <- nma.additive(p1$TE[o], p1$weights[o], p1$studlab[o],
                        p1$treat1[o], p1$treat2[o], x$level.comb,
                        X, C.matrix, B.matrix,
                        Q, df.Q.additive, df.Q.diff)
  
  
  res <- list(k = x$k, n = n, m = m, c = c,
              ##
              comparisons.fixed = res.f$comparisons,
              components.fixed = res.f$components,
              combinations.fixed = res.f$combinations,
              ##
              comparisons.random = res.r$comparisons,
              components.random = res.r$components,
              combinations.random = res.r$combinations,
              ##
              tau = tau,
              I2 = I2,
              ##
              sm = x$sm,
              level.comb = x$level.comb,
              comb.fixed = x$comb.fixed,
              comb.random = x$comb.random,
              ##
              Q.additive = Q.additive,
              df.Q.additive = df.Q.additive,
              pval.Q.additive = res.f$pval.Q.additive,
              ##
              Q.standard = x$Q,
              df.Q.standard = x$df.Q,
              pval.Q.standard = x$pval.Q,
              ##
              Q.diff = res.f$Q.diff,
              df.Q.diff = df.Q.diff,
              pval.Q.diff = res.f$pval.Q.diff, 
              ##
              C.matrix = C.matrix, B.matrix = B.matrix, X = X,
              ##
              trts = x$trts,
              seq = x$seq,
              ##
              tau.preset = tau.preset,
              ##
              sep.trts = x$sep.trts,
              ##
              backtransf = x$backtransf,
              nchar.trts = x$nchar.trts,
              ##
              title = x$title,
              ##
              x = x,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- "netcomb"
  ##
  res
}
