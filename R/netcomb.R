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
                        Q, df.Q.additive, df.Q.diff,
                        x$n, x$sep.trts)
  
  
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
                        Q, df.Q.additive, df.Q.diff,
                        x$n, x$sep.trts)
  
  
  res <- list(studlab = x$studlab,
              treat1 = x$treat1,
              treat2 = x$treat2,
              ##
              TE = x$TE,
              seTE = x$seTE,
              seTE.adj = x$seTE.adj,
              ##
              event1 = x$event1,
              event2 = x$event2,
              n1 = x$n1,
              n2 = x$n2,
              ##
              k = x$k,
              m = m,
              n = n,
              d = x$d,
              c = c,
              ##
              trts = trts,
              k.trts = x$k.trts,
              n.trts = x$n.trts,
              events.trts = x$events.trts,
              ##
              studies = x$studies,
              narms = x$narms,
              ##
              designs = x$designs,
              ##
              comps = names(res.f$components$TE),
              k.comps = NA,
              n.comps = NA,
              events.comps = NA,
              ##
              TE.nma.fixed = x$TE.nma.fixed,
              seTE.nma.fixed = x$seTE.nma.fixed,
              lower.nma.fixed = x$lower.nma.fixed,
              upper.nma.fixed = x$upper.nma.fixed,
              zval.nma.fixed = x$zval.nma.fixed,
              pval.nma.fixed = x$pval.nma.fixed,
              ##
              leverage.fixed = x$leverage,
              w.fixed = x$w.fixed,
              Q.fixed = x$Q.fixed,
              ##
              TE.cnma.fixed = res.f$comparisons$TE,
              seTE.cnma.fixed = res.f$comparisons$seTE,
              lower.cnma.fixed = res.f$comparisons$lower,
              upper.cnma.fixed = res.f$comparisons$upper,
              zval.cnma.fixed = res.f$comparisons$z,
              pval.cnma.fixed = res.f$comparisons$p,
              ##
              TE.fixed = res.f$all.comparisons$TE,
              seTE.fixed = res.f$all.comparisons$seTE,
              lower.fixed = res.f$all.comparisons$lower,
              upper.fixed = res.f$all.comparisons$upper,
              zval.fixed = res.f$all.comparisons$z,
              pval.fixed = res.f$all.comparisons$p,
              ##
              TE.nma.random = x$TE.nma.random,
              seTE.nma.random = x$seTE.nma.random,
              lower.nma.random = x$lower.nma.random,
              upper.nma.random = x$upper.nma.random,
              zval.nma.random = x$zval.nma.random,
              pval.nma.random = x$pval.nma.random,
              ##
              w.random = x$w.random,
              ##
              TE.cnma.random = res.r$comparisons$TE,
              seTE.cnma.random = res.r$comparisons$seTE,
              lower.cnma.random = res.r$comparisons$lower,
              upper.cnma.random = res.r$comparisons$upper,
              zval.cnma.random = res.r$comparisons$z,
              pval.cnma.random = res.r$comparisons$p,
              ##
              TE.random = res.r$all.comparisons$TE,
              seTE.random = res.r$all.comparisons$seTE,
              lower.random = res.r$all.comparisons$lower,
              upper.random = res.r$all.comparisons$upper,
              zval.random = res.r$all.comparisons$z,
              pval.random = res.r$all.comparisons$p,
              ##
              Comp.fixed = unname(res.f$components$TE),
              seComp.fixed = unname(res.f$components$seTE),
              lower.Comp.fixed = unname(res.f$components$lower),
              upper.Comp.fixed = unname(res.f$components$upper),
              zval.Comp.fixed = unname(res.f$components$z),
              pval.Comp.fixed = unname(res.f$components$p),
              ##
              Comp.random = unname(res.r$components$TE),
              seComp.random = unname(res.r$components$seTE),
              lower.Comp.random = unname(res.r$components$lower),
              upper.Comp.random = unname(res.r$components$upper),
              zval.Comp.random = unname(res.r$components$z),
              pval.Comp.random = unname(res.r$components$p),
              ##
              Comb.fixed = unname(res.f$combinations$TE),
              seComb.fixed = unname(res.f$combinations$seTE),
              lower.Comb.fixed = unname(res.f$combinations$lower),
              upper.Comb.fixed = unname(res.f$combinations$upper),
              zval.Comb.fixed = unname(res.f$combinations$z),
              pval.Comb.fixed = unname(res.f$combinations$p),
              ##
              Comb.random = unname(res.r$combinations$TE),
              seComb.random = unname(res.r$combinations$seTE),
              lower.Comb.random = unname(res.r$combinations$lower),
              upper.Comb.random = unname(res.r$combinations$upper),
              zval.Comb.random = unname(res.r$combinations$z),
              pval.Comb.random = unname(res.r$combinations$p),
              ##
              Q.additive = Q.additive,
              df.Q.additive = df.Q.additive,
              pval.Q.additive = res.f$pval.Q.additive,
              tau = tau,
              I2 = I2,
              ##
              Q.standard = x$Q,
              df.Q.standard = x$df.Q,
              pval.Q.standard = x$pval.Q,
              ##
              Q.diff = res.f$Q.diff,
              df.Q.diff = df.Q.diff,
              pval.Q.diff = res.f$pval.Q.diff, 
              ##
              B.matrix = B.matrix,
              C.matrix = C.matrix,
              X = X,
              ##
              n.matrix = x$n.matrix,
              events.matrix = x$events.matrix,
              ##
              sm = x$sm,
              method = "Inverse",
              level = x$level,
              level.comb = x$level.comb,
              comb.fixed = x$comb.fixed,
              comb.random = x$comb.random,
              ##
              seq = x$seq,
              ##
              tau.preset = tau.preset,
              ##
              tol.multiarm = x$tol.multiarm,
              details.chkmultiarm = x$details.chkmultiarm,
              ##
              sep.trts = x$sep.trts,
              nchar.trts = x$nchar.trts,
              ##
              backtransf = x$backtransf,
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
