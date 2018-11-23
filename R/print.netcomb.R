print.netcomb <- function(x,
                          comb.fixed = x$comb.fixed,
                          comb.random = x$comb.random,
                          backtransf = x$backtransf,
                          nchar.trts = x$nchar.trts,
                          ##
                          digits = gs("digits"),
                          digits.zval = gs("digits.zval"),
                          digits.pval = gs("digits.pval"),
                          digits.pval.Q = max(gs("digits.pval.Q"), 2),
                          digits.Q = gs("digits.Q"),
                          scientific.pval = gs("scientific.pval"),
                          big.mark = gs("big.mark"),
                          ...) {
  
  
  meta:::chkclass(x, "netcomb")
  
  
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  meta:::chklogical(backtransf)
  meta:::chknumeric(nchar.trts, min = 1, single = TRUE)
  ##
  meta:::chknumeric(digits, min = 0, single = TRUE)
  meta:::chknumeric(digits.zval, min = 0, single = TRUE)
  meta:::chknumeric(digits.pval, min = 1, single = TRUE)
  meta:::chknumeric(digits.pval.Q, min = 1, single = TRUE)
  meta:::chknumeric(digits.Q, min = 0, single = TRUE)
  ##
  meta:::chklogical(scientific.pval)
  
  
  trts <- x$trts
  trts.abbr <- treats(trts, nchar.trts)
  ##
  cnma.f <- data.frame(studlab = x$studlab,
                       treat1 = x$treat1,
                       treat2 = x$treat2,
                       TE = x$TE.cnma.fixed,
                       lower = x$lower.cnma.fixed,
                       upper = x$upper.cnma.fixed,
                       z = x$zval.cnma.fixed,
                       p = x$pval.cnma.fixed,
                       stringsAsFactors = FALSE)
  ##
  dat.f <- prcomps(cnma.f,
                   backtransf, x$sm, x$level.comb,
                   trts, trts.abbr,
                   digits, digits.zval, digits.pval.Q,
                   scientific.pval, big.mark)
  ##
  #dat.f <- dat.f[, !(colnames(dat.f) %in% c("z", "p"))]
  ##
  cnma.r <- data.frame(studlab = x$studlab,
                       treat1 = x$treat1,
                       treat2 = x$treat2,
                       TE = x$TE.cnma.random,
                       lower = x$lower.cnma.random,
                       upper = x$upper.cnma.random,
                       z = x$zval.cnma.random,
                       p = x$pval.cnma.random,
                       stringsAsFactors = FALSE)
  ##
  dat.r <- prcomps(cnma.r,
                   backtransf, x$sm, x$level.comb,
                   trts, trts.abbr,
                   digits, digits.zval, digits.pval.Q,
                   scientific.pval, big.mark)
  ##
  #dat.r <- dat.r[, !(colnames(dat.r) %in% c("z", "p"))]
  ##
  if (comb.fixed) {
    cat("Additive model (fixed effects model):\n")
    prmatrix(dat.f, quote = FALSE, right = TRUE, ...)
    cat("\n")
  }
  ##
  if (comb.random) {
    cat("Additive model (random effects model):\n")
    prmatrix(dat.r, quote = FALSE, right = TRUE, ...)
    cat("\n")
  }
  

  if (comb.fixed | comb.random)
    print(summary(x),
          comb.fixed = comb.fixed,
          comb.random = comb.random,
          backtransf = backtransf,
          nchar.trts = nchar.trts,
          ##
          digits = digits,
          digits.zval = digits.zval,
          digits.pval = digits.pval,
          digits.pval.Q = digits.pval.Q,
          digits.Q = digits.Q,
          scientific.pval = scientific.pval,
          big.mark = big.mark, ...)
  else
    cat("Please use argument 'comb.fixed = TRUE' or",
        "'comb.random = TRUE' to print meta-analysis results.\n",
        sep = "")
  
  
  invisible(NULL)
}
