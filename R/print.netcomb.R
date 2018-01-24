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
  
  
  trts <- rownames(x$x$TE.fixed)
  trts.abbr <- treats(trts, nchar.trts)
  ##
  dat.f <- prcomps(x$comparisons.fixed,
                   backtransf, x$sm, x$level.comb,
                   trts, trts.abbr,
                   digits, digits.zval, digits.pval.Q,
                   scientific.pval, big.mark)
  ##
  dat.r <- prcomps(x$comparisons.random,
                   backtransf, x$sm, x$level.comb,
                   trts, trts.abbr,
                   digits, digits.zval, digits.pval.Q,
                   scientific.pval, big.mark)
  ##
  if (comb.fixed) {
    cat("Componentwise analysis (fixed effects model):\n")
    prmatrix(dat.f, quote = FALSE, right = TRUE, ...)
    cat("\n")
  }
  ##
  if (comb.random) {
    cat("Componentwise analysis (random effects model):\n")
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
