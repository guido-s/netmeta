print.summary.netcomb <- function(x,
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
                                  digits.tau2 = gs("digits.tau2"),
                                  digits.I2 = gs("digits.I2"),
                                  scientific.pval = gs("scientific.pval"),
                                  big.mark = gs("big.mark"),
                                  ...) {
  
  
  ##
  ##
  ## (1) Check class and arguments
  ##
  ##
  meta:::chkclass(x, "summary.netcomb")
  ##  
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  formatN <- meta:::formatN
  formatPT <- meta:::formatPT
  ##  
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(backtransf)
  chknumeric(nchar.trts, min = 1, single = TRUE)
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.zval, min = 0, single = TRUE)
  chknumeric(digits.pval, min = 1, single = TRUE)
  chknumeric(digits.pval.Q, min = 1, single = TRUE)
  chknumeric(digits.Q, min = 0, single = TRUE)
  chknumeric(digits.tau2, min = 0, single = TRUE)
  chknumeric(digits.I2, min = 0, single = TRUE)
  ##
  chklogical(scientific.pval)
  
  
  I2 <- round(100 * x$I2, digits.I2)
  
  
  if (comb.fixed | comb.random) {
    cat(paste("Number of studies: k = ", x$k, "\n", sep = ""))
    cat(paste("Number of treatments: n = ", x$n, "\n", sep = ""))
    cat(paste("Number of active components: c = ", x$c, "\n", sep = ""))
    cat(paste("Number of pairwise comparisons: m = ", x$m, "\n", sep = ""))
    ##
    cat("\n")
  }
  
  
  trts <- x$trts
  trts.abbr <- treats(trts, nchar.trts)
  ##
  comps <- x$comps
  comps.abbr <- treats(comps, nchar.trts)
  
  
  dat1.f <- formatCC(x$combinations.fixed,
                     backtransf, x$sm, x$level, trts.abbr,
                     digits, digits.zval, digits.pval.Q,
                     scientific.pval, big.mark, x$seq)
  ##
  dat1.r <- formatCC(x$combinations.random,
                     backtransf, x$sm, x$level, trts.abbr,
                     digits, digits.zval, digits.pval.Q,
                     scientific.pval, big.mark, x$seq)
  ##
  if (comb.fixed) {
    cat("Results for combinations (additive model, fixed effects model):\n")
    print(dat1.f)
    cat("\n")
  }
  ##
  if (comb.random) {
    cat("Results for combinations (additive model, random effects model):\n")
    print(dat1.r)
    cat("\n")
  }
  
  
  dat2.f <- formatCC(x$components.fixed,
                     backtransf, x$sm, x$level, comps.abbr,
                     digits, digits.zval, digits.pval.Q,
                     scientific.pval, big.mark)
  ##
  dat2.r <- formatCC(x$components.random,
                     backtransf, x$sm, x$level, comps.abbr,
                     digits, digits.zval, digits.pval.Q,
                     scientific.pval, big.mark)
  ##
  if (comb.fixed) {
    cat("Results for components (fixed effects model):\n")
    print(dat2.f)
    cat("\n")
  }
  ##
  if (comb.random) {
    cat("Results for components (random effects model):\n")
    print(dat2.r)
  }
  
  
  cat(paste("\nQuantifying heterogeneity / inconsistency:\n",
            formatPT(x$tau^2,
                     lab = TRUE, labval = "tau^2",
                     digits = digits.tau2,
                     lab.NA = "NA", big.mark = big.mark),
            if (!is.na(I2))
              paste("; I^2 = ", round(I2, digits.I2), "%", "", sep = ""),
            "\n", sep = ""))
  
  
  cat("\nHeterogeneity statistics:\n")
  
  print(data.frame(Q = formatN(c(x$Q.additive,
                                 x$Q.standard,
                                 x$Q.diff),
                               digits.Q),
                   df.Q = formatN(c(x$df.Q.additive,
                                    x$df.Q.standard,
                                    x$df.Q.diff), 0),
                   pval = formatPT(c(x$pval.Q.additive,
                                     x$pval.Q.standard,
                                     x$pval.Q.diff),
                                   digits = digits.pval.Q,
                                   scientific = scientific.pval),
                   row.names = c("Additive model", "Standard model", "Difference")))
  
  
  if ((comb.fixed | comb.random)) {
    any.trts <- any(trts != trts.abbr)
    any.comps <- any(comps != comps.abbr)
    ##
    if (any.trts | any.comps)
      cat("\nLegend", if (any.trts & any.comps) "s", ":", sep = "")
    ##
    if (any.trts) {
      ##
      tmat <- data.frame(trts.abbr, trts)
      names(tmat) <- c("Abbreviation", "Treatment name")
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat("\n")
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(trts.abbr))) 
    }
    ##
    if (any.comps) {
      ##
      tmat <- data.frame(comps.abbr, comps)
      names(tmat) <- c("Abbreviation", " Component name")
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat("\n")
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(comps.abbr))) 
    }
  }
  
  
  invisible(NULL)
}
