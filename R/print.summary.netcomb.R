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
                                  scientific.pval = gs("scientific.pval"),
                                  big.mark = gs("big.mark"),
                                  ...) {
  
  
  meta:::chkclass(x, "summary.netcomb")
  ##  
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  formatN <- meta:::formatN
  formatPT <- meta:::formatPT
  
  
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
  ##
  chklogical(scientific.pval)
  
  
  if (comb.fixed | comb.random) {
    cat(paste("Number of studies: k = ", x$k, "\n", sep = ""))
    cat(paste("Number of treatments: n = ", x$n, "\n", sep = ""))
    cat(paste("Number of components: c = ", x$c, "\n", sep = ""))
    cat(paste("Number of pairwise comparisons: m = ", x$m, "\n", sep = ""))
    ##
    cat("\n")
  }
  
  
  trts <- rownames(x$x$TE.fixed)
  trts.abbr <- treats(trts, nchar.trts)
  ##
  comps <- names(x$components.fixed$TE)
  comps.abbr <- treats(comps, nchar.trts)
  
  
  dat1.f <- prcombs(x$combinations.fixed,
                    backtransf, x$sm, x$level,
                    trts, trts.abbr,
                    digits, digits.zval, digits.pval.Q,
                    scientific.pval, big.mark)
  ##
  dat1.r <- prcombs(x$combinations.random,
                    backtransf, x$sm, x$level,
                    trts, trts.abbr,
                    digits, digits.zval, digits.pval.Q,
                    scientific.pval, big.mark)
  ##
  if (comb.fixed) {
    cat("Results for combinations (componentwise analysis, fixed effects model):\n")
    print(dat1.f)
    cat("\n")
  }
  ##
  if (comb.random) {
    cat("Results for combinations (componentwise analysis, random effects model):\n")
    print(dat1.r)
    cat("\n")
  }
  
  
  dat2.f <- prcombs(x$components.fixed,
                    backtransf, x$sm, x$level,
                    comps, comps.abbr,
                    digits, digits.zval, digits.pval.Q,
                    scientific.pval, big.mark)
  ##
  dat2.r <- prcombs(x$components.random,
                    backtransf, x$sm, x$level,
                    comps, comps.abbr,
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
    cat("\n")
  }
  
  
  cat("Heterogeneity statistics:\n")
  
  
  print(data.frame(Q = formatN(c(x$Q.comp.fixed, x$Q, x$Q.diff.fixed),
                               digits.Q),
                   df.Q = c(x$df.Q.comp, x$df.Q, x$df.Q.diff),
                   pval = formatPT(c(x$pval.Q.comp.fixed, x$pval.Q,
                                     x$pval.Q.diff.fixed),
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
