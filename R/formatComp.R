formatComp <- function(x,
                       backtransf, sm, level,
                       comps, comps.abbr, sep.comps,
                       digits, digits.stat, digits.pval,
                       scientific.pval, zero.pval, JAMA.pval,
                       big.mark) {
  
  
  relative <- is.relative.effect(sm)
  
  
  sm.lab <- sm
  ##
  if (!backtransf & relative)
    sm.lab <- paste0("log", if (sm == "VE") "VR" else sm)
  ##  
  ci.lab <- paste(round(100 * level, 1), "%-CI", sep = "")
  
  
  ## First column contains row names
  ##
  rnam <- x[, 1]
  res <- x[, -1]
  ##
  if (backtransf) {
    res$TE <- backtransf(res$TE, sm)
    res$lower <- backtransf(res$lower, sm)
    res$upper <- backtransf(res$upper, sm)
    #
    # Switch lower and upper limit for VE if results have been
    # backtransformed
    #
    if (sm == "VE") {
      tmp.l <- res$lower
      res$lower <- res$upper
      res$upper <- tmp.l
    }
  }
  ##
  res$treat1 <- compos(res$treat1, comps, comps.abbr, sep.comps)
  res$treat2 <- compos(res$treat2, comps, comps.abbr, sep.comps)
  ##
  res$TE <- formatN(res$TE, digits, "NA", big.mark)
  res$lower <- formatCI(formatN(round(res$lower, digits),
                                       digits, "NA", big.mark),
                               formatN(round(res$upper, digits),
                                       digits, "NA", big.mark))
  res$statistic <- formatN(res$statistic, digits.stat, big.mark = big.mark)
  res$p <- formatPT(res$p,
                           digits = digits.pval,
                           scientific = scientific.pval,
                           zero = zero.pval,
                           JAMA = JAMA.pval)
  ##
  res$upper <- res$z <- NULL
  ##
  sel <- names(res) == "TE"
  names(res)[sel] <- sm.lab
  ##
  sel <- names(res) == "lower"
  names(res)[sel] <- ci.lab
  ##
  sel <- names(res) == "statistic"
  names(res)[sel] <- "z"
  ##
  sel <- names(res) == "p"
  names(res)[sel] <- "p-value"
  ##
  colnames <- names(res)
  ##
  res <- as.matrix(res)
  ##
  dimnames(res) <- list(rnam, colnames)
  
  
  res
}
