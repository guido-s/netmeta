formatCC <- function(x,
                     backtransf, sm, level,
                     comps, comps.abbr, sep.comps,
                     digits, digits.stat, digits.pval.Q,
                     scientific.pval, zero.pval, JAMA.pval,
                     big.mark) {
  
  
  relative <- is_relative_effect(sm)
  
  sm.lab <- sm
  ##
  if (sm != "") {
    sm.lab <- paste0("i", if (sm == "VE") "VR" else sm)
    if (!backtransf & relative)
      sm.lab <- paste0("log(", sm.lab, ")")
  }
  ##  
  ci.lab <- paste0(round(100 * level, 1), "%-CI")
  
  
  ## First column contains row names
  ##
  res <- x
  ##
  rownames(res) <- compos(rownames(res), comps, comps.abbr, sep.comps)
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
  res$TE <- formatN(res$TE, digits, "NA", big.mark)
  res$lower <- formatCI(formatN(round(res$lower, digits),
                                digits, "NA", big.mark),
                        formatN(round(res$upper, digits),
                                digits, "NA", big.mark))
  res$statistic <- formatN(round(res$statistic, digits.stat),
                           digits.stat, big.mark = big.mark)
  res$p <- formatPT(res$p,
                    digits = digits.pval.Q,
                    scientific = scientific.pval,
                    zero = zero.pval,
                    JAMA = JAMA.pval)
  ##
  res$seTE <- res$upper <- res$z <- NULL
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
  
  
  res
}
