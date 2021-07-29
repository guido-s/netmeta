formatComp <- function(x,
                       backtransf, sm, level,
                       comps, comps.abbr, sep.comps,
                       digits, digits.stat, digits.pval.Q,
                       scientific.pval, big.mark) {
  
  
  formatN <- meta:::formatN
  relative <- meta:::is.relative.effect(sm)
  
  
  sm.lab <- sm
  ##
  if (!backtransf & relative)
    sm.lab <- paste("log", sm, sep = "")
  ##  
  ci.lab <- paste(round(100 * level, 1), "%-CI", sep = "")
  
  
  ## First column contains row names
  ##
  rnam <- x[, 1]
  res <- x[, -1]
  ##
  if (backtransf & relative) {
    res$TE <- exp(res$TE)
    res$lower <- exp(res$lower)
    res$upper <- exp(res$upper)
  }
  ##
  res$treat1 <- compos(res$treat1, comps, comps.abbr, sep.comps)
  res$treat2 <- compos(res$treat2, comps, comps.abbr, sep.comps)
  ##
  res$TE <- formatN(res$TE, digits, "NA", big.mark)
  res$lower <- meta:::formatCI(formatN(round(res$lower, digits),
                                       digits, "NA", big.mark),
                               formatN(round(res$upper, digits),
                                       digits, "NA", big.mark))
  res$statistic <- formatN(res$statistic, digits.stat, big.mark = big.mark)
  res$p <- meta:::formatPT(res$p,
                           digits = digits.pval.Q,
                           scientific = scientific.pval)
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
