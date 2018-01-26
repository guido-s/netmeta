prcombs <- function(x,
                    backtransf, sm, level,
                    trts, trts.abbr,
                    digits, digits.zval, digits.pval.Q,
                    scientific.pval, big.mark,
                    seq = NULL) {
  
  
  formatN <- meta:::formatN
  
  
  sm.lab <- sm
  ##
  relative <- meta:::is.relative.effect(sm)
  ##
  if (!backtransf & relative)
    sm.lab <- paste("log", sm, sep = "")
  ##  
  ci.lab <- paste(round(100 * level, 1), "%-CI", sep = "")
  
  
  res <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!is.null(seq))
    res <- res[seq, ]
  ##
  if (backtransf & relative) {
    res$TE <- exp(res$TE)
    res$lower <- exp(res$lower)
    res$upper <- exp(res$upper)
  }
  ##
  rownames(res) <- as.character(factor(rownames(res),
                                       levels = trts, labels = trts.abbr))
  ##
  res$TE <- formatN(res$TE, digits, "NA", big.mark)
  res$lower <- meta:::formatCI(formatN(round(res$lower, digits),
                                       digits, "NA", big.mark),
                               formatN(round(res$upper, digits),
                                       digits, "NA", big.mark))
  res$z <- formatN(res$z, digits.zval, big.mark = big.mark)
  res$p <- meta:::formatPT(res$p,
                           digits = digits.pval.Q,
                           scientific = scientific.pval)
  ##
  res$seTE <- NULL
  res$upper <- NULL
  res$level <- NULL
  res$df <- NULL
  res$null.effect <- NULL
  ##
  sel <- names(res) == "TE"
  names(res)[sel] <- sm.lab
  ##
  sel <- names(res) == "lower"
  names(res)[sel] <- ci.lab
  
  
  res
}
