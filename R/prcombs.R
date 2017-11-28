prcombs <- function(x,
                    backtransf, sm, level,
                    trts, trts.abbr,
                    digits, digits.zval, digits.pval.Q,
                    scientific.pval, big.mark) {
  
  
  sm.lab <- sm
  ##
  relative <- meta:::is.relative.effect(sm)
  ##
  if (!backtransf & relative)
    sm.lab <- paste("log", sm, sep = "")
  ##  
  ci.lab <- paste(round(100 * level, 1), "%-CI", sep = "")
  
  
  res <- as.data.frame(x, stringsAsFactors = FALSE)
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
  res$TE <- meta:::format.NA(res$TE, digits, "NA", big.mark)
  res$lower <- meta:::p.ci(meta:::format.NA(round(res$lower, digits),
                                            digits, "NA", big.mark),
                           meta:::format.NA(round(res$upper, digits),
                                            digits, "NA", big.mark))
  res$z <- meta:::format.NA(res$z, digits.zval, big.mark = big.mark)
  res$p <- meta:::format.p(res$p,
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
