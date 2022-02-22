nettable_internal <- function(x, method,
                              ##
                              upper, reference.group, baseline.reference,
                              order, tol.direct, backtransf,
                              ##
                              digits, digits.I2, digits.pval,
                              ##
                              scientific.pval, zero.pval, JAMA.pval,
                              big.mark, text.NA,
                              bracket, separator, lower.blank, upper.blank,
                              ##
                              writexl,
                              ##
                              warn, verbose) {
  
  
  ##
  ##
  ## (1) Some checks and assignments
  ##
  ##
  is.bin <- inherits(x, "netmetabin")
  ##
  if (is.null(method))
    if (is.bin)
      method <- "SIDDE"
    else
      method <- "Back-calculation"
  ##
  if (is.null(reference.group))
    reference.group <- x$reference.group
  else
    reference.group <- setref(reference.group, x$trts)
  ##
  if (is.null(baseline.reference))
    baseline.reference <- x$baseline.reference
  ##
  if (is.null(backtransf))
    backtransf <- x$backtransf
  ##
  if (!is.null(order)) {
    order <- setseq(order, x$trts, equal.length = FALSE)
    baseline.reference <- FALSE
    reference.group <- ""
  }
  
  
  ##
  ##
  ## (2) Create dat.trts
  ##
  ##
  dat.trts <- comptrts(x, upper, reference.group, baseline.reference, order)
  
  
  ##
  ##
  ## (4) Change order of prop.direct.fixed and prop.direct.random
  ##
  ##
  if (!(is.bin & method == "SIDDE")) {
    prop.fixed <- sortprop(x$prop.direct.fixed, dat.trts, x$sep.trts)
    prop.random <- sortprop(x$prop.direct.random, dat.trts, x$sep.trts)
  }
  else
    prop.fixed <- prop.random <- NULL
  
  
  ##
  ##
  ## (5) Calculate / extract indirect estimates
  ##
  ##
  x.direct.indirect <- x
  ##
  if (method == "SIDDE") {
    sid <- sidde(x.direct.indirect, x$sep.trts, verbose, warn, FALSE)
    ##    
    x.direct.indirect$TE.indirect.fixed <- sid$TE.indirect.fixed
    x.direct.indirect$seTE.indirect.fixed <- sid$seTE.indirect.fixed
    ##
    if (!is.bin) {
      x.direct.indirect$TE.indirect.random <- sid$TE.indirect.random
      x.direct.indirect$seTE.indirect.random <- sid$seTE.indirect.random
    } 
  }
  ##
  direct.indirect <- direct.indirect(x.direct.indirect, tol.direct)
  
  
  ##
  ##
  ## (6) Transform matrices to data frames
  ##
  ##
  fixed <- mat2dat.table(direct.indirect, "fixed", dat.trts,
                         backtransf,
                         digits, digits.I2, digits.pval,
                         scientific.pval, zero.pval, JAMA.pval,
                         big.mark, text.NA,
                         writexl)
  ##
  random <- mat2dat.table(direct.indirect, "random", dat.trts,
                          backtransf,
                          digits, digits.I2, digits.pval,
                          scientific.pval, zero.pval, JAMA.pval,
                          big.mark, text.NA,
                          writexl)
  
  
  res <- list(fixed = fixed, random = random)
  ##
  res
}
