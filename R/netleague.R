netleague <- function(x, y,
                      comb.fixed = x$comb.fixed, comb.random = x$comb.random,
                      seq = x$seq, ci = TRUE, backtransf = TRUE,
                      direct = FALSE,
                      digits = gs("digits"),
                      bracket = gs("CIbracket"),
                      separator = gs("CIseparator"),
                      text.NA = ".",
                      big.mark = gs("big.mark")) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  meta:::chkclass(x, "netmeta")
  ##
  chkchar <- meta:::chkchar
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  formatCI <- meta:::formatCI
  formatN <- meta:::formatN
  is.relative.effect <- meta:::is.relative.effect
  ##
  if (!missing(y)) {
    meta:::chkclass(y, "netmeta")
    ##
    if (length(x$seq) != length(y$seq))
      stop("Arguments 'x' and 'y' must have the same number of treatments.")
    if (any(sort(x$seq) != sort(y$seq)))
      stop("Arguments 'x' and 'y' must have the same treatments.")
    ##
    x.is.y <- all.equal(x, y)
    x.is.y <- is.logical(x.is.y) && x.is.y
  }
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  ##
  if (missing(seq)) {
    seq.f <- seq.r <- seq
  }
  else {
    if (is.null(seq))
      stop("Argument 'seq' must be not NULL.")
    else if (inherits(seq, "netrank")) {
      pscore.f <- seq$Pscore.fixed
      pscore.r <- seq$Pscore.random
      ##
      seq.f <- setseq(names(pscore.f)[rev(order(pscore.f))], x$seq)
      seq.r <- setseq(names(pscore.r)[rev(order(pscore.r))], x$seq)
    }
    else
      seq.f <- seq.r <- setseq(seq, x$seq)
  }
  ##
  chklogical(ci)
  chklogical(backtransf)
  chklogical(direct)
  chknumeric(digits, min = 0, single = TRUE)
  ##
  bracket.old <- gs("CIbracket")
  separator.old <- gs("CIseparator")
  cilayout(bracket, separator)
  on.exit(cilayout(bracket.old, separator.old))
  ##
  chkchar(text.NA)
  chkchar(big.mark)
  
  
  ##
  ##
  ## (2) Back-transform log odds ratios & Co
  ##
  ##
  if (!missing(y) & direct) {
    TE.fixed.x    <- x$TE.direct.fixed
    lower.fixed.x <- x$lower.direct.fixed
    upper.fixed.x <- x$upper.direct.fixed
    ##
    TE.random.x    <- x$TE.direct.random
    lower.random.x <- x$lower.direct.random
    upper.random.x <- x$upper.direct.random
  }
  else {
    TE.fixed.x    <- x$TE.fixed
    lower.fixed.x <- x$lower.fixed
    upper.fixed.x <- x$upper.fixed
    ##
    TE.random.x    <- x$TE.random
    lower.random.x <- x$lower.random
    upper.random.x <- x$upper.random
  }
  ##
  if (backtransf & is.relative.effect(x$sm)) {
    TE.fixed.x    <- exp(TE.fixed.x)
    lower.fixed.x <- exp(lower.fixed.x)
    upper.fixed.x <- exp(upper.fixed.x)
    ##
    TE.random.x    <- exp(TE.random.x)
    lower.random.x <- exp(lower.random.x)
    upper.random.x <- exp(upper.random.x)
  }
  ##
  if (!missing(y)) {
    if (direct) {
      TE.fixed.y    <- y$TE.direct.fixed
      lower.fixed.y <- y$lower.direct.fixed
      upper.fixed.y <- y$upper.direct.fixed
      ##
      TE.random.y    <- y$TE.direct.random
      lower.random.y <- y$lower.direct.random
      upper.random.y <- y$upper.direct.random
    }
    else {
      TE.fixed.y    <- y$TE.fixed
      lower.fixed.y <- y$lower.fixed
      upper.fixed.y <- y$upper.fixed
      ##
      TE.random.y    <- y$TE.random
      lower.random.y <- y$lower.random
      upper.random.y <- y$upper.random
    }
    ##
    if (backtransf & is.relative.effect(y$sm)) {
      TE.fixed.y    <- exp(TE.fixed.y)
      lower.fixed.y <- exp(lower.fixed.y)
      upper.fixed.y <- exp(upper.fixed.y)
      ##
      TE.random.y    <- exp(TE.random.y)
      lower.random.y <- exp(lower.random.y)
      upper.random.y <- exp(upper.random.y)
    }
    ##
    if (x.is.y) {
      TE.fixed.y    <- t(TE.fixed.y)
      lower.fixed.y <- t(lower.fixed.y)
      upper.fixed.y <- t(upper.fixed.y)
      ##
      TE.random.y    <- t(TE.random.y)
      lower.random.y <- t(lower.random.y)
      upper.random.y <- t(upper.random.y)
    }
  }
  else {
    if (backtransf & is.relative.effect(x$sm)) {
      TE.fixed.y    <- exp(x$TE.direct.fixed)
      lower.fixed.y <- exp(x$lower.direct.fixed)
      upper.fixed.y <- exp(x$upper.direct.fixed)
      ##
      TE.random.y    <- exp(x$TE.direct.random)
      lower.random.y <- exp(x$lower.direct.random)
      upper.random.y <- exp(x$upper.direct.random)
    }
    else {
      TE.fixed.y    <- x$TE.direct.fixed
      lower.fixed.y <- x$lower.direct.fixed
      upper.fixed.y <- x$upper.direct.fixed
      ##
      TE.random.y    <- x$TE.direct.random
      lower.random.y <- x$lower.direct.random
      upper.random.y <- x$upper.direct.random
    }
  }
  ##
  ## Comparisons are column versus row
  ##
  TE.fixed.x <- t(TE.fixed.x)
  lower.fixed.x <- t(lower.fixed.x)
  upper.fixed.x <- t(upper.fixed.x)
  ##
  TE.random.x <- t(TE.random.x)
  lower.random.x <- t(lower.random.x)
  upper.random.x <- t(upper.random.x)
  ##
  TE.fixed.y <- t(TE.fixed.y)
  lower.fixed.y <- t(lower.fixed.y)
  upper.fixed.y <- t(upper.fixed.y)
  ##
  TE.random.y <- t(TE.random.y)
  lower.random.y <- t(lower.random.y)
  upper.random.y <- t(upper.random.y)
  
  
  ##
  ##
  ## (3) Print league table for fixed effect model
  ##
  ##
  TE.fixed.x    <- round(   TE.fixed.x[seq.f, seq.f], digits)
  lower.fixed.x <- round(lower.fixed.x[seq.f, seq.f], digits)
  upper.fixed.x <- round(upper.fixed.x[seq.f, seq.f], digits)
  ##
  if (ci) {
    nl.NA <- is.na(TE.fixed.x)
    nl.f <- paste(formatN(TE.fixed.x, text.NA = text.NA, big.mark = big.mark),
                  formatCI(lower.fixed.x, upper.fixed.x, lab.NA = text.NA,
                           big.mark = big.mark))
    nl.f[nl.NA] <- text.NA
  }
  else
    nl.f <- formatN(TE.fixed.x, text.NA = text.NA, big.mark = big.mark)
  ##
  nl.f <- matrix(nl.f, nrow = nrow(TE.fixed.x), ncol = ncol(TE.fixed.x))
  diag(nl.f) <- rownames(TE.fixed.x)
  ##
  TE.fixed.y    <- round(   TE.fixed.y[seq.f, seq.f], digits)
  lower.fixed.y <- round(lower.fixed.y[seq.f, seq.f], digits)
  upper.fixed.y <- round(upper.fixed.y[seq.f, seq.f], digits)
  ##
  if (ci) {
    nl.NA <- is.na(TE.fixed.y)
    nl.f.y <- paste(formatN(TE.fixed.y, text.NA = text.NA, big.mark = big.mark),
                    formatCI(lower.fixed.y, upper.fixed.y, lab.NA = text.NA,
                             big.mark = big.mark))
    nl.f.y[nl.NA] <- text.NA
  }
  else
    nl.f.y <- formatN(TE.fixed.y, text.NA = text.NA, big.mark = big.mark)
  ##
  nl.f.y <- matrix(nl.f.y, nrow = nrow(TE.fixed.y), ncol = ncol(TE.fixed.y))
  ##
  nl.f[upper.tri(nl.f)] <- t(nl.f.y)[upper.tri(nl.f)]
  ##
  nl.f <- as.data.frame(nl.f, stringsAsFactors = FALSE)
  
  
  ##
  ##
  ## (4) Print league table for random effects model
  ##
  ##
  TE.random.x    <- round(   TE.random.x[seq.r, seq.r], digits)
  lower.random.x <- round(lower.random.x[seq.r, seq.r], digits)
  upper.random.x <- round(upper.random.x[seq.r, seq.r], digits)
  ##
  if (ci) {
    nl.NA <- is.na(TE.random.x)
    nl.r <- paste(formatN(TE.random.x, text.NA = text.NA, big.mark = big.mark),
                  formatCI(lower.random.x, upper.random.x, lab.NA = text.NA,
                           big.mark = big.mark))
    nl.r[nl.NA] <- text.NA
  }
  else
    nl.r <- formatN(TE.random.x, text.NA = text.NA, big.mark = big.mark)
  ##
  nl.r <- matrix(nl.r, nrow = nrow(TE.random.x), ncol = ncol(TE.random.x))
  diag(nl.r) <- rownames(TE.random.x)
  ##
  TE.random.y    <- round(   TE.random.y[seq.r, seq.r], digits)
  lower.random.y <- round(lower.random.y[seq.r, seq.r], digits)
  upper.random.y <- round(upper.random.y[seq.r, seq.r], digits)
  ##
  if (ci) {
    nl.NA <- is.na(TE.random.y)
    nl.r.y <- paste(formatN(TE.random.y, text.NA = text.NA),
                    formatCI(lower.random.y, upper.random.y, lab.NA = text.NA))
    nl.r.y[nl.NA] <- text.NA
  }
  else
    nl.r.y <- formatN(TE.random.y, text.NA = text.NA, big.mark = big.mark)
  ##
  nl.r.y <- matrix(nl.r.y, nrow = nrow(TE.random.y), ncol = ncol(TE.random.y))
  ##
  nl.r[upper.tri(nl.r)] <- t(nl.r.y)[upper.tri(nl.f)]
  ##
  nl.r <- as.data.frame(nl.r, stringsAsFactors = FALSE)
  
  
  res <- list(fixed = nl.f,
              random = nl.r,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              seq = seq, ci = ci, backtransf = backtransf,
              digits = digits)
  ##
  class(res) <- "netleague"
  
  
  res
}
