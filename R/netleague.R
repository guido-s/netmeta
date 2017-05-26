netleague <- function(x, y,
                      comb.fixed = x$comb.fixed, comb.random = x$comb.random,
                      seq = x$seq, ci = TRUE, backtransf = TRUE,
                      digits = gs("digits")) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  meta:::chkclass(x, "netmeta")
  ##
  if (!missing(y)) {
    meta:::chkclass(y, "netmeta")
    ##
    if (length(x$seq) != length(y$seq))
      stop("Arguments 'x' and 'y' must have the same number of treatments.")
    if (any(sort(x$seq) != sort(y$seq)))
      stop("Arguments 'x' and 'y' must have the same treatments.")
  }
  ##
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
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
  meta:::chklogical(ci)
  meta:::chklogical(backtransf)
  meta:::chknumeric(digits, min = 0, single = TRUE)
  
  
  ##
  ##
  ## (2) Back-transform log odds ratios & Co
  ##
  ##
  if (backtransf & meta:::is.relative.effect(x$sm)) {
    TE.fixed.x    <- exp(x$TE.fixed)
    lower.fixed.x <- exp(x$lower.fixed)
    upper.fixed.x <- exp(x$upper.fixed)
    ##
    TE.random.x    <- exp(x$TE.random)
    lower.random.x <- exp(x$lower.random)
    upper.random.x <- exp(x$upper.random)
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
  if (!missing(y)) {
    if (backtransf & meta:::is.relative.effect(y$sm)) {
      TE.fixed.y    <- exp(y$TE.fixed)
      lower.fixed.y <- exp(y$lower.fixed)
      upper.fixed.y <- exp(y$upper.fixed)
      ##
      TE.random.y    <- exp(y$TE.random)
      lower.random.y <- exp(y$lower.random)
      upper.random.y <- exp(y$upper.random)
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
  }
  
  
  ##
  ##
  ## (3) Print league table for fixed effect model
  ##
  ##
  if (comb.fixed) {
    cat("League table (fixed effect model):\n")
    TE.fixed.x    <- format(round(TE.fixed.x[seq.f, seq.f], digits))
    lower.fixed.x <- round(lower.fixed.x[seq.f, seq.f], digits)
    upper.fixed.x <- round(upper.fixed.x[seq.f, seq.f], digits)
    ##
    if (ci)
      nl.f <- paste(TE.fixed.x, meta:::p.ci(lower.fixed.x, upper.fixed.x))
    else
      nl.f <- TE.fixed.x
    ##
    nl.f <- matrix(nl.f, nrow = nrow(TE.fixed.x), ncol = ncol(TE.fixed.x))
    diag(nl.f) <- rownames(TE.fixed.x)
    ##
    if (!missing(y)) {
      TE.fixed.y    <- format(round(TE.fixed.y[seq.f, seq.f], digits))
      lower.fixed.y <- round(lower.fixed.y[seq.f, seq.f], digits)
      upper.fixed.y <- round(upper.fixed.y[seq.f, seq.f], digits)
      ##
      if (ci)
        nl.f.y <- paste(TE.fixed.y, meta:::p.ci(lower.fixed.y, upper.fixed.y))
      else
        nl.f.y <- TE.fixed.y
      ##
      nl.f.y <- matrix(nl.f.y, nrow = nrow(TE.fixed.y), ncol = ncol(TE.fixed.y))
      ##
      nl.f[upper.tri(nl.f)] <- nl.f.y[upper.tri(nl.f)]
    }
    ##
    prmatrix(nl.f, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(TE.fixed.x)),
             collab = rep("", ncol(TE.fixed.x)))
    if (comb.random)
      cat("\n")
  }
  
  
  ##
  ##
  ## (4) Print league table for random effects model
  ##
  ##
  if (comb.random) {
    cat("League table (random effects model):\n")
    TE.random.x    <- format(round(TE.random.x[seq.r, seq.r], digits))
    lower.random.x <- round(lower.random.x[seq.r, seq.r], digits)
    upper.random.x <- round(upper.random.x[seq.r, seq.r], digits)
    ##
    if (ci)
      nl.r <- paste(TE.random.x, meta:::p.ci(lower.random.x, upper.random.x))
    else
      nl.r <- TE.random.x
    ##
    nl.r <- matrix(nl.r, nrow = nrow(TE.random.x), ncol = ncol(TE.random.x))
    diag(nl.r) <- rownames(TE.random.x)
    ##
    if (!missing(y)) {
      TE.random.y    <- format(round(TE.random.y[seq.r, seq.r], digits))
      lower.random.y <- round(lower.random.y[seq.r, seq.r], digits)
      upper.random.y <- round(upper.random.y[seq.r, seq.r], digits)
      ##
      if (ci)
        nl.r.y <- paste(TE.random.y, meta:::p.ci(lower.random.y, upper.random.y))
      else
        nl.r.y <- TE.random.y
      ##
      nl.r.y <- matrix(nl.r.y, nrow = nrow(TE.random.y), ncol = ncol(TE.random.y))
      ##
      nl.r[upper.tri(nl.r)] <- nl.r.y[upper.tri(nl.r)]
    }
    ##
    prmatrix(nl.r, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(TE.random.x)),
             collab = rep("", ncol(TE.random.x)))
  }
  
  
  invisible(NULL)
}
