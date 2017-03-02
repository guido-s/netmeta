print.netrank <- function(x,
                          comb.fixed = x$x$comb.fixed,
                          comb.random = x$x$comb.random,
                          sort = TRUE,
                          digits = max(4, .Options$digits - 3),
                          ...) {
  
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  ##
  if (is.character(sort))
    sort <- meta:::setchar(sort, c("fixed", "random"))
  else
    meta:::chklogical(sort)
  ##
  meta:::chknumeric(digits, single = TRUE)
  
  
  both <- (comb.fixed + comb.random) == 2
  ##
  if (!both & is.character(sort)) {
    if (comb.fixed & sort == "random") {
      warning("Argument 'sort=\"random\"' ignored for P-scores of fixed effect model.",
              call. = FALSE)
      sort <- TRUE
    }
    if (comb.random & sort == "fixed") {
      warning("Argument 'sort=\"fixed\"' ignored for P-scores of random effects model.",
              call. = FALSE)
      sort <- TRUE
    }
  }
  else if (both & !is.character(sort) && sort)
    sort <- "random"
  
  
  if (both) {
    if (is.character(sort)) {
      res.both <- data.frame(fixed = round(x$Pscore.fixed, digits),
                             random = round(x$Pscore.random, digits))
      res.both <- res.both[order(-res.both[sort]), ]
    }
    else if (!sort) {
      res.both <- data.frame(fixed = round(x$Pscore.fixed[x$x$seq], digits),
                             random = round(x$Pscore.random[x$x$seq], digits))
    }
    ##
    colnames(res.both)  <- c("P-score (fixed)", "P-score (random)")
  }
  else {
    if (sort) {
      res.fixed <- as.data.frame(round(x$Pscore.fixed[order(-x$Pscore.fixed)],
                                       digits))
      res.random <- as.data.frame(round(x$Pscore.random[order(-x$Pscore.random)],
                                        digits))
    }
    else {
      res.fixed <- as.data.frame(round(x$Pscore.fixed[x$x$seq], digits))
      res.random <- as.data.frame(round(x$Pscore.random[x$x$seq], digits))
    }
    ##
    colnames(res.fixed)  <- "P-score"
    colnames(res.random) <- "P-score"
  }
  ##
  matitle(x)
  ##
  if (both)
    prmatrix(res.both, quote = FALSE, ...)
  ##
  else if (comb.fixed) {
    prmatrix(res.fixed, quote = FALSE, ...)
    if (comb.random)
      cat("\n")
  }
  ##
  else if (comb.random) {
    prmatrix(res.random, quote = FALSE, ...)
  }
  
  invisible(NULL)
}
