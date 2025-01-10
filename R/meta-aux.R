## Auxiliary functions
##
## Package: meta
## Author: Guido Schwarzer <guido.schwarzer@uniklinik-freiburg.de>
## License: GPL (>= 2)
##
bylevs <- function(x) {
  if (is.factor(x))
    res <- levels(factor(x))
  else
    res <- unique(x)
  res
}
byvarname <- function(argname, matchcall) {
  ##
  ## Determine name of subgroup variable
  ##
  res <- as.character(matchcall[[match(argname, names(matchcall))]])
  ##
  if (length(res) > 1 & res[1] == "$")
    res <- res[length(res)]
  ##
  if (length(res) == 0 || length(res) > 1)
    res <- "subgroup"
  ##
  res
}
catch <- function(argname, matchcall, data, encl) {
  ##
  ## Catch value for argument
  ##
  eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)
}
int2num <- function(x) {
  ##
  ## Convert integer to numeric
  ##
  if (is.integer(x))
    res <- as.numeric(x)
  else
    res <- x
  ##
  res
}
npn <- function(x) {
  ##
  ## Check for non-positive values in vector
  ##
  selNA <- is.na(x)
  res <- selNA
  if (sum(!selNA) > 0)
    x[!selNA] <- x[!selNA] <= 0
  ##
  res
}
replaceNULL <- function(x, replace = NA) {
  if (is.null(x))
    return(replace)
  x
}
replaceNA <- function(x, replace = NA) {
  if (is.null(x))
    return(x)
  else
    x[is.na(x)] <- replace
  x
}
warnarg <- function(x, y, fun, cl, otherarg) {
  if (x %in% y)
    if (!missing(cl))
      warning("Argument '", x, "' has been removed from R function ", fun,
              ".\nThis argument can be used in R function ", cl, ".",
              call. = FALSE)
    else if (!missing(otherarg))
      warning("Argument '", x, "' has been replaced by argument '", otherarg,
              "' in R function ", fun, ".\nSee help page of R function ",
              fun, " for information on the use of the new argument.",
              call. = FALSE)
  ##
  invisible(NULL)
}
catchvar <- function(varname, x, mf) {
  res <- NULL
  error <-
    try(res <- eval(mf[[match(varname, names(mf))]],
                    x,
                    enclos = sys.frame(sys.parent())),
        silent = TRUE)
  ##
  if (inherits(error, "try-error")) {
    res <- eval(mf[[match(varname, names(mf))]],
                x$data, enclos = NULL)
  }
  ##
  res
}
augment <- function(x, len, fun) {
  if (length(x) > 1)
    chklength(x, len, fun)
  else
    x <- rep(x, len)
  x
}
stoponly <- function(arg, val, func)
  stop("Argument ", arg, " =\"", val, "\"",
       " only defined for meta-analysis conducted with ",
       func, ".",
       call. = FALSE)
deprecated <- function(newvar, newmiss, args, old, warn = TRUE) {
  ##
  new <- deparse(substitute(newvar))
  ##
  if (length(args) == 0)
    return(newvar)
  ##
  if (is.list(args[[1]]))
    args <- args[[1]]
  ##
  additional.arguments <- names(args)
  ##
  if (!is.na(charmatch(old, additional.arguments)))
    if (!newmiss) {
      if (warn)
        warning("Deprecated argument '", old, "' ignored as ",
                "'", new, "' is also provided.",
                call. = FALSE)
      return(newvar)
    }
    else {
      if (warn)
        warning("Use argument '", new, "' instead of '",
                old, "' (deprecated).",
                call. = FALSE)
      return(args[[charmatch(old, additional.arguments)]])
    }
  else
    return(newvar)
}
deprecated2 <- function(newvar, newmiss, oldvar, oldmiss, warn = TRUE) {
  ##
  new <- deparse(substitute(newvar))
  old <- deparse(substitute(oldvar))
  ##
  if (newmiss & oldmiss)
    return(newvar)
  else if (!newmiss & oldmiss)
    return(newvar)
  else if (!newmiss & !oldmiss) {
    if (warn)
      warning("Deprecated argument '", old, "' ignored as ",
              "'", new, "' is also provided.",
              call. = FALSE)
    return(newvar)
  }
  else if (newmiss & !oldmiss) {
    if (warn)
      warning("Use argument '", new, "' instead of '",
              old, "' (deprecated).",
              call. = FALSE)
    return(oldvar)
  }
}
runNN <- function(func, args, warn = TRUE) {
  args <- args[!sapply(args, is.null)]
  if (warn)
    do.call(func, args)
  else
    suppressWarnings(do.call(func, args))
}

# Estimate the heterogeneity parameter phi using the
# modified version of Pearson's statistic.
#
phi <- function(x) {
  # Extract number of trials
  n.trials <- x$prior.weights
  #
  if (identical(unique(n.trials), 1))
    stop("The number of successes must be summarized for valid computation of ",
         "c-hat.")
  
  # Pearson chi-square
  chisq <- sum(residuals(x, type = "pearson")^2)
  
  # Extract raw residuals
  raw.res <- residuals(x, type = "response")
  
  # Extract fitted values
  fit.vals <- fitted(x)
  
  # Estimate s.bar
  s.bar <- mean((1 - 2 * fit.vals) / ((n.trials * fit.vals) * (1 - fit.vals)))
  
  # Calculate estimate based on Fletcher estimator
  phi <- (chisq / x$df.residual) / (1 + s.bar)
  
  # Set phi = 1 if phi < 1 to remain consistent with common effect model
  phi <- max(phi, 1)
  
  phi
}

# Gets as arguments
# - 'trts' = treatment names,
# - 'TE.basic' = basic treatment effect estimates,
# - 'vcov.basic' = variance-covariance matrix of the estimates.
#
# Returns all NMA estimates (basic + functional) and their standard errors.

basic2all <- function(trts, TE.basic, vcov.basic) {
  
  n.trts <- length(trts)
  
  # Get rid of warning 'no visible binding for global variable'
  treat1 <- treat2 <- NULL
  
  # Identify all comparisons
  #
  comps <- expand.grid(treat1 = seq_len(n.trts), treat2 = seq_len(n.trts))
  comps <- subset(comps, treat2 < treat1)
  #
  comps.basic <- subset(comps, treat2 == 1)
  comps.other <- subset(comps, treat2 > 1)
  #
  n.comps.basic <- nrow(comps.basic)
  n.comps.other <- nrow(comps.other)
  
  # Construct X matrix, i.e., the treatment comparison matrix
  #
  X.basic <- matrix(0, nrow = n.comps.other, ncol = n.comps.basic)
  #
  for (i in seq_len(n.comps.other)) {
    for (j in seq_len(n.comps.basic)) {
      if (comps.other$treat1[i] == comps.basic$treat1[j])
        X.basic[i, j] <- 1
      #
      if (comps.other$treat2[i] == comps.basic$treat1[j])
        X.basic[i, j] <- -1
    }
  }
  
  dat.basic <-
    data.frame(comps.basic,
               TE = TE.basic,
               seTE = sqrt(diag(vcov.basic)))
  #
  dat.other <-
    data.frame(comps.other,
               TE = X.basic %*% TE.basic,
               seTE = sqrt(diag(X.basic %*% vcov.basic %*% t(X.basic))))
  #
  dat.comps <- rbind(dat.basic, dat.other)
  
  TEs <- matrix(NA, ncol = n.trts, nrow = n.trts)
  seTEs <- matrix(NA, ncol = n.trts, nrow = n.trts)
  #
  for (i in seq_len(nrow(dat.comps))) {
    TEs[dat.comps$treat1[i], dat.comps$treat2[i]] <-  dat.comps$TE[i]
    TEs[dat.comps$treat2[i], dat.comps$treat1[i]] <- -dat.comps$TE[i]
    #
    seTEs[dat.comps$treat1[i], dat.comps$treat2[i]] <- dat.comps$seTE[i]
    seTEs[dat.comps$treat2[i], dat.comps$treat1[i]] <- dat.comps$seTE[i]
  }
  
  diag(TEs) <- diag(seTEs) <- 0
  #
  rownames(TEs) <- rownames(seTEs) <- colnames(TEs) <- colnames(seTEs) <- trts
  
  res <- list(TEs = TEs, seTEs = seTEs)
  #
  res
}
