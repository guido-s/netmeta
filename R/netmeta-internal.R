.onAttach <- function (libname, pkgname) {
  msg <-
    paste0("Loading 'netmeta' package (version ",
           packageDescription("netmeta")$Version,
           ").",
           "\nType 'help(\"netmeta-package\")' for a brief overview.")
  packageStartupMessage(msg)
}


.onLoad <- function(libname, pkgname) {
  #
  # Default settings for R package netmeta
  #
  argslist <- c("baseline.reference", "small.values",
                "all.treatments", "seq",
                "method.tau.netmeta",
                "drop.reference.group", "equal.size",
                "show",
                #
                "nsim", "lump.comparator",
                #
                "plastic", "col.netgraph",
                "number.of.studies", "thickness",
                "multiarm",
                #
                "tol.multiarm", "tol.multiarm.se",
                "details.chkmultiarm",
                #
                "na.unident",
                "sep.trts", "sep.comps", "sep.ia",
                "nchar.trts", "nchar.studlab",
                #
                "sort.distance",
                #
                "legend")
  #
  suppressWarnings(
    settings.meta(baseline.reference = TRUE, small.values = "desirable",
                  all.treatments = NULL, seq = NULL,
                  method.tau.netmeta = "DL",
                  drop.reference.group = TRUE, equal.size = TRUE,
                  show = "both",
                  #
                  nsim = 1000, lump.comparator = FALSE,
                  #
                  plastic = FALSE, col.netgraph = NULL,
                  number.of.studies = TRUE, thickness = "number.of.studies",
                  multiarm = FALSE,
                  #
                  tol.multiarm = 0.001, tol.multiarm.se = NULL,
                  details.chkmultiarm = FALSE,
                  #
                  na.unident = TRUE,
                  sep.trts = ":", sep.comps = "+", sep.ia = "*",
                  nchar.trts = 666, nchar.studlab = 666,
                  #
                  sort.distance = TRUE,
                  #
                  legend = TRUE,
                  #
                  .argslist.netmeta = argslist))
}


calcV <- function(x, sm) {
  p2 <- (x$event2[1] + x$incr[1]) / (x$n2[1] + 2 * x$incr[1])
  n2 <- x$n2[1] + 2 * x$incr[1]
  ##
  if (sm == "OR")
    V <- matrix(1 / (x$event2[1] + x$incr[1]) +
                1 / (x$n2[1] - x$event2[1] + x$incr[1]),
                nrow = nrow(x), ncol = nrow(x))
  else if (sm == "RR")
    V <- matrix((1 - p2) / (n2 * p2),
                nrow = nrow(x), ncol = nrow(x))
  else if (sm == "RD")
    V <- matrix(p2 * (1 - p2) / n2,
                nrow = nrow(x), ncol = nrow(x))
  else if (sm == "ASD")
    V <- matrix(0.25 / x$n2[1],
                nrow = nrow(x), ncol = nrow(x))
  else if (sm == "MD")
    V <- matrix(x$sd2[1]^2 / x$n2[1],
                nrow = nrow(x), ncol = nrow(x))
  else if (sm == "IRR")
    V <- matrix(1 / (x$event2[1] + x$incr[1]),
                nrow = nrow(x), ncol = nrow(x))
  else if (sm == "IRD")
    V <- matrix((x$event2[1] + x$incr[1]) / x$time2[1]^2,
                nrow = nrow(x), ncol = nrow(x))
  else if (sm == "IRSD")
    V <- matrix(0.25 / x$time2[1],
                nrow = nrow(x), ncol = nrow(x))
  ##
  diag(V) <- x$seTE^2
  ##
  V
}


##
## Abbreviated component labels
##
compos <- function(x, lev, abbr, split, add) {
  x.list <- compsplit(x, split = split)
  if (!missing(add))
    add <- ifelse(add, " ", "")
  else
    add <- ifelse(attr(x.list, "withspace"), " ", "")
  ##
  x.list <- lapply(x.list, charfac, levels = lev, labels = abbr)
  x.list <- lapply(x.list, paste, collapse = paste0(add, split, add))
  ##
  unlist(x.list)
}


charfac <- function(x, ...)
  as.character(factor(x, ...))


legendabbr <- function(full, abbr, condition,
                       text = "Treatment name",
                       header = "\nLegend:\n") {
  diff <- full != abbr
  if (condition) {
    if (any(diff)) {
      tmat <- data.frame(abbr, full)
      names(tmat) <- c("Abbreviation", text)
      tmat <- tmat[diff, ]
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat(header)
      ##
      prmatrix(unique(tmat), quote = FALSE, right = TRUE,
               rowlab = rep("", nrow(unique(tmat))))
    }
  }
  invisible(condition & any(diff))
}


setsv <- function(x, add = NULL) {
  if (is.null(x))
    res <- "desirable"
  else {
    res <- setchar(x, c("good", "bad"), stop.at.error = FALSE)
    ##
    if (!is.null(res))
      res <- switch(res, good = "desirable", bad = "undesirable")
    else
      res <- x
  }
  ##
  setchar(res, c("desirable", "undesirable", add))
}

first <- function(x) x[1]
second <- function(x) x[2]


# For forest plots
#
setHet <- function(meta, netmeta) {
  #
  meta$method.tau <- netmeta$method.tau
  meta$tau2 <- netmeta$tau2
  meta$se.tau2 <- NA
  meta$lower.tau2 <- NA
  meta$upper.tau2 <- NA
  #
  meta$tau <- netmeta$tau
  meta$lower.tau <- NA
  meta$upper.tau <- NA
  #
  meta$method.tau.ci <- ""
  meta$sign.lower.tau <- ""
  meta$sign.upper.tau <- ""
  #
  meta$Q <- netmeta$Q
  meta$df.Q <- netmeta$df.Q
  meta$pval.Q <- netmeta$pval.Q
  #
  meta$method.I2 <- "Q"
  meta$I2 <- netmeta$I2
  meta$lower.I2 <- netmeta$.lower.I2
  meta$upper.I2 <- netmeta$.upper.I2
  #
  meta$H <- netmeta$H
  meta$lower.H <- netmeta$.lower.H
  meta$upper.H <- netmeta$.upper.H
  #
  meta$tau2.resid <- NA
  meta$se.tau2.resid <- NA
  meta$lower.tau2.resid <- NA
  meta$upper.tau2.resid <- NA
  #
  meta$tau.resid <- NA
  meta$lower.tau.resid <- NA
  meta$upper.tau.resid <- NA
  #
  meta$Q.resid <- netmeta$Q.resid
  meta$df.Q.resid <- netmeta$df.Q.resid
  meta$pval.Q.resid <- netmeta$pval.Q.resid
  #
  meta$H.resid <- netmeta$H.resid
  meta$lower.H.resid <- netmeta$.lower.H.resid
  meta$upper.H.resid <- netmeta$.upper.H.resid
  #
  meta$I2.resid <- netmeta$I2.resid
  meta$lower.I2.resid <- netmeta$.lower.I2.resid
  meta$upper.I2.resid <- netmeta$.upper.I2.resid
  #
  meta
}


is_identical <- function(x) {
  duplicated(as.data.frame(t(x))) |
    duplicated(as.data.frame(t(x)), fromLast = TRUE)
}
