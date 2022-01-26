.onAttach <- function (libname, pkgname) {
  msg <-
    paste0("Loading 'netmeta' package (version ",
           packageDescription("netmeta")$Version,
           ").",
           "\nType 'help(\"netmeta-package\")' for a brief overview.",
           "\nReaders of 'Meta-Analysis with R (Use R!)' should install",
           "\nolder version of 'netmeta' package: ",
           "https://tinyurl.com/kyz6wjbb")
  packageStartupMessage(msg)
}


.special.characters <- c("+", ".", "&", "$", "#", "|", "*", "^")


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
