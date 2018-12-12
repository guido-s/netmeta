print.netsplit <- function(x,
                           comb.fixed = x$comb.fixed,
                           comb.random = x$comb.random,
                           show = "all",
                           overall = TRUE,
                           ci = FALSE,
                           test = show %in% c("all", "with.direct", "both"),
                           digits = gs("digits"),
                           digits.zval = gs("digits.zval"),
                           digits.pval = gs("digits.pval"),
                           digits.prop = max(gs("digits.pval") - 2, 2),
                           text.NA = ".",
                           backtransf = x$backtransf,
                           scientific.pval = gs("scientific.pval"),
                           big.mark = gs("big.mark"),
                           legend = TRUE,
                           ...) {


  meta:::chkclass(x, "netsplit")
  ##
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  formatCI <- meta:::formatCI
  formatN <- meta:::formatN
  formatPT <- meta:::formatPT
  is.relative.effect <- meta:::is.relative.effect
  rmSpace <- meta:::rmSpace
  setchar <- meta:::setchar


  ## All individual results in a single row - be on the save side:
  ##
  oldopts <- options(width = 200)
  on.exit(options(oldopts))


  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(overall)
  chklogical(ci)
  chklogical(test)
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.zval, min = 0, single = TRUE)
  chknumeric(digits.pval, min = 1, single = TRUE)
  chknumeric(digits.prop, min = 0, single = TRUE)
  ##
  if (is.null(backtransf))
    backtransf <- TRUE
  chklogical(backtransf)
  chklogical(scientific.pval)
  chklogical(legend)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  ## Check whether first argument is a list. In this case only use
  ## this list as input.
  if (length(args) > 0 && is.list(args[[1]]))
    args <- args[[1]]
  ##
  additional.arguments <- names(args)
  ##
  if (length(additional.arguments) > 0) {
    if (!is.na(charmatch("showa", additional.arguments)))
      if (!missing(show))
        warning("Deprecated argument 'showall' ignored as argument 'show' is also provided.")
      else {
        warning("Deprecated argument 'showall' has been replaced by argument 'show'.")
        show <- args[[charmatch("showa", additional.arguments)]]
        if (show)
          show <- "all"
        else
          show <- "both"
      }
  }
  ##
  show <- setchar(show, c("all", "both", "with.direct", "direct.only", "indirect.only"))


  sm <- x$sm
  sm.lab <- sm
  ##
  relative <- is.relative.effect(sm)
  ##
  if (!backtransf & relative)
    sm.lab <- paste("log", sm, sep = "")
  ##
  if (!(sm.lab == "" | sm.lab == "log"))
    sm.lab <- paste("(", sm.lab, ") ", sep = "")
  else
    sm.lab <- ""


  level.comb <- x$level.comb
  ci.lab <- paste(100 * level.comb, "%-CI", sep ="")
  
  
  random.available <- !is.null(x$random)
  ##
  if (!random.available & comb.random) {
    warning("No results for random effects model available. ",
            "Argument 'comb.random' set to FALSE.",
            call. = FALSE)
    ##
    comb.random <- FALSE
  }
  
  
  if (show == "all")
    sel <- rep_len(TRUE, length(x$direct.fixed$TE))
  else if (show == "with.direct")
    sel <- (!is.na(x$direct.fixed$TE) & !is.na(x$direct.random$TE))
  else if (show == "both")
    sel <- (!is.na(x$direct.fixed$TE)  & !is.na(x$indirect.fixed$TE) &
            !is.na(x$direct.random$TE) & !is.na(x$indirect.random$TE))
  else if (show == "direct.only")
    sel <- (!is.na(x$direct.fixed$TE)  & is.na(x$indirect.fixed$TE) &
            !is.na(x$direct.random$TE) & is.na(x$indirect.random$TE))
  else if (show == "indirect.only")
    sel <- (is.na(x$direct.fixed$TE)  & !is.na(x$indirect.fixed$TE) &
            is.na(x$direct.random$TE) & !is.na(x$indirect.random$TE))
  ##
  comp <- x$comparison[sel]
  ##
  k <- x$k[sel]
  ##
  prop.fixed <- x$prop.fixed[sel]
  ##
  TE.fixed <- x$fixed$TE[sel]
  lower.fixed <- x$fixed$lower[sel]
  upper.fixed <- x$fixed$upper[sel]
  ##
  TE.direct.fixed <- x$direct.fixed$TE[sel]
  lower.direct.fixed <- x$direct.fixed$lower[sel]
  upper.direct.fixed <- x$direct.fixed$upper[sel]
  ##
  TE.indirect.fixed <- x$indirect.fixed$TE[sel]
  lower.indirect.fixed <- x$indirect.fixed$lower[sel]
  upper.indirect.fixed <- x$indirect.fixed$upper[sel]
  ##
  TE.compare.fixed <- x$compare.fixed$TE[sel]
  lower.compare.fixed <- x$compare.fixed$lower[sel]
  upper.compare.fixed <- x$compare.fixed$upper[sel]
  zval.compare.fixed <- x$compare.fixed$z[sel]
  pval.compare.fixed <- x$compare.fixed$p[sel]
  ##
  if (random.available) {
    prop.random <- x$prop.random[sel]
    ##
    TE.random <- x$random$TE[sel]
    lower.random <- x$random$lower[sel]
    upper.random <- x$random$upper[sel]
    ##
    TE.direct.random <- x$direct.random$TE[sel]
    lower.direct.random <- x$direct.random$lower[sel]
    upper.direct.random <- x$direct.random$upper[sel]
    ##
    TE.indirect.random <- x$indirect.random$TE[sel]
    lower.indirect.random <- x$indirect.random$lower[sel]
    upper.indirect.random <- x$indirect.random$upper[sel]
    ##
    TE.compare.random <- x$compare.random$TE[sel]
    lower.compare.random <- x$compare.random$lower[sel]
    upper.compare.random <- x$compare.random$upper[sel]
    zval.compare.random <- x$compare.random$z[sel]
    pval.compare.random <- x$compare.random$p[sel]
  }
  
  
  if (backtransf & relative) {
    TE.fixed <- exp(TE.fixed)
    lower.fixed <- exp(lower.fixed)
    upper.fixed <- exp(upper.fixed)
    ##
    TE.direct.fixed <- exp(TE.direct.fixed)
    lower.direct.fixed <- exp(lower.direct.fixed)
    upper.direct.fixed <- exp(upper.direct.fixed)
    ##
    TE.indirect.fixed <- exp(TE.indirect.fixed)
    lower.indirect.fixed <- exp(lower.indirect.fixed)
    upper.indirect.fixed <- exp(upper.indirect.fixed)
    ##
    TE.compare.fixed <- exp(TE.compare.fixed)
    lower.compare.fixed <- exp(lower.compare.fixed)
    upper.compare.fixed <- exp(upper.compare.fixed)
    ##
    if (random.available) {
      TE.random <- exp(TE.random)
      lower.random <- exp(lower.random)
      upper.random <- exp(upper.random)
      ##
      TE.direct.random <- exp(TE.direct.random)
      lower.direct.random <- exp(lower.direct.random)
      upper.direct.random <- exp(upper.direct.random)
      ##
      TE.indirect.random <- exp(TE.indirect.random)
      lower.indirect.random <- exp(lower.indirect.random)
      upper.indirect.random <- exp(upper.indirect.random)
      ##
      TE.compare.random <- exp(TE.compare.random)
      lower.compare.random <- exp(lower.compare.random)
      upper.compare.random <- exp(upper.compare.random)
    }
  }
  
  
  fixed <- list(comp = comp,
                k = k,
                prop = formatPT(prop.fixed, digits = digits.prop))
  names.fixed <- c("comparison", "k", "prop")
  ##
  if (overall) {
    fixed$TE.fixed <- formatN(TE.fixed, digits, text.NA = text.NA,
                              big.mark = big.mark)
    names.fixed <- c(names.fixed, "nma")
    if (ci) {
      fixed$ci.fixed <- formatCI(round(lower.fixed, digits),
                                 round(upper.fixed, digits))
      fixed$ci.fixed[is.na(fixed$ci.fixed)] <- text.NA
      names.fixed <- c(names.fixed, ci.lab)
    }
  }
  ##
  fixed$TE.direct.fixed <- formatN(TE.direct.fixed, digits, text.NA = text.NA,
                                   big.mark = big.mark)
  names.fixed <- c(names.fixed, "direct")
  if (ci) {
    fixed$ci.direct.fixed <- formatCI(round(lower.direct.fixed, digits),
                                      round(upper.direct.fixed, digits))
    fixed$ci.direct.fixed[is.na(fixed$ci.direct.fixed)] <- text.NA
    names.fixed <- c(names.fixed, ci.lab)
  }
  ##
  fixed$TE.indirect.fixed <- formatN(TE.indirect.fixed, digits,
                                     text.NA = text.NA, big.mark = big.mark)
  names.fixed <- c(names.fixed, "indir.")
  ##
  if (ci) {
    fixed$ci.indirect.fixed <- formatCI(round(lower.indirect.fixed, digits),
                                        round(upper.indirect.fixed, digits))
    fixed$ci.indirect.fixed[is.na(fixed$ci.indirect.fixed)] <- text.NA
    names.fixed <- c(names.fixed, ci.lab)
  }
  ##
  if (test) {
    fixed$diff <- formatN(TE.compare.fixed, digits, text.NA = text.NA,
                          big.mark = big.mark)
    names.fixed <- c(names.fixed, if (backtransf & relative) "RoR" else "Diff")
    if (ci) {
      fixed$ci.diff <- formatCI(round(lower.compare.fixed, digits),
                                round(upper.compare.fixed, digits))
      fixed$ci.diff[is.na(fixed$ci.diff)] <- text.NA
      names.fixed <- c(names.fixed, ci.lab)
    }
    ##
    fixed$z <- formatN(zval.compare.fixed, digits.zval,
                       big.mark = big.mark)
    fixed$z[fixed$z == "--"] <- text.NA
    fixed$p <- formatPT(pval.compare.fixed, digits = digits.pval,
                        scientific = scientific.pval)
    fixed$p[rmSpace(fixed$p) == "--"] <- text.NA
    names.fixed <- c(names.fixed, c("z", "p-value"))
  }
  fixed <- as.data.frame(fixed)
  names(fixed) <- names.fixed
  ##
  if (is.null(prop.fixed))
    fixed <- fixed[, !(names(fixed) %in% "prop")]
  
  
  if (random.available) {
    random <- list(comp = comp,
                   k = k,
                   prop = formatPT(prop.random, digits = digits.prop))
    names.random <- c("comparison", "k", "prop")
    ##
    if (overall) {
      random$TE.random <- formatN(TE.random, digits, text.NA = text.NA,
                                  big.mark = big.mark)
      names.random <- c(names.random, "nma")
      if (ci) {
        random$ci.random <- formatCI(round(lower.random, digits),
                                     round(upper.random, digits))
        random$ci.random[is.na(random$ci.random)] <- text.NA
        names.random <- c(names.random, ci.lab)
      }
    }
    ##
    random$TE.direct.random <- formatN(TE.direct.random, digits,
                                       text.NA = text.NA,
                                       big.mark = big.mark)
    names.random <- c(names.random, "direct")
    if (ci) {
      random$ci.direct.random <- formatCI(round(lower.direct.random, digits),
                                          round(upper.direct.random, digits))
      random$ci.direct.random[is.na(random$ci.direct.random)] <- text.NA
      names.random <- c(names.random, ci.lab)
    }
    ##
    random$TE.indirect.random <- formatN(TE.indirect.random, digits,
                                         text.NA = text.NA,
                                         big.mark = big.mark)
    names.random <- c(names.random, "indir.")
    if (ci) {
      random$ci.indirect.random <- formatCI(round(lower.indirect.random, digits),
                                            round(upper.indirect.random, digits))
      random$ci.indirect.random[is.na(random$ci.indirect.random)] <- text.NA
      names.random <- c(names.random, ci.lab)
    }
    ##
    if (test) {
      random$diff <- formatN(TE.compare.random, digits, text.NA = text.NA,
                             big.mark = big.mark)
      names.random <- c(names.random, if (backtransf & relative) "RoR" else "Diff")
      if (ci) {
        random$ci.diff <- formatCI(round(lower.compare.random, digits),
                                   round(upper.compare.random, digits))
        random$ci.diff[is.na(random$ci.diff)] <- text.NA
        names.random <- c(names.random, ci.lab)
      }
      ##
      random$z <- formatN(zval.compare.random, digits.zval,
                          big.mark = big.mark)
      random$z[random$z == "--"] <- text.NA
      random$p <- formatPT(pval.compare.random, digits = digits.pval,
                           scientific = scientific.pval)
      random$p[rmSpace(random$p) == "--"] <- text.NA
      names.random <- c(names.random, c("z", "p-value"))
    }
    random <- as.data.frame(random)
    names(random) <- names.random
  }
  
  
  cat(x$method, "method to split direct and indirect evidence\n\n")
  
  
  if (comb.fixed) {
    cat("Fixed effect model: \n\n")
    fixed[is.na(fixed)] <- text.NA
    prmatrix(fixed, quote = FALSE, right = TRUE,
             rowlab = rep("", dim(fixed)[1]))
    if (comb.random)
      cat("\n")
  }
  ##
  if (comb.random) {
    cat("Random effects model: \n\n")
    random[is.na(random)] <- text.NA
    prmatrix(random, quote = FALSE, right = TRUE,
             rowlab = rep("", dim(random)[1]))
  }
  ##
  if (legend) {
    cat("\nLegend:\n")
    cat(" comparison - Treatment comparison\n")
    cat(" k          - Number of studies providing direct evidence\n")
    if (!is.null(prop.fixed))
      cat(" prop       - Direct evidence proportion\n")
    if (overall)
      cat(paste(" nma        - Estimated treatment effect ", sm.lab,
                "in network meta-analysis\n", sep = ""))
    cat(paste(" direct     - Estimated treatment effect ", sm.lab,
              "derived from direct evidence\n", sep = ""))
    cat(paste(" indir.     - Estimated treatment effect ", sm.lab,
              "derived from indirect evidence\n", sep = ""))
    if (test) {
      if (backtransf & relative)
        cat(" RoR        - Ratio of Ratios (direct versus indirect)\n")
      else
        cat(" Diff       - Difference between direct and indirect treatment estimates\n")
      cat(" z          - z-value of test for disagreement (direct versus indirect)\n")
      cat(" p-value    - p-value of test for disagreement (direct versus indirect)\n")
    }
  }


  invisible(NULL)
}
