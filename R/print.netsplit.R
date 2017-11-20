print.netsplit <- function(x,
                           comb.fixed = x$comb.fixed,
                           comb.random = x$comb.random,
                           showall = TRUE,
                           overall = TRUE,
                           ci = FALSE,
                           test = TRUE,
                           digits = gs("digits"),
                           digits.zval = gs("digits.zval"),
                           digits.pval = gs("digits.pval"),
                           text.NA = ".",
                           backtransf = x$backtransf,
                           ...) {
  
  
  meta:::chkclass(x, "netsplit")
  
  
  ## All individual results in a single row - be on the save side:
  ##
  oldopts <- options(width = 200)
  on.exit(options(oldopts))
  
  
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  meta:::chklogical(showall)
  meta:::chklogical(overall)
  meta:::chklogical(ci)
  meta:::chklogical(test)
  ##
  if (is.null(backtransf))
    backtransf <- TRUE
  meta:::chklogical(backtransf)
  
  
  sm <- x$sm
  sm.lab <- sm
  ##
  relative <- meta:::is.relative.effect(sm)
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
  
  
  if (!showall) {
    sel <- !is.na(x$direct.fixed$TE) & !is.na(x$indirect.fixed$TE) &
      !is.na(x$direct.random$TE) & !is.na(x$indirect.random$TE)
  }
  else
    sel <- rep_len(TRUE, length(x$direct.fixed$TE))
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
    TE.compare.fixed <- exp(TE.compare.fixed)
    lower.compare.fixed <- exp(lower.compare.fixed)
    upper.compare.fixed <- exp(upper.compare.fixed)
    ##
    TE.compare.random <- exp(TE.compare.random)
    lower.compare.random <- exp(lower.compare.random)
    upper.compare.random <- exp(upper.compare.random) 
  }
  
  
  fixed <- list(comp = comp,
                k = k,
                prop = format(round(prop.fixed, 2)))
  names.fixed <- c("comparison", "k", "prop")
  ##
  if (overall) {
    fixed$TE.fixed <- meta:::format.NA(TE.fixed, digits, text.NA = text.NA)
    names.fixed <- c(names.fixed, "nma")
    if (ci) {
      fixed$ci.fixed <- meta:::p.ci(round(lower.fixed, digits),
                                    round(upper.fixed, digits))
      fixed$ci.fixed[is.na(fixed$ci.fixed)] <- text.NA
      names.fixed <- c(names.fixed, ci.lab)
    }
  }
  ##
  fixed$TE.direct.fixed <- meta:::format.NA(TE.direct.fixed, digits, text.NA = text.NA)
  names.fixed <- c(names.fixed, "direct")
  if (ci) {
    fixed$ci.direct.fixed <- meta:::p.ci(round(lower.direct.fixed, digits),
                                         round(upper.direct.fixed, digits))
    fixed$ci.direct.fixed[is.na(fixed$ci.direct.fixed)] <- text.NA
    names.fixed <- c(names.fixed, ci.lab)
  }
  ##
  fixed$TE.indirect.fixed <- meta:::format.NA(TE.indirect.fixed, digits, text.NA = text.NA)
  names.fixed <- c(names.fixed, "indir.")
  if (ci) {
    fixed$ci.indirect.fixed <- meta:::p.ci(round(lower.indirect.fixed, digits),
                                           round(upper.indirect.fixed, digits))
    fixed$ci.indirect.fixed[is.na(fixed$ci.indirect.fixed)] <- text.NA
    names.fixed <- c(names.fixed, ci.lab)
  }
  ##
  if (test) {
    fixed$diff <- meta:::format.NA(TE.compare.fixed, digits, text.NA = text.NA)
    names.fixed <- c(names.fixed, if (backtransf & relative) "RoR" else "Diff")
    if (ci) {
      fixed$ci.diff <- meta:::p.ci(round(lower.compare.fixed, digits),
                                   round(upper.compare.fixed, digits))
      fixed$ci.diff[is.na(fixed$ci.diff)] <- text.NA
      names.fixed <- c(names.fixed, ci.lab)
    }
    ##
    fixed$z <- meta:::format.NA(zval.compare.fixed, digits.zval)
    fixed$z[fixed$z == "--"] <- text.NA
    fixed$p <- meta:::format.p(pval.compare.fixed, digits = digits.pval)
    fixed$p[meta:::rmSpace(fixed$p) == "--"] <- text.NA
    names.fixed <- c(names.fixed, c("z", "p-value"))
  }
  fixed <- as.data.frame(fixed)
  names(fixed) <- names.fixed
  
  
  random <- list(comp = comp,
                 k = k,
                 prop = format(round(prop.random, 2)))
  names.random <- c("comparison", "k", "prop")
  ##
  if (overall) {
    random$TE.random <- meta:::format.NA(TE.random, digits, text.NA = text.NA)
    names.random <- c(names.random, "nma")
    if (ci) {
      random$ci.random <- meta:::p.ci(round(lower.random, digits),
                                      round(upper.random, digits))
      random$ci.random[is.na(random$ci.random)] <- text.NA
      names.random <- c(names.random, ci.lab)
    }
  }
  ##
  random$TE.direct.random <- meta:::format.NA(TE.direct.random, digits, text.NA = text.NA)
  names.random <- c(names.random, "direct")
  if (ci) {
    random$ci.direct.random <- meta:::p.ci(round(lower.direct.random, digits),
                                           round(upper.direct.random, digits))
    random$ci.direct.random[is.na(random$ci.direct.random)] <- text.NA
    names.random <- c(names.random, ci.lab)
  }
  ##
  random$TE.indirect.random <- meta:::format.NA(TE.indirect.random, digits, text.NA = text.NA)
  names.random <- c(names.random, "indir.")
  if (ci) {
    random$ci.indirect.random <- meta:::p.ci(round(lower.indirect.random, digits),
                                             round(upper.indirect.random, digits))
    random$ci.indirect.random[is.na(random$ci.indirect.random)] <- text.NA
    names.random <- c(names.random, ci.lab)
  }
  ##
  if (test) {
    random$diff <- meta:::format.NA(TE.compare.random, digits, text.NA = text.NA)
    names.random <- c(names.random, if (backtransf & relative) "RoR" else "Diff")
    if (ci) {
      random$ci.diff <- meta:::p.ci(round(lower.compare.random, digits),
                                    round(upper.compare.random, digits))
      random$ci.diff[is.na(random$ci.diff)] <- text.NA
      names.random <- c(names.random, ci.lab)
    }
    ##
    random$z <- meta:::format.NA(zval.compare.random, digits.zval)
    random$z[random$z == "--"] <- text.NA
    random$p <- meta:::format.p(pval.compare.random, digits = digits.pval)
    random$p[meta:::rmSpace(random$p) == "--"] <- text.NA
    names.random <- c(names.random, c("z", "p-value"))
  }
  random <- as.data.frame(random)
  names(random) <- names.random
  
  
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
  cat("\nLegend:\n")
  cat(" comparison - Treatment comparison\n")
  cat(" k          - Number of studies providing direct evidence\n")
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
  ##
  invisible(NULL)
}
