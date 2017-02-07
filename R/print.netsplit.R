print.netsplit <- function(x,
                           comb.fixed = x$comb.fixed,
                           comb.random = x$comb.random,
                           digits = gs("digits"),
                           digits.zval = gs("digits.zval"),
                           digits.pval = gs("digits.pval"),
                           text.NA = ".", backtransf = TRUE,
                           showall = FALSE,
                           ...) {
  
  meta:::chkclass(x, "netsplit")
  
  
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
    sel <- !is.na(x$TE.direct.fixed) & !is.na(x$TE.indirect.fixed) &
      !is.na(x$TE.direct.random) & !is.na(x$TE.indirect.random)
  }
  else
    sel <- rep_len(TRUE, length(x$TE.direct.fixed))
  ## 
  comparison <- x$comparison[sel]
  ##
  prop.fixed <- x$prop.fixed[sel]
  TE.direct.fixed <- x$TE.direct.fixed[sel]
  TE.indirect.fixed <- x$TE.indirect.fixed[sel]
  TE.diff.fixed <- x$TE.diff.fixed[sel]
  lower.diff.fixed <- x$lower.diff.fixed[sel]
  upper.diff.fixed <- x$upper.diff.fixed[sel]
  zval.diff.fixed <- x$zval.diff.fixed[sel]
  pval.diff.fixed <- x$pval.diff.fixed[sel]
  ##
  prop.random <- x$prop.random[sel]
  TE.direct.random <- x$TE.direct.random[sel]
  TE.indirect.random <- x$TE.indirect.random[sel]
  TE.diff.random <- x$TE.diff.random[sel]
  lower.diff.random <- x$lower.diff.random[sel]
  upper.diff.random <- x$upper.diff.random[sel]
  zval.diff.random <- x$zval.diff.random[sel]
  pval.diff.random <- x$pval.diff.random[sel]
  
  
  if (backtransf & relative) {
    TE.direct.fixed <- exp(TE.direct.fixed)
    TE.indirect.fixed <- exp(TE.indirect.fixed)
    TE.diff.fixed <- exp(TE.diff.fixed)
    lower.diff.fixed <- exp(lower.diff.fixed)
    upper.diff.fixed <- exp(upper.diff.fixed)
    TE.direct.random <- exp(TE.direct.random)
    TE.indirect.random <- exp(TE.indirect.random)
    TE.diff.random <- exp(TE.diff.random)
    lower.diff.random <- exp(lower.diff.random)
    upper.diff.random <- exp(upper.diff.random) 
  }
  
  
  fixed <- data.frame(comparison,
                      prop = format(round(prop.fixed, 2)),
                      TE.direct.fixed = meta:::format.NA(TE.direct.fixed, digits, text.NA = text.NA),
                      TE.indirect.fixed = meta:::format.NA(TE.indirect.fixed, digits, text.NA = text.NA),
                      diff = meta:::format.NA(TE.diff.fixed, digits, text.NA = text.NA),
                      ci = meta:::p.ci(round(lower.diff.fixed, digits),
                                       round(upper.diff.fixed, digits)),
                      z = meta:::format.NA(zval.diff.fixed, digits.zval),
                      p = meta:::format.p(pval.diff.fixed, digits = digits.pval),
                      stringsAsFactors = FALSE
                      )
  ##
  fixed$z[fixed$z == "--"] <- text.NA
  fixed$p[meta:::rmSpace(fixed$p) == "--"] <- text.NA
  ##
  names(fixed) <- c("comparison", "prop", "direct", "indir.",
                    if (backtransf & relative) "RoR" else "Diff", ci.lab, "z", "p-value")
  ##
  random <- data.frame(comparison,
                       prop = format(round(prop.random, 2)),
                       TE.dir = meta:::format.NA(TE.direct.random, digits, text.NA = text.NA),
                       TE.indirect.fixed = meta:::format.NA(TE.indirect.fixed, digits, text.NA = text.NA),
                       diff = meta:::format.NA(TE.diff.random, digits, text.NA = text.NA),
                       ci = meta:::p.ci(round(lower.diff.random, digits),
                                        round(upper.diff.random, digits)),
                       z = meta:::format.NA(zval.diff.random, digits.zval),
                       p = meta:::format.p(pval.diff.random, digits = digits.pval),
                       stringsAsFactors = FALSE
                       )
  ##
  random$z[random$z == "--"] <- text.NA
  random$p[meta:::rmSpace(random$p) == "--"] <- text.NA
  ##
  names(random) <- c("comparison", "prop", "direct", "indir.",
                     if (backtransf & relative) "RoR" else "Diff", ci.lab, "z", "p-value")
  ##
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
  cat(" prop       - Direct evidence proportion\n")
  cat(paste(" direct     - Estimated treatment effect ", sm.lab,
            "derived from direct evidence\n", sep = ""))
  cat(paste(" indir.     - Estimated treatment effect ", sm.lab,
            "derived from indirect evidence\n", sep = ""))
  if (backtransf & relative)
    cat(" RoR        - Ratio of Ratios (direct versus indirect)\n")
  else
    cat(" Diff       - Difference between direct and indirect treatment estimates\n")
  ##
  invisible(NULL)
}
