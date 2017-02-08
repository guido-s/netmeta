print.netsplit <- function(x,
                           comb.fixed = x$comb.fixed,
                           comb.random = x$comb.random,
                           digits = gs("digits"),
                           digits.zval = gs("digits.zval"),
                           digits.pval = gs("digits.pval"),
                           text.NA = ".", backtransf = TRUE,
                           showall = TRUE,
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
    sel <- !is.na(x$direct.fixed$TE) & !is.na(x$indirect.fixed$TE) &
      !is.na(x$direct.random$TE) & !is.na(x$indirect.random$TE)
  }
  else
    sel <- rep_len(TRUE, length(x$direct.fixed$TE))
  ## 
  comparison <- x$comparison[sel]
  ##
  prop.fixed <- x$prop.fixed[sel]
  TE.direct.fixed <- x$direct.fixed$TE[sel]
  TE.indirect.fixed <- x$indirect.fixed$TE[sel]
  TE.compare.fixed <- x$compare.fixed$TE[sel]
  lower.compare.fixed <- x$compare.fixed$lower[sel]
  upper.compare.fixed <- x$compare.fixed$upper[sel]
  zval.compare.fixed <- x$compare.fixed$z[sel]
  pval.compare.fixed <- x$compare.fixed$p[sel]
  ##
  prop.random <- x$prop.random[sel]
  TE.direct.random <- x$direct.random$TE[sel]
  TE.indirect.random <- x$indirect.random$TE[sel]
  TE.compare.random <- x$compare.random$TE[sel]
  lower.compare.random <- x$compare.random$lower[sel]
  upper.compare.random <- x$compare.random$upper[sel]
  zval.compare.random <- x$compare.random$z[sel]
  pval.compare.random <- x$compare.random$p[sel]
  
  
  if (backtransf & relative) {
    TE.direct.fixed <- exp(TE.direct.fixed)
    TE.indirect.fixed <- exp(TE.indirect.fixed)
    TE.compare.fixed <- exp(TE.compare.fixed)
    lower.compare.fixed <- exp(lower.compare.fixed)
    upper.compare.fixed <- exp(upper.compare.fixed)
    TE.direct.random <- exp(TE.direct.random)
    TE.indirect.random <- exp(TE.indirect.random)
    TE.compare.random <- exp(TE.compare.random)
    lower.compare.random <- exp(lower.compare.random)
    upper.compare.random <- exp(upper.compare.random) 
  }
  
  
  fixed <- data.frame(comparison,
                      prop = format(round(prop.fixed, 2)),
                      TE.direct.fixed = meta:::format.NA(TE.direct.fixed, digits, text.NA = text.NA),
                      TE.indirect.fixed = meta:::format.NA(TE.indirect.fixed, digits, text.NA = text.NA),
                      diff = meta:::format.NA(TE.compare.fixed, digits, text.NA = text.NA),
                      ci = meta:::p.ci(round(lower.compare.fixed, digits),
                                       round(upper.compare.fixed, digits)),
                      z = meta:::format.NA(zval.compare.fixed, digits.zval),
                      p = meta:::format.p(pval.compare.fixed, digits = digits.pval),
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
                       diff = meta:::format.NA(TE.compare.random, digits, text.NA = text.NA),
                       ci = meta:::p.ci(round(lower.compare.random, digits),
                                        round(upper.compare.random, digits)),
                       z = meta:::format.NA(zval.compare.random, digits.zval),
                       p = meta:::format.p(pval.compare.random, digits = digits.pval),
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
