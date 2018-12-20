funnel.netmeta <- function(x,
                           order,
                           pooled = ifelse(x$comb.random, "random", "fixed"),
                           ##
                           xlab,
                           level = x$level,
                           ##
                           pch,
                           col = "black",
                           ##
                           legend = TRUE,
                           linreg = FALSE,
                           rank = FALSE,
                           mm = FALSE,
                           ##
                           pos.legend = "topright",
                           pos.tests = "topleft",
                           ##
                           text.linreg = "(Egger)",
                           text.rank = "(Begg-Mazumdar)",
                           text.mm = "(Thompson-Sharp)",
                           ##
                           sep.trts = x$sep.trts,
                           nchar.trts = x$nchar.trts,
                           ##
                           backtransf = x$backtransf,
                           digits.pval = gs("digits.pval"),
                           ...) {
  
  
  ##
  ##
  ## Generate 'comparison-adjusted' funnel plot according to
  ## Chaimani, Anna, Julian P T Higgins, Dimitris Mavridis, Panagiota
  ## Spyridonos, and Georgia Salanti. 2013. Graphical tools for network
  ## meta-analysis in STATA. PLoS One 8 (10): e76654
  ##
  ##
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  chkchar <- meta:::chkchar
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  ##
  meta:::chkclass(x, "netmeta")
  x <- upgradenetmeta(x)
  ##
  pooled <- meta:::setchar(pooled, c("fixed", "random"))
  ##
  if (missing(order)) {
    warning("In order to construct a 'comparison-adjusted' funnel plot,\n",
            "  please provide a meaningful order of treatments using argument 'order'\n",
            "  (see help page of funnel.netmeta for some examples).")
    return(invisible(NULL))
  }
  else
    order <- setseq(order, x$trts)
  ##
  chklogical(legend)
  chklogical(linreg)
  chklogical(rank)
  chklogical(mm)
  ##
  chkchar(text.linreg)
  chkchar(text.rank)
  chkchar(text.mm)
  ##
  chkchar(sep.trts)
  chknumeric(nchar.trts, min = 1, single = TRUE)
  ##
  chklogical(backtransf)
  ##
  chknumeric(digits.pval, min = 1, single = TRUE)
  
  
  ##
  ##
  ## (2) Get data
  ##
  ##
  TE <- x$TE
  seTE <- x$seTE
  treat1 <- x$treat1
  treat2 <- x$treat2
  studlab <- x$studlab
  ##
  trts.abbr <- treats(x$trts, nchar.trts)
  trt1 <- as.character(factor(treat1, levels = x$trts, labels = trts.abbr))
  trt2 <- as.character(factor(treat2, levels = x$trts, labels = trts.abbr))
  ##
  comp <- paste(trt1, trt2, sep = sep.trts)
  comp21 <- paste(trt2, trt1, sep = sep.trts)
  ##
  comparison <- paste(treat1, treat2, sep = sep.trts)
  comparison21 <- paste(treat2, treat1, sep = sep.trts)
  ##
  treat1.pos <- as.numeric(factor(treat1, levels = order))
  treat2.pos <- as.numeric(factor(treat2, levels = order))
  ##
  wo <- treat1.pos > treat2.pos
  ##
  if (any(wo)) {
    TE[wo] <- -TE[wo]
    ##
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
    ##
    ttreat1.pos <- treat1.pos
    treat1.pos[wo] <- treat2.pos[wo]
    treat2.pos[wo] <- ttreat1.pos[wo]
    ##
    comp[wo] <- comp21[wo]
    comparison[wo] <- comparison21[wo]
  }
  ##
  o <- order(treat1.pos, treat2.pos)
  ##
  TE <- TE[o]
  seTE <- seTE[o]
  treat1 <- treat1[o]
  treat2 <- treat2[o]
  studlab <- studlab[o]
  comp <- comp[o]
  comparison <- comparison[o]
  ##
  res <- data.frame(studlab,
                    treat1, treat2, comparison,
                    trt1, trt2, comp,
                    TE, TE.direct = NA, TE.adj = NA, seTE)
  ##
  if (missing(xlab)) {
    if (meta:::xlab(x$sm, backtransf) == "")
      xlab <- "Centered at comparison-specific effect"
    else
      xlab <- paste(meta:::xlab(x$sm, backtransf),
                    "centered at\ncomparison-specific effect")
  }
  
  
  ##
  ##
  ## (3) Calculate 'comparison-adjusted' treatment effects
  ##
  ##
  if (pooled == "fixed")
    for (i in seq_along(res$TE))
      res$TE.direct[i] <- x$TE.direct.fixed[treat1[i], treat2[i]]
  else
    for (i in seq_along(res$TE))
      res$TE.direct[i] <- x$TE.direct.random[treat1[i], treat2[i]]
  ##
  res$TE.adj <- res$TE - res$TE.direct
  
  
  ##
  ##
  ## (4) Calculate necessary data for funnel plot
  ##
  ##
  m.adj <- metagen(res$TE.adj, res$seTE, studlab = res$studlab, sm = x$sm)
  ##
  n.comps <- length(unique(res$comparison))
  ##
  ## Argument 'pch'
  ##
  if (missing(pch))
    pch <- seq_len(n.comps)
  else {
    bad <- FALSE
    ##
    if (length(pch) != n.comps) {
      if (length(pch) > n.comps)
        bad <- TRUE
      else if (n.comps %% length(pch) != 0)
        bad <- TRUE
      ##
      if (bad)
        stop("Length of argument 'pch' (", length(pch),
             ") does not match the number of direct pairwise comparisons (",
             n.comps, ").")
      ##
      pch <- cbind(seq_len(n.comps), pch)[, 2]
    }
  }
  ##
  ## Argument 'col'
  ##
  bad <- FALSE
  ##
  if (length(col) != n.comps) {
    if (length(col) > n.comps)
      bad <- TRUE
    else if (n.comps %% length(col) != 0)
      bad <- TRUE
    ##
    if (bad)
      stop("Length of argument 'col' (", length(col),
           ") does not match the number of direct pairwise comparisons (",
           n.comps, ").")
    ##
    col <- cbind(seq_len(n.comps), col)[, 2]
  }
  ##
  res$col <- res$pch <- as.numeric(factor(res$comparison,
                                          levels = unique(res$comparison)))
  ##
  res$pch <- pch[res$pch]
  res$col <- col[res$col]
  ##
  if (is.numeric(col))
    res$col <- as.numeric(res$col)
  ##
  if (is.numeric(pch))
    res$pch <- as.numeric(res$pch)
  
  
  ##
  ##
  ## (5) Funnel plot
  ##
  ##
  funnel(m.adj,
         pch = res$pch,
         col = res$col,
         level = level,
         comb.fixed = FALSE, comb.random = FALSE,
         xlab = xlab,
         backtransf = backtransf,
         ref.triangle = TRUE,
         ...
         )
  ##
  if (legend) {
    d1 <- unique(res[, c("comp", "pch", "col"), drop = FALSE])
    legend(pos.legend, legend = d1$comp,
           pch = d1$pch, col = d1$col)
  }
  ##
  if (linreg)
    mb.linreg <- metabias(m.adj, method = "linreg")
  if (rank)
    mb.rank <- metabias(m.adj, method = "rank")
  if (mm)
    mb.mm <- metabias(m.adj, method = "mm")
  ##
  if (linreg | rank | mm)
    legend("topleft",
           legend = paste(meta:::formatPT(c(if (linreg) mb.linreg$p.value,
                                            if (rank) mb.rank$p.value,
                                            if (mm) mb.mm$p.value),
                                          digits = digits.pval, lab = TRUE),
                          c(if (linreg) text.linreg,
                            if (rank) text.rank,
                            if (mm) text.mm)))
  
  
  ##
  if (all(res$comparison == res$comp)) {
    res$trt1 <- NULL
    res$trt2 <- NULL
    res$comp <- NULL
  }
  ##
  attr(res, "pooled") <- pooled
  attr(res, "order") <- order
  
  
  invisible(res)
}
