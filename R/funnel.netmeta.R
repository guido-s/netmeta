#' \sQuote{Comparison-adjusted} funnel plot
#' 
#' @description
#' Draw a \sQuote{comparison-adjusted} funnel plot to assess funnel
#' plot asymmetry in network meta-analysis.
#' 
#' @param x An object of class \code{netmeta}.
#' @param order A mandatory character or numerical vector specifying
#'   the order of treatments or list of comparators (see Details).
#' @param pooled A character string indicating whether results for the
#'   common (\code{"common"}) or random effects model (\code{"random"})
#'   should be plotted. Can be abbreviated.
#' @param xlab A label for the x-axis.
#' @param level The confidence level utilised in the plot. For the
#'   funnel plot, confidence limits are not drawn if \code{yaxis =
#'   "size"}.
#' @param pch The plotting symbol(s) used for individual studies
#'   within direct comparisons.
#' @param col The colour(s) used for individual studies within direct
#'   comparisons.
#' @param legend A logical indicating whether a legend with
#'   information on direct comparisons should be added to the plot.
#' @param pos.legend The position of the legend describing plotting
#'   symbols and colours for direct comparisons.
#' @param pos.tests The position of results for test(s) of funnel plot
#'   asymmetry.
#' @param lump.comparator A logical indicating whether comparators
#'   should be lumped, e.g., to specify inactive treatments.
#'   information on direct comparisons should be added to the plot.
#' @param text.comparator A character string used in the plot to label
#'   the comparator if \code{lump.comparator} is \code{TRUE}.
#' @param method.bias A character vector indicating which test(s) for
#'   funnel plot asymmetry to use. Admissible values are
#'   \code{"Begg"}, \code{"Egger"}, and \code{"Thompson"}, can be
#'   abbreviated. See function \code{\link[meta]{metabias}}.
#' @param text.linreg A character string used in the plot to label the
#'   Egger test for funnel plot asymmetry.
#' @param text.rank A character string used in the plot to label the
#'   Begg test for funnel plot asymmetry.
#' @param text.mm A character string used in the plot to label the
#'   Thompson-Sharp test for funnel plot asymmetry.
#' @param sep.trts A character used in comparison names as separator
#'   between treatment labels.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see
#'   \code{\link{netmeta}}).
#' @param backtransf A logical indicating whether results for relative
#'   summary measures (argument \code{sm} equal to \code{"OR"},
#'   \code{"RR"}, \code{"HR"}, or \code{"IRR"}) should be back
#'   transformed in funnel plots. If \code{backtransf = TRUE}, results
#'   for \code{sm = "OR"} are printed as odds ratios rather than log
#'   odds ratios, for example.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of test(s) for funnel plot asymmetry.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param linreg Deprecated argument (replaced by \code{method.bias}).
#' @param rank Deprecated argument (replaced by \code{method.bias}).
#' @param mm Deprecated argument (replaced by \code{method.bias}).
#' @param \dots Additional graphical arguments passed as arguments to
#'   \code{\link[meta]{funnel.meta}}.
#' 
#' @details
#' A \sQuote{comparison-adjusted} funnel plot (Chaimani & Salanti,
#' 2012) is drawn in the active graphics window.
#' 
#' Argument \code{order} is mandatory to determine the order of
#' treatments (Chaimani et al., 2013):
#' 
#' \emph{\dQuote{Before using this plot, investigators should order
#' the treatments in a meaningful way and make assumptions about how
#' small studies differ from large ones. For example, if they
#' anticipate that newer treatments are favored in small trials, then
#' they could name the treatments from oldest to newest so that all
#' comparisons refer to \sQuote{old versus new intervention}. Other
#' possibilities include defining the comparisons so that all refer to
#' an active treatment versus placebo or sponsored versus
#' non-sponsored intervention.}}
#' 
#' Alternatively, it is possible to only provide a single or few
#' treatment name(s) in argument \code{order} to define the
#' comparator(s). In this case only comparisons with this / these
#' treatment(s) will be considered. If argument \code{lump.comparator}
#' is \code{TRUE}, all comparators will be lumped into a single
#' group. The text for this group can be specified with argument
#' \code{text.comparator}.
#' 
#' In the funnel plot, if \code{yaxis} is \code{"se"}, the standard
#' error of the treatment estimates is plotted on the y-axis which is
#' likely to be the best choice (Sterne & Egger, 2001). Other possible
#' choices for \code{yaxis} are \code{"invvar"} (inverse of the
#' variance), \code{"invse"} (inverse of the standard error), and
#' \code{"size"} (study size).
#'
#' @return
#' A data frame with the following columns:
#' \item{studlab}{Study label.}
#' \item{treat1}{Label/Number for first treatment.}
#' \item{treat2}{Label/Number for second treatment.}
#' \item{comparison}{Treatment comparison.}
#' \item{TE}{Estimate of treatment effect, e.g., log odds ratio.}
#' \item{TE.direct}{Pooled estimate from direct evidence.}
#' \item{TE.adj}{\sQuote{Comparison-adjusted} treatment effect (TE -
#'   TE.direct).}
#' \item{seTE}{Standard error of treatment estimate.}
#' \item{pch}{Plotting symbol(s).}
#' \item{col}{Colour of plotting symbol(s).}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link[meta]{funnel.meta}},
#'   \code{\link[meta]{metabias}}
#' 
#' @references
#' Chaimani A & Salanti G (2012):
#' Using network meta-analysis to evaluate the existence of
#' small-study effects in a network of interventions.
#' \emph{Research Synthesis Methods},
#' \bold{3}, 161--76
#' 
#' Chaimani A, Higgins JP, Mavridis D, Spyridonos P, Salanti G (2013):
#' Graphical tools for network meta-analysis in STATA.
#' PLOS ONE,
#' \bold{8}, e76654
#' 
#' Sterne JAC & Egger M (2001):
#' Funnel plots for detecting bias in meta-analysis: Guidelines on
#' choice of axis.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{54}, 1046--55
#' 
#' @keywords hplot
#' 
#' @examples
#' \dontrun{
#' data(Senn2013)
#' 
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD")
#' 
#' # 'Comparison-adjusted' funnel plot not created as argument 'order'
#' # is missing
#' #
#' try(funnel(net1))
#' 
#' # Only show comparisons with placebo
#' #
#' funnel(net1, order = "pl")
#' 
#' # Add result for Egger test of funnel plot asymmetry
#' #
#' funnel(net1, order = "pl", method.bias = "Egger",
#'   digits.pval = 2)
#' 
#' # (Non-sensical) alphabetic order of treatments with placebo as
#' # last treatment
#' #
#' ord <- c("a", "b", "me", "mi", "pi", "r", "si", "su", "v", "pl")
#' funnel(net1, order = ord)
#'
#' # Add results for tests of funnel plot asymmetry and use different
#' # plotting symbols and colours
#' #
#' funnel(net1, order = ord,
#'   pch = rep(c(15:18, 1), 3), col = 1:3,
#'   method.bias = c("Egger", "Begg", "Thompson"), digits.pval = 2)
#' 
#' # Same results for tests of funnel plot asymmetry using reversed
#' # order of treatments
#' #
#' funnel(net1, order = rev(ord),
#'   pch = rep(c(15:18, 1), 3), col = 1:3,
#'   method.bias = c("Egger", "Begg", "Thompson"), digits.pval = 2)
#' 
#' # Calculate tests for funnel plot asymmetry
#' #
#' f1 <- funnel(net1, order = ord)
#' #
#' metabias(metagen(TE.adj, seTE, data = f1))
#' metabias(metagen(TE.adj, seTE, data = f1), method = "Begg")
#' metabias(metagen(TE.adj, seTE, data = f1), method = "Thompson")
#' }
#'
#' @method funnel netmeta
#' @export


funnel.netmeta <- function(x,
                           order,
                           pooled = ifelse(x$random, "random", "common"),
                           ##
                           xlab = NULL,
                           level = x$level,
                           ##
                           pch,
                           col = "black",
                           ##
                           legend = TRUE,
                           ##
                           pos.legend = "topright",
                           pos.tests = "topleft",
                           ##
                           lump.comparator = FALSE,
                           text.comparator = "comparator",
                           ##
                           method.bias,
                           text.linreg = "(Egger)",
                           text.rank = "(Begg-Mazumdar)",
                           text.mm = "(Thompson-Sharp)",
                           ##
                           sep.trts = x$sep.trts,
                           nchar.trts = x$nchar.trts,
                           ##
                           backtransf = x$backtransf,
                           digits.pval = gs("digits.pval"),
                           ##
                           warn.deprecated = gs("warn.deprecated"),
                           linreg = FALSE,
                           rank = FALSE,
                           mm = FALSE,
                           ##
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
  chkclass(x, "netmeta")
  x <- updateversion(x)
  ##
  pooled <- setchar(pooled, c("common", "random", "fixed"))
  pooled[pooled == "fixed"] <- "common"
  chklogical(lump.comparator)
  ##
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  if (missing(order))
    stop("Argument 'order' with a meaningful order of treatments ",
         "must be provided.\n  ",
         "(see help page of funnel.netmeta for some examples).")
  else {
    order <- catch("order", mc, x, sfsp)
    ##
    if (length(order) != length(x$trts)) {
      order <- setchar(order, x$trts)
      order.all <- c(x$trts[!(x$trts %in% order)], order)
    }
    else
      order.all <- order <- setseq(order, x$trts)
  }
  ##
  if ((length(order) == length(order.all)) & lump.comparator) {
    warning("Argument 'lump.comparator' ignored as full treatment order is ",
            "provided.",
            call. = FALSE)
    lump.comparator <- FALSE
  }
  ##
  chklogical(warn.deprecated)
  ##
  missing.method.bias <- missing(method.bias)
  if (missing.method.bias)
    method.bias <- character(0)
  ##
  chklogical(legend)
  ##
  deprecated2(method.bias, missing.method.bias,
              linreg, missing(linreg), warn.deprecated)
  deprecated2(method.bias, missing.method.bias,
              rank, missing(rank), warn.deprecated)
  deprecated2(method.bias, missing.method.bias,
              mm, missing(mm), warn.deprecated)
  ##
  chklogical(linreg)
  chklogical(rank)
  chklogical(mm)
  ##
  if (missing.method.bias) {
    if (linreg)
      method.bias <- "Egger"
    if (rank)
      method.bias <- c(method.bias, "Begg")
    if (mm)
      method.bias <- c(method.bias, "Thompson")
  }
  else
    method.bias <- setmethodbias(method.bias, 1:3)
  ##
  linreg <- "Egger" %in% method.bias
  rank <- "Begg" %in% method.bias
  mm <- "Thompson" %in% method.bias
  
  
  chkchar(text.comparator)
  chkchar(text.linreg)
  chkchar(text.rank)
  chkchar(text.mm)
  ##
  chkchar(sep.trts)
  chknumeric(nchar.trts, min = 1, length = 1)
  ##
  chklogical(backtransf)
  ##
  chknumeric(digits.pval, min = 1, length = 1)
  
  
  ##
  ##
  ## (2) Get data
  ##
  ##
  TE <- x$TE
  seTE <- x$seTE
  studlab <- x$studlab
  ##
  treat1 <- x$treat1
  treat2 <- x$treat2
  ##
  comparison <- paste(treat1, treat2, sep = sep.trts)
  comparison21 <- paste(treat2, treat1, sep = sep.trts)
  ##
  treat1.pos <- as.numeric(factor(treat1, levels = order.all))
  treat2.pos <- as.numeric(factor(treat2, levels = order.all))
  ##
  trts.abbr <- treats(x$trts, nchar.trts)
  trt1 <- as.character(factor(treat1, levels = x$trts, labels = trts.abbr))
  trt2 <- as.character(factor(treat2, levels = x$trts, labels = trts.abbr))
  ##
  comp <- paste(trt1, trt2, sep = sep.trts)
  comp21 <- paste(trt2, trt1, sep = sep.trts)
  ##
  wo <- treat1.pos > treat2.pos
  ##
  if (any(wo)) {
    ttreat1.pos <- treat1.pos
    treat1.pos[wo] <- treat2.pos[wo]
    treat2.pos[wo] <- ttreat1.pos[wo]
    ##
    TE[wo] <- -TE[wo]
    ##
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
    ##
    comparison[wo] <- comparison21[wo]
    ##
    ttrt1 <- trt1
    trt1[wo] <- trt2[wo]
    trt2[wo] <- ttrt1[wo]
    ##
    comp[wo] <- comp21[wo]
  }
  ##
  o <- order(treat1.pos, treat2.pos)
  ##
  TE <- TE[o]
  seTE <- seTE[o]
  ##
  studlab <- studlab[o]
  ##
  treat1 <- treat1[o]
  treat2 <- treat2[o]
  comparison <- comparison[o]
  ##
  trt1 <- trt1[o]
  trt2 <- trt2[o]
  comp <- comp[o]
  ##
  res <- data.frame(studlab,
                    treat1, treat2, comparison,
                    trt1, trt2, comp,
                    TE, TE.direct = NA, TE.adj = NA, seTE)
  ##
  if (missing(xlab)) {
    if (xlab(x$sm, backtransf) == "")
      xlab <- "Centered at comparison-specific effect"
    else
      xlab <- paste(xlab(x$sm, backtransf),
                    "centered at\ncomparison-specific effect")
  }
  
  
  ##
  ##
  ## (3) Calculate 'comparison-adjusted' treatment effects
  ##
  ##
  if (is.numeric(treat1))
    treat1 <- as.character(treat1)
  if (is.numeric(treat2))
    treat2 <- as.character(treat2)
  ##
  if (pooled == "common")
    for (i in seq_along(res$TE))
      res$TE.direct[i] <- x$TE.direct.common[treat1[i], treat2[i]]
  else
    for (i in seq_along(res$TE))
      res$TE.direct[i] <- x$TE.direct.random[treat1[i], treat2[i]]
  ##
  res$TE.adj <- res$TE - res$TE.direct
  ##
  if (length(order) != length(order.all))
    res <- subset(res, treat1 %in% order | treat2 %in% order)
  ##
  if (lump.comparator) {
    res <- subset(res, !(treat1 %in% order & treat2 %in% order))
    ##
    res$trt2[res$treat2 %in% order] <- text.comparator
    res$treat2[res$treat2 %in% order] <- text.comparator
    ##
    res$comp <- paste(res$trt1, res$trt2, sep = sep.trts)
    res$comparison <- paste(res$treat1, res$treat2, sep = sep.trts)
  }
  
  
  ##
  ##
  ## (4) Calculate necessary data for funnel plot
  ##
  ##
  m.adj <-
    suppressWarnings(metagen(res$TE.adj, res$seTE,
                             studlab = res$studlab, sm = x$sm,
                             method.tau = "DL", method.tau.ci = ""))
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
         common = FALSE, random = FALSE,
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
    mb.linreg <- suppressWarnings(metabias(m.adj, method = "Egger"))
  if (rank)
    mb.rank <- suppressWarnings(metabias(m.adj, method = "Begg"))
  if (mm)
    mb.mm <- suppressWarnings(metabias(m.adj, method = "Thompson"))
  ##
  if (linreg | rank | mm)
    legend(pos.tests,
           legend = paste(formatPT(c(if (linreg) mb.linreg$p.value,
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
  attr(res, "lump.comparator") <- lump.comparator
  
  
  invisible(res)
}
