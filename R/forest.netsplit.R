#' Forest plot for direct and indirect evidence
#' 
#' @description
#' Forest plot to show direct and indirect evidence in network
#' meta-analysis.  Furthermore, estimates from network meta-analysis
#' as well as prediction intervals can be printed.
#' 
#' @param x An object of class \code{netsplit}.
#' @param pooled A character string indicating whether results for the
#'   fixed effect (\code{"fixed"}) or random effects model
#'   (\code{"random"}) should be plotted. Can be abbreviated.
#' @param show A character string indicating which comparisons should
#'   be printed (see Details).
#' @param overall A logical indicating whether network meta-analysis
#'   estimates should be printed.
#' @param direct A logical indicating whether direct estimates should
#'   be printed.
#' @param indirect A logical indicating whether indirect estimates
#'   should be printed.
#' @param prediction A logical indicating whether prediction intervals
#'   should be printed.
#' @param subgroup A character string indicating which layout should
#'   be used in forest plot: subgroups by comparisons
#'   (\code{"comparison"}) or subgroups by estimates
#'   (\code{"estimate"}). Can be abbreviated.
#' @param text.overall A character string used in the plot to label
#'   the network estimates.
#' @param text.direct A character string used in the plot to label the
#'   direct estimates.
#' @param text.indirect A character string used in the plot to label
#'   the indirect estimates.
#' @param text.predict A character string used in the plot to label
#'   the prediction interval.
#' @param type.overall A character string specifying how to plot
#'   treatment effects and confidence intervals for the overall
#'   network evidence.
#' @param type.direct A character string specifying how to plot
#'   treatment effects and confidence intervals for the direct
#'   evidence.
#' @param type.indirect A character string specifying how to plot
#'   treatment effects and confidence intervals for the indirect
#'   evidence.
#' @param col.square The colour for squares.
#' @param col.square.lines The colour for the outer lines of squares.
#' @param col.inside The colour for results and confidence limits if
#'   confidence limits are completely within squares squares.
#' @param col.diamond The colour of diamonds.
#' @param col.diamond.lines The colour of the outer lines of diamonds.
#' @param col.predict Background colour of prediction intervals.
#' @param col.predict.lines Colour of outer lines of prediction
#'   intervals.
#' @param equal.size A logical indicating whether all squares should
#'   be of equal size. Otherwise, the square size is proportional to
#'   the precision of estimates.
#' @param leftcols A character vector specifying columns to be plotted
#'   on the left side of the forest plot (see Details).
#' @param leftlabs A character vector specifying labels for columns on
#'   left side of the forest plot.
#' @param rightcols A character vector specifying columns to be
#'   plotted on the right side of the forest plot (see Details).
#' @param rightlabs A character vector specifying labels for columns
#'   on right side of the forest plot.
#' @param digits Minimal number of significant digits for treatment
#'   effects and confidence intervals, see \code{print.default}.
#' @param digits.prop Minimal number of significant digits for the
#'   direct evidence proportion.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in forest plots. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios, for example.
#' @param lab.NA A character string to label missing values.
#' @param smlab A label printed at top of figure. By default, text
#'   indicating either fixed effect or random effects model is
#'   printed.
#' @param \dots Additional arguments for \code{\link{forest.meta}}
#'   function.
#' 
#' @details
#' A forest plot, also called confidence interval plot, is drawn in
#' the active graphics window.
#' 
#' The arguments \code{leftcols} and \code{rightcols} can be used to
#' specify columns which are plotted on the left and right side of the
#' forest plot, respectively. If argument \code{rightcols} is
#' \code{FALSE}, no columns will be plotted on the right side.
#' 
#' If direct estimates are included in the forest plot (\code{direct =
#' TRUE}, default), the following columns will be printed on the left
#' side of the forest plot: the comparisons (column \code{"studlab"}
#' in \code{\link{forest.meta}}), number of pairwise comparisons
#' (\code{"k"}), and direct evidence proportion (\code{"k"}).
#' 
#' If direct estimates are not included in the forest plot
#' (\code{direct = FALSE}), only the comparisons (\code{"studlab"})
#' are printed on the left side of the forest plot.
#' 
#' For more information see help page of \code{\link{forest.meta}}
#' function.
#' 
#' Argument \code{show} determines which comparisons are printed:
#' \tabular{ll}{
#' \dQuote{all} \tab All comparisons \cr
#' \dQuote{both} \tab Only comparisons contributing both direct and
#'   indirect evidence \cr
#' \dQuote{with.direct} \tab Comparisons providing direct evidence \cr
#' \dQuote{direct.only} \tab Comparisons providing only direct
#'   evidence \cr
#' \dQuote{indirect.only} \tab Comparisons providing only indirect
#'   evidence
#' }
#'
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{forest.meta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Senn2013)
#' #
#' net1 <- netmeta(TE, seTE, treat1.long, treat2.long,
#'                 studlab, data = Senn2013,
#'                 comb.fixed = FALSE)
#' #
#' ns1 <- netsplit(net1)
#' 
#' # Forest plot showing comparisons contributing both direct and
#' # indirect evidence
#' #
#' forest(ns1, fontsize = 6, spacing = 0.5, addrow.subgroups = FALSE)
#' 
#' \dontrun{
#' # Forest plot showing comparisons contributing direct evidence
#' #
#' forest(ns1, fontsize = 6, spacing = 0.5, addrow.subgroups = FALSE,
#'        show = "with.direct")
#' }
#' 
#' @method forest netsplit
#' @export
#' @export forest.netsplit


forest.netsplit <- function(x,
                            pooled = ifelse(x$comb.random, "random", "fixed"),
                            show = "both",
                            ##
                            subgroup = "comparison",
                            ##
                            overall = TRUE,
                            direct = TRUE,
                            indirect = TRUE,
                            prediction = x$prediction,
                            ##
                            text.overall = "Network estimate",
                            text.direct = "Direct estimate",
                            text.indirect = "Indirect estimate",
                            text.predict = "Prediction interval",
                            ##
                            type.overall,
                            type.direct,
                            type.indirect,
                            ##
                            col.square = "gray",
                            col.square.lines = col.square,
                            col.inside = "white",
                            col.diamond = "gray",
                            col.diamond.lines = "black",
                            col.predict = "red",
                            col.predict.lines = "black",
                            ##
                            equal.size = FALSE,
                            ##
                            leftcols,
                            leftlabs,
                            rightcols = c("effect", "ci"),
                            rightlabs = NULL,
                            ##
                            digits = gs("digits.forest"),
                            digits.prop = max(gs("digits.pval") - 2, 2),
                            ##
                            backtransf = x$backtransf,
                            lab.NA = "",
                            smlab,
                            ...) {


  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  meta:::chkclass(x, "netsplit")
  ##
  x <- upgradenetmeta(x)
  ##
  chkchar <- meta:::chkchar
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  formatPT <- meta:::formatPT
  setchar <- meta:::setchar
  ##
  pooled <- setchar(pooled, c("fixed", "random"))
  ##
  subgroup <- setchar(subgroup, c("comparison", "estimate"))
  ##
  chklogical(overall)
  chklogical(direct)
  chklogical(indirect)
  chklogical(prediction)
  ##
  chkchar(text.overall)
  chkchar(text.direct)
  chkchar(text.indirect)
  chkchar(text.predict)
  ##
  missing.type.overall <- missing(type.overall)
  if (missing.type.overall)
    type.overall <- "diamond"
  else
    type.overall <- setchar(type.overall, c("diamond", "square"))
  ##
  if (missing(type.direct))
    type.direct <- "square"
  else
    type.direct <- setchar(type.direct, c("diamond", "square"))
  if (missing(type.indirect))
    type.indirect <- "square"
  else
    type.indirect <- setchar(type.indirect, c("diamond", "square"))
  ##
  chkchar(col.square)
  chkchar(col.square.lines)
  chkchar(col.inside)
  chkchar(col.diamond)
  chkchar(col.diamond.lines)
  chkchar(col.predict)
  chkchar(col.predict.lines)
  ##
  chklogical(equal.size)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)
  chklogical(backtransf)
  ##
  chkchar(lab.NA)
  ##
  if (pooled == "fixed") {
    if (!(missing(prediction)) & prediction)
      warning("Prediction intervals not shown for estimates from fixed effect model.")
    prediction <- FALSE
  }
  ##
  if (!any(c(overall, direct, indirect)))
    stop("At least, one of the following estimates must be included in forest plot:\n- network estimates (argument 'overall')\n- direct estimates (argument 'direct')\n- indirect estimates (argument 'indirect')")
  ##
  if (missing(leftcols))
    if (direct)
      leftcols <- c("studlab", "k", "prop")
    else
      leftcols <- "studlab"
  ##
  if (missing(leftlabs)) {
    leftlabs <- rep(NA, length(leftcols))
    leftlabs[leftcols == "studlab"] <- "Comparison"
    leftlabs[leftcols == "k"] <- "Number of\nStudies"
    leftlabs[leftcols == "prop"] <- "Direct\nEvidence"
  }
  ##
  n.subgroup <- direct + indirect + overall + prediction
  missing.smlab <- missing(smlab)
  ##
  if (n.subgroup == 1 & overall & missing.type.overall)
    type.overall <- "square"
  ##
  if (missing(text.predict))
    if (!(length(x$level.predict) == 0) &&
        x$level.comb != x$level.predict)
      text.predict <- paste(text.predict, " (",
                            round(x$level.predict * 100), "%-PI)", sep = "")
  ##
  if (overall & n.subgroup > 1) {
    if (text.overall == text.predict)
      stop("Text must be different for arguments 'text.overall' and 'text.predict'.")
    if (text.overall == text.direct)
      stop("Text must be different for arguments 'text.overall' and 'text.direct'.")
    if (text.overall == text.indirect)
      stop("Text must be different for arguments 'text.overall' and 'text.indirect'.")
  }
  ##
  if (prediction & n.subgroup > 1) {
    if (text.predict == text.overall)
      stop("Text must be different for arguments 'text.predict' and 'text.overall'.")
    if (text.predict == text.direct)
      stop("Text must be different for arguments 'text.predict' and 'text.direct'.")
    if (text.predict == text.indirect)
      stop("Text must be different for arguments 'text.predict' and 'text.indirect'.")
  }
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


  ##
  ##
  ## (2) Extract results for fixed effect and random effects model
  ##
  ##
  if (pooled == "fixed") {
    dat.direct <- x$direct.fixed
    dat.indirect <- x$indirect.fixed
    dat.overall <- x$fixed
    ##
    dat.direct$prop <- formatPT(x$prop.fixed, digits = digits.prop)
    dat.indirect$prop <- NA
    dat.overall$prop <- NA
    ##
    if (missing.smlab)
      smlab <- "Fixed effect model"
  }
  else {
    dat.direct <- x$direct.random
    dat.indirect <- x$indirect.random
    dat.overall <- x$random
    ##
    dat.direct$prop <- formatPT(x$prop.random, digits = digits.prop)
    dat.indirect$prop <- NA
    dat.overall$prop <- NA
    ##
    if (missing.smlab)
      smlab <- "Random effects model"
  }
  ##
  if (missing.smlab & n.subgroup == 1)
    smlab <- paste(if (direct)
                     paste(text.direct, "\n", sep = ""),
                   if (indirect)
                     paste(text.indirect, "\n", sep = ""),
                   if (overall)
                     paste(text.overall, "\n", sep = ""),
                   "(",
                   tolower(smlab),
                   ")",
                   sep = "")
  ##
  dat.predict <- x$predict
  ##
  if ( (subgroup == "estimate") | (prediction & !overall) )
    dat.predict$TE <- dat.overall$TE
  else
    dat.predict$TE <- NA
  ##
  dat.predict$seTE <- dat.predict$statistic <- dat.predict$p <-
    dat.predict$z <- dat.predict$prop <- NA
  ##
  dat.predict <- dat.predict[, c("comparison", "TE", "seTE",
                                 "lower", "upper", "statistic", "p", "prop")]
  ##
  dat.direct$comps <- dat.indirect$comps <-
    dat.overall$comps <- dat.predict$comps <- x$comparison
  ##
  dat.direct$k <- x$k
  dat.indirect$k <- dat.overall$k <- dat.predict$k <- NA
  ##
  dat.direct$evidence   <- text.direct
  dat.indirect$evidence <- text.indirect
  dat.overall$evidence  <- text.overall
  dat.predict$evidence  <- text.predict
  ##
  dat.direct$type.study <- type.direct
  dat.indirect$type.study <- type.indirect
  dat.overall$type.study <- type.overall
  dat.predict$type.study <- "predict"
  ##
  dat.direct$col.estimate <- if (type.direct == "square")
                               col.square
                             else
                               col.diamond
  dat.indirect$col.estimate <- if (type.indirect == "square")
                                 col.square
                               else
                                 col.diamond
  dat.overall$col.estimate <- if (type.overall == "square")
                                col.square
                              else
                                col.diamond
  ##
  dat.direct$col.lines <- if (type.direct == "square")
                            col.square.lines
                          else
                            col.diamond.lines
  dat.indirect$col.lines <- if (type.indirect == "square")
                              col.square.lines
                            else
                              col.diamond
  dat.overall$col.lines <- if (type.overall == "square")
                             col.square.lines
                           else
                             col.diamond.lines
  ##
  dat.predict$col.estimate <- col.predict
  dat.predict$col.lines <- col.predict.lines
  ##
  ## col.square.lines = col.square,
  ## col.inside = "white",
  ## col.diamond.lines = "black",
  ## col.predict.lines = "black",
  
  
  ##
  ##
  ## (3) Select treatment comparisons to show in forest plot
  ##
  ##
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
  dat.direct <- dat.direct[sel, ]
  dat.indirect <- dat.indirect[sel, ]
  dat.overall <- dat.overall[sel, ]
  dat.predict <- dat.predict[sel, ]


  ##
  ##
  ## (4) Forest plot
  ##
  ##
  if (subgroup == "comparison") {
    dat <- rbind(
      if (direct) dat.direct,
      if (indirect) dat.indirect,
      if (overall) dat.overall,
      if (prediction) dat.predict
    )
    ##
    if (nrow(dat) == 0) {
      warning("No comparison(s) selected. Consider using argument ",
              "'show = \"all\"'.")
      return(invisible(NULL))
    }
    ##
    if (n.subgroup > 1)
      m <-
        suppressWarnings(metagen(dat$TE, dat$seTE,
                                 studlab = dat$evidence, data = dat,
                                 sm = x$sm,
                                 byvar = dat$comps, print.byvar = FALSE))
    else
      m <-
        suppressWarnings(metagen(dat$TE, dat$seTE,
                                 studlab = dat$comps, data = dat, sm = x$sm))
    ##
    if (overall) {
      m$w.fixed[m$studlab == text.overall] <- max(m$w.fixed, na.rm = TRUE)
      m$w.random[m$studlab == text.overall] <- max(m$w.random, na.rm = TRUE)
    }
    ##
    if (prediction) {
      m$lower[m$studlab == text.predict] <- dat.predict$lower
      m$upper[m$studlab == text.predict] <- dat.predict$upper
      ##
      m$w.fixed[m$studlab == text.predict] <- max(m$w.fixed, na.rm = TRUE)
      m$w.random[m$studlab == text.predict] <- max(m$w.random, na.rm = TRUE)
    }
    ##
    forest(m,
           digits = digits,
           comb.fixed = FALSE, comb.random = FALSE,
           hetstat = FALSE,
           leftcols = leftcols,
           leftlabs = leftlabs,
           rightcols = rightcols,
           rightlabs = rightlabs,
           lab.NA = lab.NA,
           smlab = smlab,
           backtransf = backtransf,
           type.study = dat$type.study,
           col.square = dat$col.estimate,
           col.square.lines = dat$col.lines,
           weight.study = if (equal.size) "same" else "fixed",
           ...)
  }
  else {
    dat <- rbind(
      if (direct) dat.direct,
      if (indirect) dat.indirect,
      if (overall) dat.overall,
      if (prediction) dat.predict
    )
    ##
    if (nrow(dat) == 0) {
      warning("No comparison(s) selected. Consider using argument ",
              "'show = \"all\"'.")
      return(invisible(NULL))
    }
    ##
    if (n.subgroup > 1)
      m <-
        suppressWarnings(metagen(dat$TE, dat$seTE,
                                 studlab = dat$comps, data = dat,
                                 sm = x$sm,
                                 byvar = dat$evidence, print.byvar = FALSE))
    else
      m <-
        suppressWarnings(metagen(dat$TE, dat$seTE,
                                 studlab = dat$comps, data = dat, sm = x$sm))
    ##
    if (overall) {
      m$w.fixed[m$byvar == text.overall] <- max(m$w.fixed, na.rm = TRUE)
      m$w.random[m$byvar == text.overall] <- max(m$w.random, na.rm = TRUE)
    }
    ##
    if (prediction) {
      m$lower[m$byvar == text.predict] <- dat.predict$lower
      m$upper[m$byvar == text.predict] <- dat.predict$upper
      ##
      m$w.fixed[m$byvar == text.predict] <- max(m$w.fixed, na.rm = TRUE)
      m$w.random[m$byvar == text.predict] <- max(m$w.random, na.rm = TRUE)
    }
    ##
    forest(m,
           digits = digits,
           comb.fixed = FALSE, comb.random = FALSE,
           hetstat = FALSE,
           leftcols = leftcols,
           leftlabs = leftlabs,
           rightcols = rightcols,
           rightlabs = rightlabs,
           lab.NA = lab.NA,
           backtransf = backtransf,
           smlab = smlab,
           type.study = dat$type.study,
           col.square = dat$col.estimate,
           col.square.lines = dat$col.lines,
           weight.study = if (equal.size) "same" else "fixed",
           ...)
  }


  invisible(NULL)
}
