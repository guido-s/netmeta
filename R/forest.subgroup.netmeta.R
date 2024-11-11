#' Forest plot showing results of network meta-analysis with subgroups
#' 
#' @description
#' Forest plot to show subgroup estimates of network meta-analysis.
#'
#' @aliases forest.subgroup.netmeta plot.subgroup.netmeta
#' 
#' @param x An object of class \code{subgroup.netmeta}.
#' @param pooled A character string indicating whether results for the
#'   common (\code{"common"}) or random effects model (\code{"random"})
#'   should be plotted. Can be abbreviated.
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
#' @param calcwidth.subgroup A logical indicating whether text with
#'   comparison labels should be considered to calculate width of the
#'   column with treatment labels.
#' @param digits Minimal number of significant digits for treatment
#'   effects and confidence intervals, see \code{print.default}.
#' @param digits.Q Minimal number of significant digits for
#'   heterogeneity statistic Q, see \code{print.default}.
#' @param digits.pval.Q Minimal number of significant digits for
#'   p-value of heterogeneity test, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance \eqn{\tau^2}, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for
#'   \eqn{\tau}, the square root of the between-study variance
#'   \eqn{\tau^2}.
#' @param sep.trts A character string used to label treatment comparisons.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in forest plots. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios, for example.
#' @param lab.NA A character string to label missing values.
#' @param smlab A label printed at top of figure. By default, text
#'   indicating either common or random effects model is printed.
#' @param col.subgroup The colour to print information on subgroups.
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
#' For more information see help page of \code{\link{forest.meta}}
#' function.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{subgroup.netmeta}}, \code{\link{netmeta}},
#'   \code{\link{forest.meta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' # Add examples
#' 
#' @method forest subgroup.netmeta
#' @export


forest.subgroup.netmeta <- function(x,
                                    pooled =
                                      ifelse(x$x$random, "random", "common"),
                                    #
                                    equal.size = gs("equal.size"),
                                    #
                                    leftcols =
                                      c("studlab",
                                        "Q", "df.Q", "pval.Q"),
                                    leftlabs =
                                      c("Comparison /\nSubgroup",
                                        "Q", "d.f.", "p-value"),
                                    rightcols =
                                      c("effect", "ci", "k",
                                        if (pooled == "random") "tau"),
                                    rightlabs =
                                      c(NA, NA, "Number of\nStudies",
                                        if (pooled == "random") "Tau"),
                                    #
                                    calcwidth.subgroup =
                                      gs("calcwidth.subgroup"),
                                    #
                                    digits = gs("digits.forest"),
                                    digits.Q = gs("digits.Q"),
                                    digits.pval.Q = gs("digits.pval.Q"),
                                    digits.tau2 = gs("digits.tau2"),
                                    digits.tau = gs("digits.tau"),
                                    #
                                    sep.trts = " vs ",
                                    backtransf = x$x$backtransf,
                                    lab.NA = ".",
                                    smlab,
                                    #
                                    col.subgroup = "black",
                                    ...) {
  
  
  #
  #
  # (1) Check and set arguments
  #
  #
  
  chkclass(x, "subgroup.netmeta")
  #
  pooled <- setchar(pooled, c("common", "random", "fixed"))
  pooled[pooled == "fixed"] <- "common"
  #
  chklogical(equal.size)
  chklogical(calcwidth.subgroup)
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.pval.Q, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  #
  chkchar(sep.trts, length = 1)
  chklogical(backtransf)
  chkchar(lab.NA, length = 1)
  
  
  #
  #
  # (2) Extract results for common and random effects model
  #
  #
  
  if (pooled == "common") {
    common <- x$common %>% filter(is.na(df.Q) | df.Q > 0)
    #
    common$Q <- formatN(common$Q, digits = digits.Q, text.NA = "")
    common$df.Q <- formatN(common$df.Q, digits = 0, text.NA = "")
    common$pval.Q <-
      formatPT(common$pval.Q, digits = digits.pval.Q, lab.NA = "")
    #
    m <-
      suppressWarnings(metagen(TE, seTE, studlab = subgroup, data = common,
                               subgroup = paste(treat1, treat2, sep = sep.trts),
                               print.subgroup.name = FALSE,
                               sm = x$x$sm,
                               overall = FALSE, overall.hetstat = FALSE,
                               common = FALSE, random = FALSE))
    #
    text.pooled <- "Common Effects Model"
  }
  else {
    random <- x$random %>% filter(is.na(df.Q) | df.Q > 0)
    #
    random$Q <- formatN(random$Q, digits = digits.Q, text.NA = "")
    random$df.Q <- formatN(random$df.Q, digits = 0, text.NA = "")
    random$pval.Q <-
      formatPT(random$pval.Q, digits = digits.pval.Q, lab.NA = "")
    random$tau2 <- formatPT(random$tau2, digits = digits.tau2, lab.NA = lab.NA)
    random$tau <- formatPT(random$tau, digits = digits.tau, lab.NA = lab.NA)
    #
    m <-
      suppressWarnings(metagen(TE, seTE, studlab = subgroup, data = random,
                               subgroup = paste(treat1, treat2, sep = " vs "),
                               print.subgroup.name = FALSE,
                               sm = x$x$sm,
                               overall = FALSE, overall.hetstat = FALSE,
                               common = FALSE, random = FALSE))
    #
    text.pooled <- "Random Effects Model"
  }
  #
  if (missing(smlab))
      smlab <- text.pooled

  
  #
  #
  # (3) Forest plot
  #
  #
  
  # Get rid of warning 'Undefined global functions or variables'
  .seTE <- .studlab <- .TE <- .treat1 <- .treat2 <- comparison <-
    df.Q <- subnet <- TE <- seTE <- treat1 <- treat2 <- NULL
  
  forest(m,
         digits = digits,
         #
         overall = FALSE, common = FALSE, random = FALSE,
         hetstat = FALSE, test.subgroup = FALSE,
         #
         subgroup.hetstat = FALSE,
         prediction.subgroup = FALSE,
         calcwidth.subgroup = calcwidth.subgroup,
         #
         leftcols = leftcols,
         leftlabs = leftlabs,
         rightcols = rightcols,
         rightlabs = rightlabs,
         #
         lab.NA = lab.NA,
         smlab = smlab,
         backtransf = backtransf,
         #
         col.subgroup = col.subgroup,
         #
         weight.study = if (equal.size) "same" else pooled,
         ...)
  
  ret <- m
  #
  ret$leftcols <- leftcols
  ret$rightcols <- rightcols
  ret$leftlabs <- leftlabs
  ret$rightlabs <- rightlabs
  #
  invisible(ret)
}





#' @rdname forest.subgroup.netmeta
#' @method plot subgroup.netmeta
#' @export
#'

plot.subgroup.netmeta <- function(x, ...)
  forest(x, ...)
