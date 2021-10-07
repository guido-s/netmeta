#' Forest plot showing results of two or more network meta-analyses
#' 
#' @description
#' Forest plot to show network estimates of two or more network
#' meta-analyses.
#'
#' @aliases forest.netbind plot.netbind
#' 
#' @param x An object of class \code{netbind}.
#' @param pooled A character string indicating whether results for the
#'   fixed (\code{"fixed"}) or random effects model (\code{"random"})
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
#' @param subset.treatments A character vector specifying treatments
#'   to show in forest plot as comparators to the reference.
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
#'   indicating either fixed or random effects model is printed.
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
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netbind}}, \code{\link{netcomb}},
#'   \code{\link{forest.meta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Linde2016)
#' 
#' # Only consider studies including Face-to-face PST (to reduce
#' # runtime of example)
#' #
#' face <- subset(Linde2016, id %in% c(16, 24, 49, 118))
#' 
#' # Standard random effects NMA model (with placebo as reference
#' # treatment)
#' #
#' net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'                 data = face, reference.group = "placebo",
#'                 sm = "OR", comb.fixed = FALSE)
#' 
#' # Additive CNMA model with placebo as inactive component and
#' # reference
#' #
#' nc1 <- netcomb(net1, inactive = "placebo")
#' 
#' # Combine results of standard NMA and CNMA
#' #
#' nb1 <- netbind(nc1, net1,
#'                name = c("Additive CNMA", "Standard NMA"),
#'                col.study = c("red", "black"),
#'                col.square = c("red", "black"))
#' forest(nb1,
#'        col.by = "black", addrow.subgroups = FALSE,
#'        fontsize = 10, spacing = 0.7, squaresize = 0.9,
#'        label.left = "Favours Placebo",
#'        label.right = "Favours other")
#' 
#' @method forest netbind
#' @export


forest.netbind <- function(x,
                           pooled = ifelse(x$comb.random, "random", "fixed"),
                           ##
                           equal.size = TRUE,
                           ##
                           leftcols = "studlab",
                           leftlabs = "Treatment",
                           rightcols = c("effect", "ci"),
                           rightlabs = NULL,
                           ##
                           subset.treatments,
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
  meta:::chkclass(x, "netbind")
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
  chklogical(equal.size)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)
  chklogical(backtransf)
  ##
  chkchar(lab.NA)
  
  
  ##
  ##
  ## (2) Extract results for fixed and random effects model
  ##
  ##
  if (pooled == "fixed") {
    if (!missing(subset.treatments)) {
      subset.treatments <- setchar(subset.treatments, unique(x$fixed$treat))
      sel <- x$fixed$treat %in% subset.treatments
    }
    else
      sel <- x$fixed$treat != x$reference.group
    ##
    m <-
      suppressWarnings(metagen(x$fixed$TE, x$fixed$seTE,
                               studlab = x$fixed$name,
                               sm = x$sm,
                               comb.fixed = FALSE, comb.random = FALSE,
                               byvar = x$fixed$treat, print.byvar = FALSE,
                               subset = sel))
    ##
    m$studlab <- x$fixed$name[sel]
    m$TE <- x$fixed$TE[sel]
    m$seTE <- x$fixed$seTE[sel]
    m$lower <- x$fixed$lower[sel]
    m$upper <- x$fixed$upper[sel]
    m$statistic <- x$fixed$statistic[sel]
    m$pval <- x$fixed$pval[sel]
    m$zval <- x$fixed$statistic[sel]
    ##
    m$col.study <- x$fixed$col.study[sel]
    m$col.square <- x$fixed$col.square[sel]
    m$col.square.lines <- x$fixed$col.square.lines[sel]
    m$col.inside <- x$fixed$col.inside[sel]
    ##
    text.pooled <- "Fixed Effects Model"
  }
  else {
    if (!missing(subset.treatments)) {
      subset.treatments <- setchar(subset.treatments, unique(x$random$treat))
      sel <- x$random$treat %in% subset.treatments
    }
    else
      sel <- x$random$treat != x$reference.group
    ##
    m <-
      suppressWarnings(metagen(x$random$TE, x$random$seTE,
                               studlab = x$random$name,
                               sm = x$sm,
                               comb.fixed = FALSE, comb.random = FALSE,
                               byvar = x$random$treat, print.byvar = FALSE,
                               subset = sel))
    ##
    m$studlab <- x$random$name[sel]
    m$TE <- x$random$TE[sel]
    m$seTE <- x$random$seTE[sel]
    m$lower <- x$random$lower[sel]
    m$upper <- x$random$upper[sel]
    m$statistic <- x$random$statistic[sel]
    m$pval <- x$random$pval[sel]
    m$zval <- x$random$statistic[sel]
    ##
    m$col.study <- x$random$col.study[sel]
    m$col.square <- x$random$col.square[sel]
    m$col.square.lines <- x$random$col.square.lines[sel]
    m$col.inside <- x$random$col.inside[sel]
    ##
    text.pooled <- "Random Effects Model"
  }
  ##
  if (missing(smlab))
    if (x$baseline.reference)
      smlab <- paste0("Comparison: other vs '",
                      x$reference.group, "'\n(",
                      text.pooled,
                      ")")
    else
      smlab <- paste0("Comparison: '",
                      x$reference.group, "' vs other\n(",
                      text.pooled,
                      ")")
  
  
  ##
  ##
  ## (3) Forest plot
  ##
  ##
  forest(m,
         digits = digits,
         comb.fixed = FALSE, comb.random = FALSE,
         overall = FALSE, hetstat = FALSE, test.subgroup = FALSE,
         leftcols = leftcols,
         leftlabs = leftlabs,
         rightcols = rightcols,
         rightlabs = rightlabs,
         lab.NA = lab.NA,
         smlab = smlab,
         backtransf = backtransf,
         ##
         col.study = m$col.study,
         col.square = m$col.square,
         col.square.lines = m$col.square.lines,
         col.inside = m$col.inside,
         col.inside.fixed = "black",
         col.inside.random = "black",
         ##
         weight.study = if (equal.size) "same" else pooled,
         calcwidth.subgroup = TRUE,
         ...)
  
  
  invisible(NULL)
}





#' @rdname forest.netbind
#' @method plot netbind
#' @export
#'

plot.netbind <- function(x, ...)
  forest(x, ...)
