#' Print detailed information for component network meta-analysis
#' 
#' @description
#' Print detailed information  for component network meta-analysis.
#' 
#' @param x An object of class \code{summary.netcomb}
#' @param common A logical indicating whether results for the common
#'   effects model should be printed.
#' @param random A logical indicating whether results for the random
#'   effects model should be printed.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf=TRUE}, results for \code{sm="OR"} are presented
#'   as odds ratios rather than log odds ratios, for example.
#' @param nchar.comps A numeric defining the minimum number of
#'   characters used to create unique component names.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-value, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param digits.pval.Q Minimal number of significant digits for
#'   p-value of heterogeneity tests, see \code{print.default}.
#' @param digits.Q Minimal number of significant digits for
#'   heterogeneity statistics, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   effects should be printed according to JAMA reporting standards.
#' @param big.mark A character used as thousands separator.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netcomb}}, \code{\link{discomb}},
#'   \code{\link{summary.netcomb}}
#' 
#' @keywords print
#' 
#' @examples
#' data(Linde2016)
#' 
#' # Only consider studies including Face-to-face PST (to reduce
#' # runtime of example)
#' #
#' face <- subset(Linde2016, id %in% c(16, 24, 49, 118))
#' 
#' # Conduct random effects network meta-analysis
#' #
#' net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'   data = face, reference.group = "placebo",
#'   sm = "OR", common = FALSE)
#' 
#' # Additive model for treatment components
#' #
#' nc1 <- netcomb(net1)
#' print(summary(nc1), digits = 2)
#' 
#' @method print summary.netcomb
#' @export


print.summary.netcomb <- function(x,
                                  common = x$x$common,
                                  random = x$x$random,
                                  backtransf = x$backtransf,
                                  nchar.comps = x$nchar.comps,
                                  ##
                                  digits = gs("digits"),
                                  digits.stat = gs("digits.stat"),
                                  digits.pval = gs("digits.pval"),
                                  digits.pval.Q = max(gs("digits.pval.Q"), 2),
                                  digits.Q = gs("digits.Q"),
                                  ##
                                  scientific.pval = gs("scientific.pval"),
                                  zero.pval = gs("zero.pval"),
                                  JAMA.pval = gs("JAMA.pval"),
                                  ##
                                  big.mark = gs("big.mark"),
                                  ##
                                  legend = TRUE,
                                  ##
                                  warn.deprecated = gs("warn.deprecated"),
                                  ##
                                  ...) {
  
  
  ##
  ##
  ## (1) Check for summary.netcomb object and upgrade object
  ##
  ##
  chkclass(x, "summary.netcomb")
  updateversion(x)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.pval.Q, min = 1, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  ##
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  ##
  chklogical(legend)
  ##
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                      warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  ##
  nchar.comps <-
    deprecated(nchar.comps, missing(nchar.comps), args, "nchar.trts")
  nchar.comps <- replaceNULL(nchar.comps, 666)
  chknumeric(nchar.comps, min = 1, length = 1)
  
  
  ##
  ##
  ## (3) Print summary results
  ##
  ##
  comps <- sort(c(x$comps, x$inactive))
  comps.abbr <- treats(comps, nchar.comps)
  ##
  comp.f <- x$comparison.cnma.common
  comp.f$seTE <- NULL
  ##
  dat.f <- formatComp(comp.f,
                      backtransf, x$sm, x$level.ma,
                      comps, comps.abbr, x$sep.comps,
                      digits, digits.stat, digits.pval.Q,
                      scientific.pval, zero.pval, JAMA.pval,
                      big.mark)
  ##
  comp.r <- x$comparison.cnma.random
  comp.r$seTE <- NULL
  ##
  dat.r <- formatComp(comp.r,
                      backtransf, x$sm, x$level.ma,
                      comps, comps.abbr, x$sep.comps,
                      digits, digits.stat, digits.pval.Q,
                      scientific.pval, zero.pval, JAMA.pval,
                      big.mark)
  ##
  if (common) {
    cat("Additive model (common effects model):\n")
    prmatrix(dat.f, quote = FALSE, right = TRUE, ...)
    cat("\n")
  }
  ##
  if (random) {
    cat("Additive model (random effects model):\n")
    prmatrix(dat.r, quote = FALSE, right = TRUE, ...)
    cat("\n")
  }
  ##  
  if (common | random)
    print.netcomb(x$x,
                  common = common,
                  random = random,
                  backtransf = backtransf,
                  nchar.comps = nchar.comps,
                  ##
                  digits = digits,
                  digits.stat = digits.stat,
                  digits.pval = digits.pval,
                  digits.pval.Q = digits.pval.Q,
                  digits.Q = digits.Q,
                  ##
                  scientific.pval = scientific.pval,
                  zero.pval = zero.pval,
                  JAMA.pval = JAMA.pval,
                  ##
                  big.mark = big.mark,
                  ##
                  legend = legend,
                  ##
                  ...)
  else
    cat("Please use argument 'common = TRUE' or 'random = TRUE'",
        "to print network meta-analysis results.\n")
  
  
  invisible(NULL)
}
