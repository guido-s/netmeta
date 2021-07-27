#' Print objects of class netcomb
#' 
#' @description
#' Print method for objects of class \code{netcomb}.
#' 
#' @param x An object of class \code{netcomb}
#' @param comb.fixed A logical indicating whether results for the
#'   fixed effects (common effects) model should be printed.
#' @param comb.random A logical indicating whether results for the
#'   random effects model should be printed.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf=TRUE}, results for \code{sm="OR"} are presented
#'   as odds ratios rather than log odds ratios, for example.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
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
#' @param big.mark A character used as thousands separator.
#' @param \dots Additional arguments.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
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
#'                 data = face, reference.group = "placebo",
#'                 sm = "OR", comb.fixed = FALSE)
#' 
#' # Additive model for treatment components
#' #
#' nc1 <- netcomb(net1)
#' print(nc1, digits = 2)
#' 
#' @method print netcomb
#' @export
#' @export print.netcomb


print.netcomb <- function(x,
                          comb.fixed = x$comb.fixed,
                          comb.random = x$comb.random,
                          backtransf = x$backtransf,
                          nchar.trts = x$nchar.trts,
                          ##
                          digits = gs("digits"),
                          digits.stat = gs("digits.stat"),
                          digits.pval = gs("digits.pval"),
                          digits.pval.Q = max(gs("digits.pval.Q"), 2),
                          digits.Q = gs("digits.Q"),
                          scientific.pval = gs("scientific.pval"),
                          big.mark = gs("big.mark"),
                          ...) {
  
  
  meta:::chkclass(x, "netcomb")
  ##  
  x <- upgradenetmeta(x)
  
  
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  meta:::chklogical(backtransf)
  meta:::chknumeric(nchar.trts, min = 1, length = 1)
  ##
  meta:::chknumeric(digits, min = 0, length = 1)
  meta:::chknumeric(digits.stat, min = 0, length = 1)
  meta:::chknumeric(digits.pval, min = 1, length = 1)
  meta:::chknumeric(digits.pval.Q, min = 1, length = 1)
  meta:::chknumeric(digits.Q, min = 0, length = 1)
  ##
  meta:::chklogical(scientific.pval)
  
  
  trts <- x$trts
  trts.abbr <- treats(trts, nchar.trts)
  ##
  cnma.f <- data.frame(studlab = x$studlab,
                       treat1 = x$treat1,
                       treat2 = x$treat2,
                       TE = x$TE.cnma.fixed,
                       lower = x$lower.cnma.fixed,
                       upper = x$upper.cnma.fixed,
                       statistic = x$statistic.cnma.fixed,
                       p = x$pval.cnma.fixed,
                       stringsAsFactors = FALSE)
  ##
  dat.f <- formatComp(cnma.f,
                      backtransf, x$sm, x$level.comb,
                      trts, trts.abbr,
                      digits, digits.stat, digits.pval.Q,
                      scientific.pval, big.mark)
  ##
  cnma.r <- data.frame(studlab = x$studlab,
                       treat1 = x$treat1,
                       treat2 = x$treat2,
                       TE = x$TE.cnma.random,
                       lower = x$lower.cnma.random,
                       upper = x$upper.cnma.random,
                       statistic = x$statistic.cnma.random,
                       p = x$pval.cnma.random,
                       stringsAsFactors = FALSE)
  ##
  dat.r <- formatComp(cnma.r,
                      backtransf, x$sm, x$level.comb,
                      trts, trts.abbr,
                      digits, digits.stat, digits.pval.Q,
                      scientific.pval, big.mark)
  ##
  if (comb.fixed) {
    cat("Additive model (fixed effects model):\n")
    prmatrix(dat.f, quote = FALSE, right = TRUE, ...)
    cat("\n")
  }
  ##
  if (comb.random) {
    cat("Additive model (random effects model):\n")
    prmatrix(dat.r, quote = FALSE, right = TRUE, ...)
    cat("\n")
  }
  

  if (comb.fixed | comb.random)
    print(summary(x),
          comb.fixed = comb.fixed,
          comb.random = comb.random,
          backtransf = backtransf,
          nchar.trts = nchar.trts,
          ##
          digits = digits,
          digits.stat = digits.stat,
          digits.pval = digits.pval,
          digits.pval.Q = digits.pval.Q,
          digits.Q = digits.Q,
          scientific.pval = scientific.pval,
          big.mark = big.mark, ...)
  else
    cat("Please use argument 'comb.fixed = TRUE' or",
        "'comb.random = TRUE' to print meta-analysis results.\n",
        sep = "")
  
  
  invisible(NULL)
}
