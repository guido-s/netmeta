#' Print method for objects of class netcomb
#' 
#' @description
#' Print method for objects of class \code{netcomb}.
#' 
#' @param x An object of class \code{netcomb} or
#'   \code{summary.netcomb}.
#' @param fixed A logical indicating whether results for the
#'   fixed effects (common effects) model should be printed.
#' @param random A logical indicating whether results for the
#'   random effects model should be printed.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param nchar.comps A numeric defining the minimum number of
#'   characters used to create unique names for components.
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
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for square
#'   root of between-study variance, see \code{print.default}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   statistic, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   component or combination effect should be printed according to
#'   JAMA reporting standards.
#' @param big.mark A character used as thousands separator.
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param text.I2 Text printed to identify heterogeneity statistic
#'   I\eqn{^2}.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#'
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netcomb}}, \code{\link{discomb}}
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
#'                 sm = "OR", fixed = FALSE)
#' 
#' # Additive model for treatment components
#' #
#' nc1 <- netcomb(net1)
#' nc1
#' print(nc1, digits = 2, digits.stat = 3)
#' 
#' \dontrun{
#' # Conduct random effects network meta-analysis
#' #
#' net2 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'                 data = Linde2016, reference.group = "placebo",
#'                 sm = "OR", fixed = FALSE)
#' 
#' # Additive model for treatment components
#' #
#' nc2 <- netcomb(net2)
#' nc2
#' print(nc2, digits = 2, digits.stat = 3)
#' }
#' 
#' @method print netcomb
#' @export


print.netcomb <- function(x,
                          fixed = x$fixed,
                          random = x$random,
                          backtransf = x$backtransf,
                          nchar.comps = x$nchar.comps,
                          ##
                          digits = gs("digits"),
                          digits.stat = gs("digits.stat"),
                          digits.pval = gs("digits.pval"),
                          digits.pval.Q = max(gs("digits.pval.Q"), 2),
                          digits.Q = gs("digits.Q"),
                          digits.tau2 = gs("digits.tau2"),
                          digits.tau = gs("digits.tau"),
                          digits.I2 = gs("digits.I2"),
                          ##
                          scientific.pval = gs("scientific.pval"),
                          zero.pval = gs("zero.pval"),
                          JAMA.pval = gs("JAMA.pval"),
                          ##
                          big.mark = gs("big.mark"),
                          ##
                          text.tau2 = gs("text.tau2"),
                          text.tau = gs("text.tau"),
                          text.I2 = gs("text.I2"),
                          ##
                          legend = TRUE,
                          ##
                          warn.deprecated = gs("warn.deprecated"),
                          ##
                          ...) {
  
  
  ##
  ##
  ## (1) Check for netcomb object and upgrade object
  ##
  ##
  chkclass(x, "netcomb")
  x <- updateversion(x)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chklogical(backtransf)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.pval.Q, min = 1, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  ##
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  ##
  chkchar(text.tau2)
  chkchar(text.tau)
  chkchar(text.I2)
  ##
  chklogical(legend)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  fixed <- deprecated(fixed, missing(fixed), args, "comb.fixed",
                      warn.deprecated)
  chklogical(fixed)
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  ##
  nchar.comps <-
    deprecated(nchar.comps, missing(nchar.comps), args, "nchar.trts")
  
  
  ##
  ##
  ## (3) Print results
  ##
  ##
  I2 <- round(100 * x$I2, digits.I2)
  lower.I2 <- round(100 * x$lower.I2, digits.I2)
  upper.I2 <- round(100 * x$upper.I2, digits.I2)
  ##  
  if (fixed | random) {
    cat(paste("Number of studies: k = ", x$k, "\n", sep = ""))
    cat(paste("Number of treatments: n = ", x$n, "\n", sep = ""))
    cat(paste("Number of active components: c = ", x$c, "\n", sep = ""))
    cat(paste("Number of pairwise comparisons: m = ", x$m, "\n", sep = ""))
    if (!is.null(x$d))
      cat(paste("Number of designs: d = ", x$d, "\n", sep = ""))
    if (inherits(x, "discomb"))
      cat(paste("Number of subnetworks: s = ", x$s, "\n", sep = ""))
    ##
    cat("\n")
  }
  ##  
  comps <- sort(c(x$comps, x$inactive))
  comps.abbr <- treats(comps, nchar.comps)
  ##  
  ci.comb.f <- data.frame(TE = x$Comb.fixed,
                          seTE = x$seComb.fixed,
                          lower = x$lower.Comb.fixed,
                          upper = x$upper.Comb.fixed,
                          statistic = x$statistic.Comb.fixed,
                          p = x$pval.Comb.fixed,
                          stringsAsFactors = FALSE)
  rownames(ci.comb.f) <- x$trts
  ##
  dat1.f <- formatCC(ci.comb.f,
                     backtransf, x$sm, x$level,
                     comps, comps.abbr, x$sep.comps,
                     digits, digits.stat, digits.pval,
                     scientific.pval, zero.pval, JAMA.pval,
                     big.mark,
                     x$seq)
  ##
  ci.comb.r <- data.frame(TE = x$Comb.random,
                          seTE = x$seComb.random,
                          lower = x$lower.Comb.random,
                          upper = x$upper.Comb.random,
                          statistic = x$statistic.Comb.random,
                          p = x$pval.Comb.random,
                          stringsAsFactors = FALSE)
  rownames(ci.comb.r) <- x$trts
  ##
  dat1.r <- formatCC(ci.comb.r,
                     backtransf, x$sm, x$level,
                     comps, comps.abbr, x$sep.comps,
                     digits, digits.stat, digits.pval,
                     scientific.pval, zero.pval, JAMA.pval,
                     big.mark,
                     x$seq)
  ##
  if (fixed) {
    cat("Results for combinations (additive model, fixed effects model):\n")
    print(dat1.f)
    cat("\n")
  }
  ##
  if (random) {
    cat("Results for combinations (additive model, random effects model):\n")
    print(dat1.r)
    cat("\n")
  }
  ##  
  ci.comp.f <- data.frame(TE = x$Comp.fixed,
                          seTE = x$seComp.fixed,
                          lower = x$lower.Comp.fixed,
                          upper = x$upper.Comp.fixed,
                          statistic = x$statistic.Comp.fixed,
                          p = x$pval.Comp.fixed,
                          stringsAsFactors = FALSE)
  rownames(ci.comp.f) <- x$comps
  ##
  dat2.f <- formatCC(ci.comp.f,
                     backtransf, x$sm, x$level,
                     comps, comps.abbr, x$sep.comps,
                     digits, digits.stat, digits.pval,
                     scientific.pval, zero.pval, JAMA.pval,
                     big.mark)
  ##
  ci.comp.r <- data.frame(TE = x$Comp.random,
                          seTE = x$seComp.random,
                          lower = x$lower.Comp.random,
                          upper = x$upper.Comp.random,
                          statistic = x$statistic.Comp.random,
                          p = x$pval.Comp.random,
                          stringsAsFactors = FALSE)
  rownames(ci.comp.r) <- x$comps
  ##
  dat2.r <- formatCC(ci.comp.r,
                     backtransf, x$sm, x$level,
                     comps, comps.abbr, x$sep.comps,
                     digits, digits.stat, digits.pval,
                     scientific.pval, zero.pval, JAMA.pval,
                     big.mark)
  ##
  if (fixed) {
    cat("Results for components (fixed effects model):\n")
    print(dat2.f)
    cat("\n")
  }
  ##
  if (random) {
    cat("Results for components (random effects model):\n")
    print(dat2.r)
  }
  ##  
  cat(paste0("\nQuantifying heterogeneity / inconsistency:\n",
             formatPT(x$tau^2,
                      lab = TRUE, labval = text.tau2,
                      digits = digits.tau2,
                      lab.NA = "NA", big.mark = big.mark),
             "; ",
             formatPT(x$tau,
                      lab = TRUE, labval = text.tau,
                      digits = digits.tau,
                      lab.NA = "NA", big.mark = big.mark),
             if (!is.na(I2))
               paste0("; ", text.I2, " = ", round(I2, digits.I2), "%"),
             if (!(is.na(lower.I2) | is.na(upper.I2)))
               pasteCI(lower.I2, upper.I2, digits.I2, big.mark, unit = "%"),
             "\n")
      )
  ##  
  cat("\nHeterogeneity statistics:\n")
  ##
  hetdat <- 
    data.frame(Q = formatN(c(x$Q.additive,
                             x$Q.standard,
                             x$Q.diff),
                           digits.Q),
               df.Q = formatN(c(x$df.Q.additive,
                                x$df.Q.standard,
                                x$df.Q.diff), 0),
               pval = formatPT(c(x$pval.Q.additive,
                                 x$pval.Q.standard,
                                 x$pval.Q.diff),
                               digits = digits.pval.Q,
                               scientific = scientific.pval),
               row.names = c("Additive model", "Standard model",
                             "Difference"))
  ##
  names(hetdat) <- c("Q", "df", "p-value")
  ##
  print(hetdat)
  ##  
  if (legend && (fixed | random)) {
    diff.comps <- comps != comps.abbr
    if (any(diff.comps)) {
      tmat <- data.frame(comps.abbr, comps)
      tmat <- tmat[diff.comps, ]
      names(tmat) <- c("Abbreviation", " Component name")
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat("\nLegend:\n")
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(comps.abbr))) 
    }
  }
  
  invisible(NULL)
}
