#' Summary method for objects of class netcomb
#' 
#' @description
#' Summary method for objects of class \code{netcomb}.
#' 
#' @param x An object of class \code{netcomb} or
#'   \code{summary.netcomb}.
#' @param object An object of class \code{netcomb}.
#' @param comb.fixed A logical indicating whether results for the
#'   fixed effects (common effects) model should be printed.
#' @param comb.random A logical indicating whether results for the
#'   random effects model should be printed.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
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
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for square
#'   root of between-study variance, see \code{print.default}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   statistic, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param text.I2 Text printed to identify heterogeneity statistic
#'   I\eqn{^2}.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param \dots Additional arguments.
#'
#' @return
#' A list is returned with the same elements as a
#' \code{\link{netcomb}} object.
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
#'                 sm = "OR", comb.fixed = FALSE)
#' 
#' # Additive model for treatment components
#' #
#' nc1 <- netcomb(net1)
#' summary(nc1)
#' print(summary(nc1), digits = 2, digits.stat = 3)
#' 
#' \dontrun{
#' # Conduct random effects network meta-analysis
#' #
#' net2 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'                 data = Linde2016, reference.group = "placebo",
#'                 sm = "OR", comb.fixed = FALSE)
#' 
#' # Additive model for treatment components
#' #
#' nc2 <- netcomb(net2)
#' summary(nc2)
#' print(summary(nc2), digits = 2, digits.stat = 3)
#' }
#' 
#' @rdname summary.netcomb
#' @method summary netcomb
#' @export
#' @export summary.netcomb


summary.netcomb <- function(object,
                            comb.fixed = object$comb.fixed,
                            comb.random = object$comb.random,
                            ...) {
  
  ##
  ##
  ## (1) Check for netcomb object
  ##
  ##
  meta:::chkclass(object, "netcomb")
  ##  
  object <- upgradenetmeta(object)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  
  
  ##
  ##
  ## (3) Summarise results for individual studies and network
  ##     meta-analyses
  ##
  ##
  keepvars <- c("TE", "seTE", "lower", "upper", "statistic", "p")
  ##
  ci.comp <- data.frame(studlab = object$studlab,
                        treat1 = object$treat1, treat2 = object$treat2,
                        ci(object$TE, object$seTE, object$level)[keepvars],
                        stringsAsFactors = FALSE)
  ##
  ci.nma.fixed <- data.frame(studlab = object$studlab,
                             treat1 = object$treat1,
                             treat2 = object$treat2,
                             TE = object$TE.nma.fixed,
                             seTE = object$seTE.nma.fixed,
                             lower = object$lower.nma.fixed,
                             upper = object$upper.nma.fixed,
                             statistic = object$statistic.nma.fixed,
                             p = object$pval.nma.fixed,
                             stringsAsFactors = FALSE)
  ##
  ci.cnma.fixed <- data.frame(studlab = object$studlab,
                              treat1 = object$treat1,
                              treat2 = object$treat2,
                              TE = object$TE.cnma.fixed,
                              seTE = object$seTE.cnma.fixed,
                              lower = object$lower.cnma.fixed,
                              upper = object$upper.cnma.fixed,
                              statistic = object$statistic.cnma.fixed,
                              p = object$pval.cnma.fixed,
                              stringsAsFactors = FALSE)
  ##
  ci.nma.random <- data.frame(studlab = object$studlab,
                              treat1 = object$treat1,
                              treat2 = object$treat2,
                              TE = object$TE.nma.random,
                              seTE = object$seTE.nma.random,
                              lower = object$lower.nma.random,
                              upper = object$upper.nma.random,
                              statistic = object$statistic.nma.random,
                              p = object$pval.nma.random,
                              stringsAsFactors = FALSE)
  ##
  ci.cnma.random <- data.frame(studlab = object$studlab,
                               treat1 = object$treat1,
                               treat2 = object$treat2,
                               TE = object$TE.cnma.random,
                               seTE = object$seTE.cnma.random,
                               lower = object$lower.cnma.random,
                               upper = object$upper.cnma.random,
                               statistic = object$statistic.cnma.random,
                               p = object$pval.cnma.random,
                               stringsAsFactors = FALSE)
  ##
  ci.f <- list(TE = object$TE.fixed,
               seTE = object$seTE.fixed,
               lower = object$lower.fixed,
               upper = object$upper.fixed,
               statistic = object$statistic.fixed,
               p = object$pval.fixed)
  ##
  ci.r <- list(TE = object$TE.random,
               seTE = object$seTE.random,
               lower = object$lower.random,
               upper = object$upper.random,
               statistic = object$statistic.random,
               p = object$pval.random)
  ##
  ci.comp.f <- data.frame(TE = object$Comp.fixed,
                          seTE = object$seComp.fixed,
                          lower = object$lower.Comp.fixed,
                          upper = object$upper.Comp.fixed,
                          statistic = object$statistic.Comp.fixed,
                          p = object$pval.Comp.fixed,
                          stringsAsFactors = FALSE)
  rownames(ci.comp.f) <- object$comps
  ##
  ci.comp.r <- data.frame(TE = object$Comp.random,
                          seTE = object$seComp.random,
                          lower = object$lower.Comp.random,
                          upper = object$upper.Comp.random,
                          statistic = object$statistic.Comp.random,
                          p = object$pval.Comp.random,
                          stringsAsFactors = FALSE)
  rownames(ci.comp.r) <- object$comps
  ##
  ci.comb.f <- data.frame(TE = object$Comb.fixed,
                          seTE = object$seComb.fixed,
                          lower = object$lower.Comb.fixed,
                          upper = object$upper.Comb.fixed,
                          statistic = object$statistic.Comb.fixed,
                          p = object$pval.Comb.fixed,
                          stringsAsFactors = FALSE)
  rownames(ci.comb.f) <- object$trts
  ##
  ci.comb.r <- data.frame(TE = object$Comb.random,
                          seTE = object$seComb.random,
                          lower = object$lower.Comb.random,
                          upper = object$upper.Comb.random,
                          statistic = object$statistic.Comb.random,
                          p = object$pval.Comb.random,
                          stringsAsFactors = FALSE)
  rownames(ci.comb.r) <- object$trts
  
  
  ##
  ##
  ## (4) Create summary.netmeta object
  ##
  ##
  res <- list(k = object$k,
              m = object$m,
              n = object$n,
              d = object$d,
              c = object$c,
              s = object$s,
              ##
              trts = object$trts,
              k.trts = object$k.trts,
              n.trts = object$n.trts,
              events.trts = object$events.trts,
              ##
              studies = object$studies,
              narms = object$narms,
              ##
              designs = object$designs,
              ##
              comps = object$comps,
              k.comps = object$k.comps,
              n.comps = object$n.comps,
              events.comps = object$events.comps,
              ##
              comparison = ci.comp,
              comparison.nma.fixed = ci.nma.fixed,
              comparison.nma.random = ci.nma.random,
              comparison.cnma.fixed = ci.cnma.fixed,
              comparison.cnma.random = ci.cnma.random,
              ##
              components.fixed = ci.comp.f,
              components.random = ci.comp.r,
              ##
              combinations.fixed = ci.comb.f,
              combinations.random = ci.comb.r,
              ##
              fixed = ci.f, random = ci.r,
              ##
              Q.additive = object$Q.additive,
              df.Q.additive = object$df.Q.additive,
              pval.Q.additive = object$pval.Q.additive,
              tau = object$tau,
              I2 = object$I2,
              lower.I2 = object$lower.I2, upper.I2 = object$upper.I2,
              ##
              Q.standard = object$Q.standard,
              df.Q.standard = object$df.Q.standard,
              pval.Q.standard = object$pval.Q.standard,
              ##
              Q.diff = object$Q.diff,
              df.Q.diff = object$df.Q.diff,
              pval.Q.diff = object$pval.Q.diff, 
              ##
              sm = object$sm,
              method = object$method,
              level = object$level,
              level.comb = object$level.comb,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              ##
              ci.lab = paste0(round(100 * object$level.comb, 1),"%-CI"),
              ##
              reference.group = object$reference.group,
              baseline.reference = object$baseline.reference,
              all.treatments = object$all.treatments,
              seq = object$seq,
              ##
              tau.preset = object$tau.preset,
              ##
              sep.trts = object$sep.trts,
              nchar.trts = object$nchar.trts,
              ##
              inactive = object$inactive,
              sep.comps = object$sep.comps,
              ##
              backtransf = object$backtransf,
              ##
              title = object$title,
              ##
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- c(if (inherits(object, "discomb")) "summary.discomb",
                  "summary.netcomb")
  
  res
}





#' @rdname summary.netcomb
#' @method print summary.netcomb
#' @export
#' @export print.summary.netcomb


print.summary.netcomb <- function(x,
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
                                  digits.tau2 = gs("digits.tau2"),
                                  digits.tau = gs("digits.tau"),
                                  digits.I2 = gs("digits.I2"),
                                  scientific.pval = gs("scientific.pval"),
                                  big.mark = gs("big.mark"),
                                  ##
                                  text.tau2 = gs("text.tau2"),
                                  text.tau = gs("text.tau"),
                                  text.I2 = gs("text.I2"),
                                  ##
                                  legend = TRUE,
                                  ...) {
  
  
  ##
  ##
  ## (1) Check class and arguments
  ##
  ##
  meta:::chkclass(x, "summary.netcomb")
  ##  
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  chkchar <- meta:::chkchar
  formatN <- meta:::formatN
  formatPT <- meta:::formatPT
  pasteCI <- meta:::pasteCI
  ##  
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(backtransf)
  chknumeric(nchar.trts, min = 1, length = 1)
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
  ##
  chkchar(text.tau2)
  chkchar(text.tau)
  chkchar(text.I2)
  ##
  chklogical(legend)
  
  
  I2 <- round(100 * x$I2, digits.I2)
  lower.I2 <- round(100 * x$lower.I2, digits.I2)
  upper.I2 <- round(100 * x$upper.I2, digits.I2)
  
  
  if (comb.fixed | comb.random) {
    cat(paste("Number of studies: k = ", x$k, "\n", sep = ""))
    cat(paste("Number of treatments: n = ", x$n, "\n", sep = ""))
    cat(paste("Number of active components: c = ", x$c, "\n", sep = ""))
    cat(paste("Number of pairwise comparisons: m = ", x$m, "\n", sep = ""))
    if (!is.null(x$d))
      cat(paste("Number of designs: d = ", x$d, "\n", sep = ""))
    if (inherits(x, "summary.discomb"))
      cat(paste("Number of subnetworks: s = ", x$s, "\n", sep = ""))
    ##
    cat("\n")
  }
  
  
  trts <- x$trts
  trts.abbr <- treats(trts, nchar.trts)
  ##
  comps <- x$comps
  comps.abbr <- treats(comps, nchar.trts)
  
  
  dat1.f <- formatCC(x$combinations.fixed,
                     backtransf, x$sm, x$level, trts.abbr,
                     digits, digits.stat, digits.pval,
                     scientific.pval, big.mark, x$seq)
  ##
  dat1.r <- formatCC(x$combinations.random,
                     backtransf, x$sm, x$level, trts.abbr,
                     digits, digits.stat, digits.pval,
                     scientific.pval, big.mark, x$seq)
  ##
  if (comb.fixed) {
    cat("Results for combinations (additive model, fixed effects model):\n")
    print(dat1.f)
    cat("\n")
  }
  ##
  if (comb.random) {
    cat("Results for combinations (additive model, random effects model):\n")
    print(dat1.r)
    cat("\n")
  }
  
  
  dat2.f <- formatCC(x$components.fixed,
                     backtransf, x$sm, x$level, comps.abbr,
                     digits, digits.stat, digits.pval,
                     scientific.pval, big.mark)
  ##
  dat2.r <- formatCC(x$components.random,
                     backtransf, x$sm, x$level, comps.abbr,
                     digits, digits.stat, digits.pval,
                     scientific.pval, big.mark)
  ##
  if (comb.fixed) {
    cat("Results for components (fixed effects model):\n")
    print(dat2.f)
    cat("\n")
  }
  ##
  if (comb.random) {
    cat("Results for components (random effects model):\n")
    print(dat2.r)
  }
  
  
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
  
  
  cat("\nHeterogeneity statistics:\n")

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
  
  
  if (legend && (comb.fixed | comb.random)) {
    diff.trts <- trts != trts.abbr
    any.trts <- any(diff.trts)
    diff.comps <- comps != comps.abbr
    any.comps <- any(diff.comps)
    ##
    if (any.trts | any.comps)
      cat("\nLegend", if (any.trts & any.comps) "s", ":", sep = "")
    ##
    if (any.trts) {
      ##
      tmat <- data.frame(trts.abbr, trts)
      tmat <- tmat[diff.trts, ]
      names(tmat) <- c("Abbreviation", "Treatment name")
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat("\n")
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(trts.abbr))) 
    }
    ##
    if (any.comps) {
      ##
      tmat <- data.frame(comps.abbr, comps)
      tmat <- tmat[diff.comps, ]
      names(tmat) <- c("Abbreviation", " Component name")
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat("\n")
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(comps.abbr))) 
    }
  }
  
  
  invisible(NULL)
}
