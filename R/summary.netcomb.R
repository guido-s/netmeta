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
#' @param digits.zval Minimal number of significant digits for z- or
#'   t-value, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param digits.pval.Q Minimal number of significant digits for
#'   p-value of heterogeneity tests, see \code{print.default}.
#' @param digits.Q Minimal number of significant digits for
#'   heterogeneity statistics, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance, see \code{print.default}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   statistic, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
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
#' print(summary(nc1), digits = 2, digits.zval = 3)
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
#' print(summary(nc2), digits = 2, digits.zval = 3)
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
  keepvars <- c("TE", "seTE", "lower", "upper", "z", "p")
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
                             z = object$zval.nma.fixed,
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
                              z = object$zval.cnma.fixed,
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
                              z = object$zval.nma.random,
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
                               z = object$zval.cnma.random,
                               p = object$pval.cnma.random,
                               stringsAsFactors = FALSE)
  ##
  ci.f <- list(TE = object$TE.fixed,
               seTE = object$seTE.fixed,
               lower = object$lower.fixed,
               upper = object$upper.fixed,
               z = object$zval.fixed,
               p = object$pval.fixed)
  ##
  ci.r <- list(TE = object$TE.random,
               seTE = object$seTE.random,
               lower = object$lower.random,
               upper = object$upper.random,
               z = object$zval.random,
               p = object$pval.random)
  ##
  ci.comp.f <- data.frame(TE = object$Comp.fixed,
                          seTE = object$seComp.fixed,
                          lower = object$lower.Comp.fixed,
                          upper = object$upper.Comp.fixed,
                          z = object$zval.Comp.fixed,
                          p = object$pval.Comp.fixed,
                          stringsAsFactors = FALSE)
  rownames(ci.comp.f) <- object$comps
  ##
  ci.comp.r <- data.frame(TE = object$Comp.random,
                          seTE = object$seComp.random,
                          lower = object$lower.Comp.random,
                          upper = object$upper.Comp.random,
                          z = object$zval.Comp.random,
                          p = object$pval.Comp.random,
                          stringsAsFactors = FALSE)
  rownames(ci.comp.r) <- object$comps
  ##
  ci.comb.f <- data.frame(TE = object$Comb.fixed,
                          seTE = object$seComb.fixed,
                          lower = object$lower.Comb.fixed,
                          upper = object$upper.Comb.fixed,
                          z = object$zval.Comb.fixed,
                          p = object$pval.Comb.fixed,
                          stringsAsFactors = FALSE)
  rownames(ci.comb.f) <- object$trts
  ##
  ci.comb.r <- data.frame(TE = object$Comb.random,
                          seTE = object$seComb.random,
                          lower = object$lower.Comb.random,
                          upper = object$upper.Comb.random,
                          z = object$zval.Comb.random,
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
  class(res) <- "summary.netcomb"
  
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
                                  digits.zval = gs("digits.zval"),
                                  digits.pval = gs("digits.pval"),
                                  digits.pval.Q = max(gs("digits.pval.Q"), 2),
                                  digits.Q = gs("digits.Q"),
                                  digits.tau2 = gs("digits.tau2"),
                                  digits.I2 = gs("digits.I2"),
                                  scientific.pval = gs("scientific.pval"),
                                  big.mark = gs("big.mark"),
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
  formatN <- meta:::formatN
  formatPT <- meta:::formatPT
  ##  
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(backtransf)
  chknumeric(nchar.trts, min = 1, single = TRUE)
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.zval, min = 0, single = TRUE)
  chknumeric(digits.pval, min = 1, single = TRUE)
  chknumeric(digits.pval.Q, min = 1, single = TRUE)
  chknumeric(digits.Q, min = 0, single = TRUE)
  chknumeric(digits.tau2, min = 0, single = TRUE)
  chknumeric(digits.I2, min = 0, single = TRUE)
  ##
  chklogical(scientific.pval)
  
  
  I2 <- round(100 * x$I2, digits.I2)
  
  
  if (comb.fixed | comb.random) {
    cat(paste("Number of studies: k = ", x$k, "\n", sep = ""))
    cat(paste("Number of treatments: n = ", x$n, "\n", sep = ""))
    cat(paste("Number of active components: c = ", x$c, "\n", sep = ""))
    cat(paste("Number of pairwise comparisons: m = ", x$m, "\n", sep = ""))
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
                     digits, digits.zval, digits.pval.Q,
                     scientific.pval, big.mark, x$seq)
  ##
  dat1.r <- formatCC(x$combinations.random,
                     backtransf, x$sm, x$level, trts.abbr,
                     digits, digits.zval, digits.pval.Q,
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
                     digits, digits.zval, digits.pval.Q,
                     scientific.pval, big.mark)
  ##
  dat2.r <- formatCC(x$components.random,
                     backtransf, x$sm, x$level, comps.abbr,
                     digits, digits.zval, digits.pval.Q,
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
  
  
  cat(paste("\nQuantifying heterogeneity / inconsistency:\n",
            formatPT(x$tau^2,
                     lab = TRUE, labval = "tau^2",
                     digits = digits.tau2,
                     lab.NA = "NA", big.mark = big.mark),
            if (!is.na(I2))
              paste("; I^2 = ", round(I2, digits.I2), "%", "", sep = ""),
            "\n", sep = ""))
  
  
  cat("\nHeterogeneity statistics:\n")
  
  print(data.frame(Q = formatN(c(x$Q.additive,
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
                   row.names = c("Additive model", "Standard model", "Difference")))
  
  
  if ((comb.fixed | comb.random)) {
    any.trts <- any(trts != trts.abbr)
    any.comps <- any(comps != comps.abbr)
    ##
    if (any.trts | any.comps)
      cat("\nLegend", if (any.trts & any.comps) "s", ":", sep = "")
    ##
    if (any.trts) {
      ##
      tmat <- data.frame(trts.abbr, trts)
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
