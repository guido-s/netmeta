#' Print and summary method for objects of class netmeta
#' 
#' @description
#' Print and summary method for objects of class \code{netmeta}.
#' 
#' @param x An object of class \code{summary.netmeta}.
#' @param object An object of class \code{netmeta}.
#' @param comb.fixed A logical indicating whether results for the
#'   fixed effects (common effects) model should be printed.
#' @param comb.random A logical indicating whether results for the
#'   random effects model should be printed.
#' @param prediction A logical indicating whether prediction intervals
#'   should be printed.
#' @param reference.group Reference treatment.
#' @param baseline.reference A logical indicating whether results
#'   should be expressed as comparisons of other treatments versus the
#'   reference treatment (default) or vice versa. This argument is
#'   only considered if \code{reference.group} has been specified.
#' @param all.treatments A logical or \code{"NULL"}. If \code{TRUE},
#'   matrices with all treatment effects, and confidence limits will
#'   be printed.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param header A logical indicating whether information on title of
#'   meta-analysis, comparison and outcome should be printed at the
#'   beginning of the printout.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
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
#' @param text.tau Text printed to identify \eqn{\tau}, the square root
#'   of the between-study variance \eqn{\tau^2}.
#' @param text.I2 Text printed to identify heterogeneity statistic
#'   I\eqn{^2}.
#' @param \dots Additional arguments.
#' 
#' @return
#'
#' A list is returned with the following elements:
#' \item{comparison}{Results for pairwise comparisons (data frame with
#'   columns studlab, treat1, treat2, TE, seTE, lower, upper, z, p).}
#' \item{comparison.nma.fixed}{Results for pairwise comparisons based
#'   on fixed effects model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p, leverage).}
#' \item{comparison.nma.random}{Results for pairwise comparisons based
#'   on random effects model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p).}
#' \item{fixed}{Results for fixed effects model (a list with elements
#'   TE, seTE, lower, upper, z, p).}
#' \item{random}{Results for random effects model (a list with
#'   elements TE, seTE, lower, upper, z, p).}
#' \item{predict}{Prediction intervals (a list with elements seTE,
#'   lower, upper).}
#' \item{studies}{Study labels coerced into a factor with its levels
#'   sorted alphabetically.}
#' \item{narms}{Number of arms for each study.}
#' \item{k}{Total number of studies.}
#' \item{m}{Total number of pairwise comparisons.}
#' \item{n}{Total number of treatments.}
#' \item{d}{Total number of designs (corresponding to the unique set
#'   of treatments compared within studies).}
#' \item{Q}{Overall heterogeneity / inconsistency statistic.}
#' \item{df.Q}{Degrees of freedom for test of heterogeneity /
#'   inconsistency.}
#' \item{pval.Q}{P-value for test of heterogeneity / inconsistency.}
#' \item{I2, lower.I2, upper.I2}{I-squared, lower and upper confidence
#'   limits.}
#' \item{tau}{Square-root of between-study variance.}
#' \item{Q.heterogeneity}{Overall heterogeneity statistic.}
#' \item{df.Q.heterogeneity}{Degrees of freedom for test of overall
#'   heterogeneity.}
#' \item{pval.Q.heterogeneity}{P-value for test of overall
#'   heterogeneity.}
#' \item{Q.inconsistency}{Overall inconsistency statistic.}
#' \item{df.Q.inconsistency}{Degrees of freedom for test of overall
#'   inconsistency.}
#' \item{pval.Q.inconsistency}{P-value for test of overall
#'   inconsistency.}
#' \item{sm}{A character string indicating underlying summary
#'   measure.}
#' \item{method}{A character string indicating which method is to be
#'   used for pooling of studies.}
#' \item{level}{The level used to calculate confidence intervals for
#'   individual studies.}
#' \item{level.comb}{The level used to calculate confidence intervals
#'   for pooled estimates.}
#' \item{comb.fixed, comb.random}{As defined above.}
#' \item{prediction, level.predict}{As defined above.}
#' \item{reference.group, baseline.reference}{As defined above.}
#' \item{all.treatments, backtransf}{As defined above.}
#' \item{ci.lab}{Label for confidence interval.}
#' \item{seq}{A character specifying the sequence of treatments.}
#' \item{tau.preset}{An optional value for the square-root of the
#'   between-study variance \eqn{\tau^2}.}
#' \item{sep.trts}{A character used in comparison names as separator
#'   between treatment labels.}
#' \item{nchar.trts}{A numeric defining the minimum number of
#'   characters used to create unique treatment names.}
#' \item{title}{Title of meta-analysis / systematic review.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}
#' 
#' @keywords print
#' 
#' @examples
#' data(Senn2013)
#' 
#' # Conduct fixed effects network meta-analysis
#' #
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'                 data = Senn2013, sm = "MD",
#'                 comb.random = FALSE)
#' print(net1, ref = "plac", digits = 3)
#' summary(net1)
#'
#' \dontrun{
#' # Conduct random effects network meta-analysis
#' #
#' net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'                 data = Senn2013, sm = "MD",
#'                 comb.fixed = FALSE)
#' print(net2, ref = "plac", digits = 3)
#' summary(net2)
#' }
#' 
#' @rdname summary.netmeta
#' @method summary netmeta
#' @export
#' @export summary.netmeta


summary.netmeta <- function(object,
                            comb.fixed = object$comb.fixed,
                            comb.random = object$comb.random,
                            prediction = object$prediction,
                            reference.group = object$reference.group,
                            baseline.reference = object$baseline.reference,
                            all.treatments = object$all.treatments,
                            ...) {
  
  ##
  ##
  ## (1) Check for netmeta object and upgrade older meta objects
  ##
  ##
  meta:::chkclass(object, "netmeta")
  ##
  is.bin <- inherits(object, "netmetabin")
  ##  
  object <- upgradenetmeta(object)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  meta:::chklogical(prediction)
  meta:::chklogical(baseline.reference)
  ##
  cl <- "netmeta()"
  addargs <- names(list(...))
  ##
  fun <- "summary.netmeta"
  ##
  meta:::warnarg("level", addargs, fun, cl)
  meta:::warnarg("level.comb", addargs, fun, cl)
  
  
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
                             TE = if (!is.bin) object$TE.nma.fixed else NA,
                             seTE = if (!is.bin) object$seTE.nma.fixed else NA,
                             lower = if (!is.bin) object$lower.nma.fixed else NA,
                             upper = if (!is.bin) object$upper.nma.fixed else NA,
                             z = if (!is.bin) object$zval.nma.fixed else NA,
                             p = if (!is.bin) object$pval.nma.fixed else NA,
                             leverage = if (!is.bin) object$leverage.fixed else NA,
                             stringsAsFactors = FALSE)
  ##
  ci.nma.random <- data.frame(studlab = object$studlab,
                              treat1 = object$treat1,
                              treat2 = object$treat2,
                              TE = if (!is.bin) object$TE.nma.random else NA,
                              seTE = if (!is.bin) object$seTE.nma.random else NA,
                              lower = if (!is.bin) object$lower.nma.random else NA,
                              upper = if (!is.bin) object$upper.nma.random else NA,
                              z = if (!is.bin) object$zval.nma.random else NA,
                              p = if (!is.bin) object$pval.nma.random else NA,
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
  ci.p <- list(seTE = object$seTE.predict,
               lower = object$lower.predict,
               upper = object$upper.predict)
  
  
  ##
  ##
  ## (4) Create summary.netmeta object
  ##
  ##
  res <- list(comparison = ci.comp,
              comparison.nma.fixed = ci.nma.fixed,
              comparison.nma.random = ci.nma.random,
              fixed = ci.f,
              random = ci.r,
              predict = ci.p,
              ##
              studies = object$studies,
              narms = object$narms,
              ##
              k = object$k, m = object$m, n = object$n, d = object$d,
              ##
              Q = object$Q,
              df.Q = object$df.Q,
              pval.Q = object$pval.Q,
              I2 = object$I2,
              lower.I2 = object$lower.I2, upper.I2 = object$upper.I2,
              tau = object$tau,
              ##
              Q.heterogeneity = object$Q.heterogeneity,
              df.Q.heterogeneity = object$df.Q.heterogeneity,
              pval.Q.heterogeneity = object$pval.Q.heterogeneity,
              ##
              Q.inconsistency = object$Q.inconsistency,
              df.Q.inconsistency = object$df.Q.inconsistency,
              pval.Q.inconsistency = object$pval.Q.inconsistency,
              ##
              Q.decomp = object$Q.decomp,
              ##
              sm = object$sm,
              method = object$method,
              level = object$level,
              level.comb = object$level.comb,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              ##
              prediction = prediction,
              level.predict = object$level.predict,
              ##
              incr = object$incr,
              allincr = object$allincr,
              addincr = object$addincr,
              allstudies = object$allstudies,
              cc.pooled = object$cc.pooled,
              ##
              ci.lab = paste0(round(100 * object$level.comb, 1),"%-CI"),
              ##
              reference.group = NA,
              baseline.reference = NA,
              all.treatments = NA,
              seq = object$seq,
              ##
              tau.preset = object$tau.preset,
              ##
              sep.trts = object$sep.trts,
              nchar.trts = object$nchar.trts,
              ##
              backtransf = object$backtransf,
              ##
              title = object$title,
              ##
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  if (reference.group != "" & missing(all.treatments))
    all.treatments <- FALSE
  ##
  if (reference.group != "")
    reference.group <- setref(reference.group, rownames(object$A.matrix))
  ##
  res$reference.group <- reference.group
  res$baseline.reference <- baseline.reference
  res$all.treatments <- all.treatments
  ##
  if (is.bin)
    class(res) <- c("summary.netmeta", "summary.netmetabin")
  else
    class(res) <- "summary.netmeta"
  
  res
}





#' @rdname summary.netmeta
#' @method print summary.netmeta
#' @export
#' @export print.summary.netmeta


print.summary.netmeta <- function(x,
                                  comb.fixed = x$comb.fixed,
                                  comb.random = x$comb.random,
                                  prediction = x$prediction,
                                  reference.group = x$reference.group,
                                  baseline.reference = x$baseline.reference,
                                  all.treatments = x$all.treatments,
                                  backtransf = x$backtransf,
                                  nchar.trts = x$nchar.trts,
                                  header = TRUE,
                                  digits = gs("digits"),
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
                                  ...) {
  
  
  meta:::chkclass(x, "summary.netmeta")
  ##
  is.bin <- inherits(x, "summary.netmetabin")
  ##
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  chkchar <- meta:::chkchar
  formatCI <- meta:::formatCI
  formatN <- meta:::formatN
  formatPT <- meta:::formatPT
  is.relative.effect <- meta:::is.relative.effect
  pasteCI <- meta:::pasteCI
  
  
  if (is.null(x$df.Q))
    oldversion <- TRUE
  else
    oldversion <- FALSE
  ##
  if (is.null(x$predict$lower))
    prediction <- FALSE
  ##
  if (is.null(x$backtransf))
    backtransf <- TRUE
  ##
  if (is.null(x$nchar.trts))
    nchar.trts <- 666
  
  
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
  chklogical(baseline.reference)
  ##
  chklogical(backtransf)
  chknumeric(nchar.trts, min = 1, length = 1)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.pval.Q, min = 1, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  chklogical(scientific.pval)
  ##
  chkchar(text.tau2)
  chkchar(text.tau)
  chkchar(text.I2)
  
  
  ##
  ## Additional arguments
  ##
  fun <- "print.summary.netmeta"
  addargs <- names(list(...))
  ##
  meta:::warnarg("logscale", addargs, fun, otherarg = "backtransf")
  
  
  k <- x$k
  m <- x$m
  n <- x$n
  sm <- x$sm
  
  sm.lab <- sm
  ##
  if (!backtransf & is.relative.effect(sm))
    sm.lab <- paste("log", sm, sep = "")

  ci.lab <- paste(round(100 * x$level.comb, 1), "%-CI", sep = "")

  
  TE.fixed    <- x$fixed$TE
  seTE.fixed  <- x$fixed$seTE
  lowTE.fixed <- x$fixed$lower
  uppTE.fixed <- x$fixed$upper
  ##
  TE.random    <- x$random$TE
  seTE.random  <- x$random$seTE
  lowTE.random <- x$random$lower
  uppTE.random <- x$random$upper
  ##
  lowTE.predict <- x$predict$lower
  uppTE.predict <- x$predict$upper
  ##
  if (!is.null(x$seq)) {
    TE.fixed <- TE.fixed[x$seq, x$seq]
    seTE.fixed <- seTE.fixed[x$seq, x$seq]
    lowTE.fixed <- lowTE.fixed[x$seq, x$seq]
    uppTE.fixed <- uppTE.fixed[x$seq, x$seq]
    ##
    if (!all(is.na(TE.random))) {
      TE.random <- TE.random[x$seq, x$seq]
      seTE.random <- seTE.random[x$seq, x$seq]
      lowTE.random <- lowTE.random[x$seq, x$seq]
      uppTE.random <- uppTE.random[x$seq, x$seq]
      ##
      lowTE.predict <- lowTE.predict[x$seq, x$seq]
      uppTE.predict <- uppTE.predict[x$seq, x$seq]
    }
  }
  
  
  noeffect <- 0
  ##
  if (backtransf & is.relative.effect(sm)) {
    noeffect <- 1
    ##
    TE.fixed    <- exp(TE.fixed)
    lowTE.fixed <- exp(lowTE.fixed)
    uppTE.fixed <- exp(uppTE.fixed)
    ##
    TE.random    <- exp(TE.random)
    lowTE.random <- exp(lowTE.random)
    uppTE.random <- exp(uppTE.random)
    ##
    if (prediction) {
      lowTE.predict <- exp(lowTE.predict)
      uppTE.predict <- exp(uppTE.predict)
    }
  }
  ##
  TE.fixed    <- round(TE.fixed, digits)
  lowTE.fixed <- round(lowTE.fixed, digits)
  uppTE.fixed <- round(uppTE.fixed, digits)
  ##
  TE.random    <- round(TE.random, digits)
  lowTE.random <- round(lowTE.random, digits)
  uppTE.random <- round(uppTE.random, digits)
  ##
  if (prediction) {
    lowTE.predict <- round(lowTE.predict, digits)
    uppTE.predict <- round(uppTE.predict, digits)
  }
  ##
  I2 <- round(100 * x$I2, digits.I2)
  lower.I2 <- round(100 * x$lower.I2, digits.I2)
  upper.I2 <- round(100 * x$upper.I2, digits.I2)
  
  
  if (header)
    matitle(x)
  
  
  if (reference.group != "" & missing(all.treatments))
    all.treatments <- FALSE
  ##
  if (reference.group != "")
    reference.group <- setref(reference.group, rownames(TE.fixed))
  
  
  if (comb.fixed | comb.random) {
    cat(paste("Number of studies: k = ", k, "\n", sep = ""))
    cat(paste("Number of treatments: n = ", n, "\n", sep = ""))
    cat(paste("Number of pairwise comparisons: m = ", m, "\n", sep = ""))
    if (!oldversion)
      cat(paste("Number of designs: d = ", x$d, "\n", sep = ""))
    
    
    if (reference.group != "")
      if (baseline.reference)
        comptext <- paste("comparison: ",
                          if (x$n == 2)
                            paste("'",
                                  treats(rownames(TE.fixed),
                                         nchar.trts)[rownames(TE.fixed)
                                                     != reference.group],
                                  "'", sep = "")
                          else
                            "other treatments",
                          " vs '",
                          treats(rownames(TE.fixed),
                                 nchar.trts)[rownames(TE.fixed)
                                             == reference.group],
                          "'", sep = "")
      else
        comptext <- paste("comparison: '",
                          treats(rownames(TE.fixed),
                                 nchar.trts)[rownames(TE.fixed)
                                             == reference.group],
                          "' vs ",
                          if (x$n == 2)
                            paste("'",
                                  treats(rownames(TE.fixed),
                                         nchar.trts)[rownames(TE.fixed)
                                                     != reference.group],
                                  "'", sep = "")
                          else
                            "other treatments", sep = "")
    
    
    if (comb.fixed) {
      if (all.treatments | reference.group != "") {
        text.fixed <- "Fixed effects model"
        ##
        if (x$method == "MH")
          text.fixed <- paste(text.fixed, "(Mantel-Haenszel method)")
        else if (x$method == "NCH")
          text.fixed <- paste(text.fixed, "(Non-central hypergeometric distribution)")
        ##
        cat(paste0("\n", text.fixed, "\n"))
      }
      ##
      if (all.treatments) {
        cat("\nTreatment estimate (sm = '", sm.lab, "'):\n", sep = "")
        ##
        TEf <- formatN(TE.fixed, digits = digits)
        rownames(TEf) <- treats(TEf, nchar.trts)
        colnames(TEf) <- treats(TEf, nchar.trts, FALSE)
        ##
        if (all(diag(TE.fixed) == noeffect))
          diag(TEf) <- "."
        ##
        prmatrix(TEf, quote = FALSE, right = TRUE)
        ##
        cat("\nLower ", 100 * x$level.comb, "%-confidence limit:\n", sep = "")
        ##
        lowTEf <- formatN(lowTE.fixed, digits = digits)
        rownames(lowTEf) <- treats(lowTEf, nchar.trts)
        colnames(lowTEf) <- treats(lowTEf, nchar.trts, FALSE)
        ##
        if (all(diag(lowTE.fixed) == noeffect))
          diag(lowTEf) <- "."
        ##
        prmatrix(lowTEf, quote = FALSE, right = TRUE)
        ##
        cat("\nUpper ", 100 * x$level.comb, "%-confidence limit:\n", sep = "")
        ##
        uppTEf <- formatN(uppTE.fixed, digits = digits)
        rownames(uppTEf) <- treats(uppTEf, nchar.trts)
        colnames(uppTEf) <- treats(uppTEf, nchar.trts, FALSE)
        ##
        if (all(diag(uppTE.fixed) == noeffect))
          diag(uppTEf) <- "."
        ##
        prmatrix(uppTEf, quote = FALSE, right = TRUE)
        ##
        ## Print prediction intervals
        ##
        if (!comb.random & prediction & x$df.Q >= 2) {
          cat("\nPrediction intervals\n")
          ##
          cat("\nLower ", 100 * x$level.predict, "%-prediction limit:\n", sep = "")
          ##
          lowTEp <- formatN(lowTE.predict, digits = digits)
          rownames(lowTEp) <- treats(lowTEp, nchar.trts)
          colnames(lowTEp) <- treats(lowTEp, nchar.trts, FALSE)
          ##
          if (all(diag(lowTE.predict) == noeffect))
            diag(lowTEp) <- "."
          ##
          prmatrix(lowTEp, quote = FALSE, right = TRUE)
          ##
          cat("\nUpper ", 100 * x$level.predict, "%-prediction limit:\n", sep = "")
          ##
          uppTEp <- formatN(uppTE.predict, digits = digits)
          rownames(uppTEp) <- treats(uppTEp, nchar.trts)
          colnames(uppTEp) <- treats(uppTEp, nchar.trts, FALSE)
          ##
          if (all(diag(uppTE.predict) == noeffect))
            diag(uppTEp) <- "."
          ##
          prmatrix(uppTEp, quote = FALSE, right = TRUE)
        }
      }
      if (reference.group != "") {
        if (all(colnames(TE.fixed) != reference.group))
          stop(paste("Argument 'reference.group' must match any of the following values: ",
                     paste(paste("'", colnames(TE.fixed), "'", sep = ""),
                           collapse = " - "), sep = ""))
        ##
        if (baseline.reference) {
          TE.fixed.b <- TE.fixed[, colnames(TE.fixed) == reference.group]
          lowTE.fixed.b <- lowTE.fixed[, colnames(lowTE.fixed) == reference.group]
          uppTE.fixed.b <- uppTE.fixed[, colnames(uppTE.fixed) == reference.group]
        }
        else {
          TE.fixed.b <- TE.fixed[rownames(TE.fixed) == reference.group, ]
          lowTE.fixed.b <- lowTE.fixed[rownames(lowTE.fixed) == reference.group, ]
          uppTE.fixed.b <- uppTE.fixed[rownames(uppTE.fixed) == reference.group, ]
        }
        ##
        ## Add prediction interval (or not)
        ##
        if (!comb.random & prediction & x$df.Q >= 2) {
          if (baseline.reference) {
            lowTE.predict.b <- lowTE.predict[, colnames(lowTE.predict) == reference.group]
            uppTE.predict.b <- uppTE.predict[, colnames(uppTE.predict) == reference.group]
          }
          else {
            lowTE.predict.b <- lowTE.predict[rownames(lowTE.predict) == reference.group, ]
            uppTE.predict.b <- uppTE.predict[rownames(uppTE.predict) == reference.group, ]
          }
          ##
          pi.lab <- paste(round(100 * x$level.predict, 1), "%-PI", sep = "")
          ##
          res <- cbind(formatN(TE.fixed.b, digits, text.NA = "NA",
                               big.mark = big.mark),
                       formatCI(formatN(round(lowTE.fixed.b, digits),
                                        digits, "NA", big.mark = big.mark),
                                formatN(round(uppTE.fixed.b, digits),
                                        digits, "NA", big.mark = big.mark)),
                       formatCI(formatN(round(lowTE.predict.b, digits),
                                        digits, "NA", big.mark = big.mark),
                                formatN(round(uppTE.predict.b, digits),
                                        digits, "NA", big.mark = big.mark)))
          dimnames(res) <- list(colnames(TE.fixed), c(sm.lab, ci.lab, pi.lab))
        }
        else {
          res <- cbind(formatN(TE.fixed.b, digits, text.NA = "NA",
                               big.mark = big.mark),
                       formatCI(formatN(round(lowTE.fixed.b, digits),
                                        digits, "NA", big.mark = big.mark),
                                formatN(round(uppTE.fixed.b, digits),
                                        digits, "NA", big.mark = big.mark))
                       )
          dimnames(res) <- list(colnames(TE.fixed), c(sm.lab, ci.lab))
        }
        ##
        if (TE.fixed.b[rownames(res) == reference.group] == noeffect)
          res[rownames(res) == reference.group, ] <- "."
        ##
        rownames(res) <- treats(rownames(res), nchar.trts)
        
        cat("\nTreatment estimate (sm = '", sm.lab,
            "', ", comptext, "):\n", sep = "")
        
        prmatrix(res, quote = FALSE, right = TRUE)
      }
    }
    
    
    if (comb.random) {
      if (all.treatments | reference.group != "")
        cat("\nRandom effects model\n")
      if (all.treatments) {
        cat("\nTreatment estimate (sm = '", sm.lab, "'):\n", sep = "")
        ##
        TEr <- formatN(TE.random, digits = digits)
        rownames(TEr) <- treats(TEr, nchar.trts)
        colnames(TEr) <- treats(TEr, nchar.trts, FALSE)
        ##
        if (all(diag(TE.random) == noeffect))
          diag(TEr) <- "."
        ##
        prmatrix(TEr, quote = FALSE, right = TRUE)
        ##
        cat("\nLower ", 100 * x$level.comb, "%-confidence limit:\n", sep = "")
        ##
        lowTEr <- formatN(lowTE.random, digits = digits)
        rownames(lowTEr) <- treats(lowTEr, nchar.trts)
        colnames(lowTEr) <- treats(lowTEr, nchar.trts, FALSE)
        ##
        if (all(diag(lowTE.random) == noeffect))
          diag(lowTEr) <- "."
        ##
        prmatrix(lowTEr, quote = FALSE, right = TRUE)
        ##
        cat("\nUpper ", 100 * x$level.comb, "%-confidence limit:\n", sep = "")
        ##
        uppTEr <- formatN(uppTE.random, digits = digits)
        rownames(uppTEr) <- treats(uppTEr, nchar.trts)
        colnames(uppTEr) <- treats(uppTEr, nchar.trts, FALSE)
        ##
        if (all(diag(uppTE.random) == noeffect))
          diag(uppTEr) <- "."
        ##
        prmatrix(uppTEr, quote = FALSE, right = TRUE)
        ##
        ## Print prediction intervals
        ##
        if (prediction & x$df.Q >= 2) {
          cat("\nPrediction intervals\n")
          ##
          cat("\nLower ", 100 * x$level.predict, "%-prediction limit:\n", sep = "")
          ##
          lowTEp <- formatN(lowTE.predict, digits = digits)
          rownames(lowTEp) <- treats(lowTEp, nchar.trts)
          colnames(lowTEp) <- treats(lowTEp, nchar.trts, FALSE)
          ##
          if (all(diag(lowTE.predict) == noeffect))
            diag(lowTEp) <- "."
          ##
          prmatrix(lowTEp, quote = FALSE, right = TRUE)
          ##
          cat("\nUpper ", 100 * x$level.predict, "%-prediction limit:\n", sep = "")
          ##
          uppTEp <- formatN(uppTE.predict, digits = digits)
          rownames(uppTEp) <- treats(uppTEp, nchar.trts)
          colnames(uppTEp) <- treats(uppTEp, nchar.trts, FALSE)
          ##
          if (all(diag(uppTE.predict) == noeffect))
            diag(uppTEp) <- "."
          ##
          prmatrix(uppTEp, quote = FALSE, right = TRUE)
        }
      }
      if (reference.group != "") {
        if (all(colnames(TE.random) != reference.group))
          stop(paste("Argument 'reference.group' must match any of the following values: ",
                     paste(paste("'", colnames(TE.random), "'", sep = ""),
                           collapse = " - "), sep = ""))
        ##
        if (baseline.reference) {
          TE.random.b <- TE.random[, colnames(TE.random) == reference.group]
          lowTE.random.b <- lowTE.random[, colnames(lowTE.random) == reference.group]
          uppTE.random.b <- uppTE.random[, colnames(uppTE.random) == reference.group]
        }
        else {
          TE.random.b <- TE.random[colnames(TE.random) == reference.group]
          lowTE.random.b <- lowTE.random[rownames(lowTE.random) == reference.group, ]
          uppTE.random.b <- uppTE.random[rownames(uppTE.random) == reference.group, ]
        }
        ##
        ## Add prediction interval (or not)
        ##
        if (prediction & x$df.Q >= 2) {
          if (baseline.reference) {
            lowTE.predict.b <- lowTE.predict[, colnames(lowTE.predict) == reference.group]
            uppTE.predict.b <- uppTE.predict[, colnames(uppTE.predict) == reference.group]
          }
          else {
            lowTE.predict.b <- lowTE.predict[rownames(lowTE.predict) == reference.group, ]
            uppTE.predict.b <- uppTE.predict[rownames(uppTE.predict) == reference.group, ]
          }
          ##
          pi.lab <- paste(round(100 * x$level.predict, 1), "%-PI", sep = "")
          ##
          res <- cbind(formatN(TE.random.b, digits, text.NA = "NA",
                               big.mark = big.mark),
                       formatCI(formatN(round(lowTE.random.b, digits),
                                        digits, "NA", big.mark = big.mark),
                                formatN(round(uppTE.random.b, digits),
                                        digits, "NA", big.mark = big.mark)),
                       formatCI(formatN(round(lowTE.predict.b, digits),
                                        digits, "NA", big.mark = big.mark),
                                formatN(round(uppTE.predict.b, digits),
                                        digits, "NA", big.mark = big.mark))
                       )
          dimnames(res) <- list(colnames(TE.fixed), c(sm.lab, ci.lab, pi.lab))
        }
        else {
          res <- cbind(formatN(TE.random.b, digits, text.NA = "NA",
                               big.mark = big.mark),
                       formatCI(formatN(round(lowTE.random.b, digits),
                                        digits, "NA", big.mark = big.mark),
                                formatN(round(uppTE.random.b, digits),
                                        digits, "NA", big.mark = big.mark)))
          dimnames(res) <- list(colnames(TE.fixed), c(sm.lab, ci.lab))
        }
        ##
        if (!is.na(TE.random.b[rownames(res) == reference.group]) &&
            TE.random.b[rownames(res) == reference.group] == noeffect)
          res[rownames(res) == reference.group, ] <- "."
        ##
        rownames(res) <- treats(rownames(res), nchar.trts)
        
        cat("\nTreatment estimate (sm = '", sm.lab,
            "', ", comptext, "):\n", sep = "")
        
        prmatrix(res, quote = FALSE, right = TRUE)
      }
    }
    ##
    zlab <- "z"
    
    
    if (!is.null(x$tau.preset))
      tau <- x$tau.preset
    else
      tau <- x$tau
    ##
    if (is.bin)
      hi.txt <- "inconsistency (between designs)"
    else if (x$d == 1)
      hi.txt <- paste0("heterogeneity")
    else
      hi.txt <- "heterogeneity / inconsistency"
    ##
    if (!is.bin)
      cat(paste0("\nQuantifying ", hi.txt, ":\n",
                 formatPT(tau^2,
                          lab = TRUE, labval = text.tau2,
                          digits = digits.tau2,
                          lab.NA = "NA", big.mark = big.mark),
                 "; ",
                 formatPT(tau,
                          lab = TRUE, labval = text.tau,
                          digits = digits.tau,
                          lab.NA = "NA", big.mark = big.mark),
                 if (!is.na(I2))
                   paste0("; ", text.I2, " = ", round(I2, digits.I2), "%"),
                 if (!(is.na(lower.I2) | is.na(upper.I2)))
                   pasteCI(lower.I2, upper.I2, digits.I2, big.mark, unit = "%"),
                 "\n")
          )
    
    
    if (m > 1) {

      if (is.bin) {
        Q.overall <- x$Q.inconsistency
        df.Q.overall <- x$df.Q.inconsistency
        pval.Q.overall <- formatPT(x$pval.Q.inconsistency,
                                   digits = digits.pval.Q,
                                   scientific = scientific.pval)
      }
      else {
        Q.overall <- x$Q
        if (oldversion) {
          df.Q.overall <- x$df
          pval.Q.overall <- ifelse(df.Q.overall == 0, "--",
                                   formatPT(x$pval.Q,
                                            digits = digits.pval.Q,
                                            scientific = scientific.pval))
        }
        else {
          df.Q.overall <- x$df.Q
          pval.Q.overall <- formatPT(x$pval.Q,
                                     digits = digits.pval.Q,
                                     scientific = scientific.pval)
        }
      }
      
      if (is.bin & x$d == 1)
        cat("")
      else if (x$d == 1 | is.bin |
               is.na(x$Q.heterogeneity) | is.na(x$Q.inconsistency)) {
        Qdata <- cbind(round(Q.overall, digits.Q), df.Q.overall,
                       pval.Q.overall)
        
        dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
        ##
        cat(paste0("\nTest of ", hi.txt, ":\n"))
        prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
      }
      else {
        Qs <- c(x$Q, x$Q.heterogeneity, x$Q.inconsistency)
        df.Qs <- c(x$df.Q, x$df.Q.heterogeneity, x$df.Q.inconsistency)
        pval.Qs <- c(x$pval.Q, x$pval.Q.heterogeneity, x$pval.Q.inconsistency)
        pval.Qs <- formatPT(pval.Qs, digits = digits.pval.Q,
                            scientific = scientific.pval)
        cat(paste0("\nTests of heterogeneity (within designs) and ",
                   "inconsistency (between designs):\n"))
        Qdata <- data.frame(Q = round(Qs, digits.Q),
                            df = df.Qs,
                            pval = pval.Qs)
        names(Qdata) <- c("Q", "d.f.", "p-value")
        rownames(Qdata) <- c("Total",
                             "Within designs",
                             "Between designs")
        prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
      }
    }
    
    
    if (!is.null(x$tau.preset)) {
      cat("\nDetails:")
      ##
      tau2 <- x$tau.preset^2
      tau2 <- formatPT(tau2, lab = TRUE, labval = text.tau2,
                       digits = digits.tau2,
                       lab.NA = "NA", big.mark = big.mark)
      ##
      cat(paste("\n- Preset between-study variance: ",
                tau2, "\n", sep = ""))
    }
  }
  
  
  if (any(rownames(TE.fixed) != treats(TE.fixed, nchar.trts))) {
    abbr <- unique(treats(TE.fixed, nchar.trts))
    full <- unique(rownames(TE.fixed))
    ##
    tmat <- data.frame(abbr, full)
    names(tmat) <- c("Abbreviation", "Treatment name")
    tmat <- tmat[order(tmat$Abbreviation), ]
    ##
    cat("\nLegend:\n")
    prmatrix(tmat, quote = FALSE, right = TRUE,
             rowlab = rep("", length(abbr))) 
  }
  
  
  invisible(NULL)
}
