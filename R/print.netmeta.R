#' Print method for objects of class netmeta
#' 
#' @description
#' Print method for objects of class \code{netmeta}.
#' 
#' @param x An object of class \code{netmeta}.
#' @param common A logical indicating whether results for the common
#'   effects model should be printed.
#' @param random A logical indicating whether results for the random
#'   effects model should be printed.
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
#' @param overall.hetstat A logical indicating whether to print heterogeneity
#'   measures.
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
#' @param digits.stat Minimal number of significant digits for tests
#'   of overall effect, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall effects, see \code{print.default}.
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
#' @param details.methods A logical specifying whether details on statistical
#'   methods should be printed.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments.
#' 
#' @rdname netmeta
#' @method print netmeta
#' @export


print.netmeta <- function(x,
                          common = x$common,
                          random = x$random,
                          prediction = x$prediction,
                          reference.group = x$reference.group,
                          baseline.reference = x$baseline.reference,
                          all.treatments = x$all.treatments,
                          overall.hetstat = x$overall.hetstat,
                          backtransf = x$backtransf,
                          nchar.trts = x$nchar.trts,
                          header = TRUE,
                          digits = gs("digits"),
                          digits.stat = gs("digits.stat"),
                          digits.pval = max(gs("digits.pval"), 2),
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
                          #
                          details.methods = gs("details.netmeta"),
                          legend = gs("legend.netmeta"),
                          #
                          warn.deprecated = gs("warn.deprecated"),
                          ##
                          ...) {
  
  
  ##
  ##
  ## (1) Check for netmeta object and upgrade object
  ##
  ##
  chkclass(x, "netmeta")
  x <- updateversion(x)
  ##
  is.bin <- inherits(x, "netmetabin")
  ##  
  if (is.null(x$df.Q))
    oldversion <- TRUE
  else
    oldversion <- FALSE
  ##
  if (is.null(x$lower.predict))
    prediction <- FALSE
  ##
  if (is.null(x$nchar.trts))
    nchar.trts <- 666
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chklogical(prediction)
  chklogical(baseline.reference)
  ##
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
  chklogical(scientific.pval)
  ##
  chkchar(text.tau2)
  chkchar(text.tau)
  chkchar(text.I2)
  ##
  chklogical(details.methods)
  chklogical(legend)
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
  #
  chklogical(overall.hetstat)
  #
  backtransf <-
    deprecated(backtransf, missing(backtransf), args, "logscale")
  if (is.untransformed(x$sm))
    backtransf <- TRUE
  backtransf <- replaceNULL(backtransf, TRUE)
  chklogical(backtransf)
  
  
  ##
  ##
  ## (3) Some additional settings
  ##
  ##
  k <- x$k
  m <- x$m
  n <- x$n
  sm <- sm.lab <- x$sm
  #
  if (!backtransf & (is.relative.effect(sm) | sm == "VE"))
    sm.lab <- paste0("log", if (sm == "VE") "VR" else sm)
  ##
  ci.lab <- paste0(round(100 * x$level.ma, 1), "%-CI")
  
  
  ##
  ##
  ## (4) Set and backtransform results of meta-analysis
  ##
  ##
  TE.common <- x$TE.common
  seTE.common <- x$seTE.common
  lowTE.common <- x$lower.common
  uppTE.common <- x$upper.common
  statistic.common <- replaceNULL(x$statistic.common, x$zval.common)
  pval.common <- x$pval.common
  ##
  TE.random <- x$TE.random
  seTE.random <- x$seTE.random
  lowTE.random <- x$lower.random
  uppTE.random <- x$upper.random
  statistic.random <- replaceNULL(x$statistic.random, x$zval.random)
  pval.random <- x$pval.random
  ##
  lower.predict <- x$lower.predict
  upper.predict <- x$upper.predict
  ##
  if (!is.null(x$seq)) {
    TE.common <- TE.common[x$seq, x$seq]
    seTE.common <- seTE.common[x$seq, x$seq]
    lowTE.common <- lowTE.common[x$seq, x$seq]
    uppTE.common <- uppTE.common[x$seq, x$seq]
    statistic.common <- statistic.common[x$seq, x$seq]
    pval.common <- pval.common[x$seq, x$seq]
    ##
    if (!all(is.na(TE.random))) {
      TE.random <- TE.random[x$seq, x$seq]
      seTE.random <- seTE.random[x$seq, x$seq]
      lowTE.random <- lowTE.random[x$seq, x$seq]
      uppTE.random <- uppTE.random[x$seq, x$seq]
      statistic.random <- statistic.random[x$seq, x$seq]
      pval.random <- pval.random[x$seq, x$seq]
      ##
      lower.predict <- lower.predict[x$seq, x$seq]
      upper.predict <- upper.predict[x$seq, x$seq]
    }
  }
  ##  
  noeffect <- 1L * (backtransf & is.relative.effect(sm))
  #
  if (backtransf) {
    TE.common    <- backtransf(TE.common, sm)
    lowTE.common <- backtransf(lowTE.common, sm)
    uppTE.common <- backtransf(uppTE.common, sm)
    #
    TE.random    <- backtransf(TE.random, sm)
    lowTE.random <- backtransf(lowTE.random, sm)
    uppTE.random <- backtransf(uppTE.random, sm)
    #
    if (prediction) {
      lower.predict <- backtransf(lower.predict, sm)
      upper.predict <- backtransf(upper.predict, sm)
    }
    #
    # Switch lower and upper limit for VE if results have been
    # backtransformed
    #
    if (sm == "VE") {
      tmp.l <- lowTE.common
      lowTE.common <- uppTE.common
      uppTE.common <- tmp.l
      #
      tmp.l <- lowTE.random
      lowTE.random <- uppTE.random
      uppTE.random <- tmp.l
      #
      if (prediction) {
        tmp.l <- lower.predict
        lower.predict <- upper.predict
        upper.predict <- tmp.l
      }
    }
  }
  ##
  TE.common <- round(TE.common, digits)
  lowTE.common <- round(lowTE.common, digits)
  uppTE.common <- round(uppTE.common, digits)
  statistic.common <- round(statistic.common, digits.stat)
  ##
  TE.random    <- round(TE.random, digits)
  lowTE.random <- round(lowTE.random, digits)
  uppTE.random <- round(uppTE.random, digits)
  statistic.random <- round(statistic.random, digits.stat)
  ##
  if (prediction) {
    lower.predict <- round(lower.predict, digits)
    upper.predict <- round(upper.predict, digits)
  }
  ##
  I2 <- round(100 * x$I2, digits.I2)
  lower.I2 <- round(100 * x$lower.I2, digits.I2)
  upper.I2 <- round(100 * x$upper.I2, digits.I2)
  
  
  ##
  ##
  ## (5) Print result for network meta-analysis
  ##
  ##
  if (header)
    matitle(x)
  ##  
  if (reference.group != "" & is.null(all.treatments))
    all.treatments <- FALSE
  ##
  if (reference.group != "")
    reference.group <- setref(reference.group, rownames(TE.common))
  ##  
  if (common | random) {
    cat("Number of studies: k = ", k, "\n", sep = "")
    cat("Number of pairwise comparisons: m = ", m, "\n", sep = "")
    if (!is.null(x$n.trts))
      cat("Number of observations: o = ",
          round(sum(x$n.trts, na.rm = TRUE), 1),
          "\n",
          sep = "")
    cat("Number of treatments: n = ", n, "\n", sep = "")
    if (!oldversion)
      cat("Number of designs: d = ", x$d, "\n", sep = "")
    ##
    if (reference.group != "") {
      if (baseline.reference)
        comptext <-
          paste0("comparison: ",
                 if (x$n == 2)
                   paste0("'",
                          treats(rownames(TE.common),
                                 nchar.trts)[rownames(TE.common)
                                             != reference.group],
                          "'")
                 else
                   "other treatments",
                 " vs '",
                 treats(rownames(TE.common),
                        nchar.trts)[rownames(TE.common)
                                    == reference.group],
                 "'")
      else
        comptext <-
          paste0("comparison: '",
                 treats(rownames(TE.common),
                        nchar.trts)[rownames(TE.common)
                                    == reference.group],
                 "' vs ",
                 if (x$n == 2)
                   paste0("'",
                          treats(rownames(TE.common),
                                 nchar.trts)[rownames(TE.common)
                                             != reference.group],
                          "'")
                 else
                   "other treatments")
    }
    ##
    if (common) {
      if (all.treatments | reference.group != "") {
        text.common <- "Common effects model"
        ##
        if (x$method == "MH")
          text.common <-
            paste(text.common, "(Mantel-Haenszel method)")
        else if (x$method == "NCH")
          text.common <-
            paste(text.common, "(Non-central hypergeometric distribution)")
        ##
        cat("\n", text.common, "\n", sep = "")
      }
      ##
      if (all.treatments) {
        cat("\nTreatment estimate (sm = '", sm.lab, "'):\n", sep = "")
        ##
        TEf <- formatN(TE.common, digits = digits)
        rownames(TEf) <- treats(TEf, nchar.trts)
        colnames(TEf) <- treats(TEf, nchar.trts, FALSE)
        ##
        if (all(diag(TE.common) == noeffect))
          diag(TEf) <- "."
        ##
        prmatrix(TEf, quote = FALSE, right = TRUE)
        ##
        cat("\nLower ", 100 * x$level.ma, "%-confidence limit:\n", sep = "")
        ##
        lowTEf <- formatN(lowTE.common, digits = digits)
        rownames(lowTEf) <- treats(lowTEf, nchar.trts)
        colnames(lowTEf) <- treats(lowTEf, nchar.trts, FALSE)
        ##
        if (all(diag(lowTE.common) == noeffect))
          diag(lowTEf) <- "."
        ##
        prmatrix(lowTEf, quote = FALSE, right = TRUE)
        ##
        cat("\nUpper ", 100 * x$level.ma, "%-confidence limit:\n", sep = "")
        ##
        uppTEf <- formatN(uppTE.common, digits = digits)
        rownames(uppTEf) <- treats(uppTEf, nchar.trts)
        colnames(uppTEf) <- treats(uppTEf, nchar.trts, FALSE)
        ##
        if (all(diag(uppTE.common) == noeffect))
          diag(uppTEf) <- "."
        ##
        prmatrix(uppTEf, quote = FALSE, right = TRUE)
        ##
        ## Print prediction intervals
        ##
        if (!random && prediction & x$df.Q >= 2) {
          cat("\nPrediction intervals\n")
          ##
          cat("\nLower ", 100 * x$level.predict, "%-prediction limit:\n",
              sep = "")
          ##
          lowTEp <- formatN(lower.predict, digits = digits)
          rownames(lowTEp) <- treats(lowTEp, nchar.trts)
          colnames(lowTEp) <- treats(lowTEp, nchar.trts, FALSE)
          ##
          if (all(diag(lower.predict) == noeffect))
            diag(lowTEp) <- "."
          ##
          prmatrix(lowTEp, quote = FALSE, right = TRUE)
          ##
          cat("\nUpper ", 100 * x$level.predict, "%-prediction limit:\n",
              sep = "")
          ##
          uppTEp <- formatN(upper.predict, digits = digits)
          rownames(uppTEp) <- treats(uppTEp, nchar.trts)
          colnames(uppTEp) <- treats(uppTEp, nchar.trts, FALSE)
          ##
          if (all(diag(upper.predict) == noeffect))
            diag(uppTEp) <- "."
          ##
          prmatrix(uppTEp, quote = FALSE, right = TRUE)
        }
      }
      if (reference.group != "") {
        if (all(colnames(TE.common) != reference.group))
          stop("Argument 'reference.group' must match any of ",
               "the following values: ",
               paste(paste0("'", colnames(TE.common), "'"),
                     collapse = " - "))
        ##
        if (baseline.reference) {
          TE.common.b <- TE.common[, colnames(TE.common) == reference.group]
          lowTE.common.b <-
            lowTE.common[, colnames(lowTE.common) == reference.group]
          uppTE.common.b <-
            uppTE.common[, colnames(uppTE.common) == reference.group]
          statistic.common.b <-
            statistic.common[, colnames(statistic.common) == reference.group]
          pval.common.b <-
            pval.common[, colnames(pval.common) == reference.group]
        }
        else {
          TE.common.b <- TE.common[rownames(TE.common) == reference.group, ]
          lowTE.common.b <-
            lowTE.common[rownames(lowTE.common) == reference.group, ]
          uppTE.common.b <-
            uppTE.common[rownames(uppTE.common) == reference.group, ]
          statistic.common.b <-
            statistic.common[rownames(statistic.common) == reference.group, ]
          pval.common.b <-
            pval.common[rownames(pval.common) == reference.group, ]
        }
        ##
        res <- cbind(formatN(TE.common.b, digits, text.NA = "NA",
                             big.mark = big.mark),
                     formatCI(formatN(round(lowTE.common.b, digits),
                                      digits, "NA", big.mark = big.mark),
                              formatN(round(uppTE.common.b, digits),
                                      digits, "NA", big.mark = big.mark)),
                     formatN(statistic.common.b, digits.stat, text.NA = "NA",
                             big.mark = big.mark),
                     formatPT(pval.common.b,
                              digits = digits.pval,
                              scientific = scientific.pval)
                     )
        dimnames(res) <-
          list(colnames(TE.common), c(sm.lab, ci.lab, "z", "p-value"))
        ##
        ## Add prediction interval (or not)
        ##
        if (!random && prediction & x$df.Q >= 2) {
          if (baseline.reference) {
            lower.predict.b <-
              lower.predict[, colnames(lower.predict) == reference.group]
            upper.predict.b <-
              upper.predict[, colnames(upper.predict) == reference.group]
          }
          else {
            lower.predict.b <-
              lower.predict[rownames(lower.predict) == reference.group, ]
            upper.predict.b <-
              upper.predict[rownames(upper.predict) == reference.group, ]
          }
          ##
          pi.lab <- paste0(round(100 * x$level.predict, 1), "%-PI")
          ##
          res <- cbind(res,
                       rep_len("", nrow(res)),
                       formatCI(formatN(round(lower.predict.b, digits),
                                        digits, "NA", big.mark = big.mark),
                                formatN(round(upper.predict.b, digits),
                                        digits, "NA", big.mark = big.mark)))
          colnames(res)[ncol(res) - 1] <- ".not.a.column."
          colnames(res)[ncol(res)] <- pi.lab
        }
        ##
        if (TE.common.b[rownames(res) == reference.group] == noeffect) {
          res[rownames(res) == reference.group, ] <- "."
          res[, colnames(res) == ".not.a.column."] <- ""
          colnames(res)[colnames(res) == ".not.a.column."] <- ""
        }
        ##
        rownames(res) <- treats(rownames(res), nchar.trts)
        ##
        cat("\nTreatment estimate (sm = '", sm.lab,
            "', ", comptext, "):\n", sep = "")
        ##
        prmatrix(res, quote = FALSE, right = TRUE)
      }
    }
    ##
    if (random) {
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
        cat("\nLower ", 100 * x$level.ma, "%-confidence limit:\n", sep = "")
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
        cat("\nUpper ", 100 * x$level.ma, "%-confidence limit:\n", sep = "")
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
          cat("\nLower ", 100 * x$level.predict, "%-prediction limit:\n",
              sep = "")
          ##
          lowTEp <- formatN(lower.predict, digits = digits)
          rownames(lowTEp) <- treats(lowTEp, nchar.trts)
          colnames(lowTEp) <- treats(lowTEp, nchar.trts, FALSE)
          ##
          if (all(diag(lower.predict) == noeffect))
            diag(lowTEp) <- "."
          ##
          prmatrix(lowTEp, quote = FALSE, right = TRUE)
          ##
          cat("\nUpper ", 100 * x$level.predict, "%-prediction limit:\n",
              sep = "")
          ##
          uppTEp <- formatN(upper.predict, digits = digits)
          rownames(uppTEp) <- treats(uppTEp, nchar.trts)
          colnames(uppTEp) <- treats(uppTEp, nchar.trts, FALSE)
          ##
          if (all(diag(upper.predict) == noeffect))
            diag(uppTEp) <- "."
          ##
          prmatrix(uppTEp, quote = FALSE, right = TRUE)
        }
      }
      if (reference.group != "") {
        if (all(colnames(TE.random) != reference.group))
          stop("Argument 'reference.group' must match any of ",
               "the following values: ",
               paste(paste0("'", colnames(TE.random), "'"),
                     collapse = " - "))
        ##
        if (baseline.reference) {
          TE.random.b <-
            TE.random[, colnames(TE.random) == reference.group]
          lowTE.random.b <-
            lowTE.random[, colnames(lowTE.random) == reference.group]
          uppTE.random.b <-
            uppTE.random[, colnames(uppTE.random) == reference.group]
          statistic.random.b <-
            statistic.random[, colnames(statistic.random) == reference.group]
          pval.random.b <-
            pval.random[, colnames(pval.random) == reference.group]
        }
        else {
          TE.random.b <-
            TE.random[colnames(TE.random) == reference.group]
          lowTE.random.b <-
            lowTE.random[rownames(lowTE.random) == reference.group, ]
          uppTE.random.b <-
            uppTE.random[rownames(uppTE.random) == reference.group, ]
          statistic.random.b <-
            statistic.random[rownames(statistic.random) == reference.group, ]
          pval.random.b <-
            pval.random[rownames(pval.random) == reference.group, ]
        }
        ##
        res <- cbind(formatN(TE.random.b, digits, text.NA = "NA",
                             big.mark = big.mark),
                     formatCI(formatN(round(lowTE.random.b, digits),
                                      digits, "NA", big.mark = big.mark),
                              formatN(round(uppTE.random.b, digits),
                                      digits, "NA", big.mark = big.mark)),
                     formatN(statistic.random.b, digits.stat, text.NA = "NA",
                             big.mark = big.mark),
                     formatPT(pval.random.b,
                              digits = digits.pval,
                              scientific = scientific.pval)
                     )
        dimnames(res) <-
          list(colnames(TE.random), c(sm.lab, ci.lab, "z", "p-value"))
        ##
        ## Add prediction interval (or not)
        ##
        if (prediction & x$df.Q >= 2) {
          if (baseline.reference) {
            lower.predict.b <-
              lower.predict[, colnames(lower.predict) == reference.group]
            upper.predict.b <-
              upper.predict[, colnames(upper.predict) == reference.group]
          }
          else {
            lower.predict.b <-
              lower.predict[rownames(lower.predict) == reference.group, ]
            upper.predict.b <-
              upper.predict[rownames(upper.predict) == reference.group, ]
          }
          ##
          pi.lab <- paste0(round(100 * x$level.predict, 1), "%-PI")
          ##
          res <- cbind(res,
                       rep_len("", nrow(res)),
                       formatCI(formatN(round(lower.predict.b, digits),
                                        digits, "NA", big.mark = big.mark),
                                formatN(round(upper.predict.b, digits),
                                        digits, "NA", big.mark = big.mark))
                       )
          colnames(res)[ncol(res) - 1] <- ".not.a.column."
          colnames(res)[ncol(res)] <- pi.lab
        }
        ##
        if (!is.na(TE.random.b[rownames(res) == reference.group]) &&
            TE.random.b[rownames(res) == reference.group] == noeffect) {
          res[rownames(res) == reference.group, ] <- "."
          res[, colnames(res) == ".not.a.column."] <- ""
          colnames(res)[colnames(res) == ".not.a.column."] <- ""
        }
        ##
        rownames(res) <- treats(rownames(res), nchar.trts)
        ##
        cat("\nTreatment estimate (sm = '", sm.lab,
            "', ", comptext, "):\n", sep = "")
        ##
        prmatrix(res, quote = FALSE, right = TRUE)
      }
    }
    ##
    zlab <- "z"
    ##
    if (!is.null(x$tau.preset))
      tau <- x$tau.preset
    else
      tau <- x$tau
    ##
    if (overall.hetstat) {
      if (is.bin)
        hi.txt <- "inconsistency (between designs)"
      else if (x$d == 1)
        hi.txt <- "heterogeneity"
      else
        hi.txt <- "heterogeneity / inconsistency"
      ##
      if (!is.bin)
        cat("\nQuantifying ", hi.txt, ":\n",
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
              pasteCI(lower.I2, upper.I2,
                      digits.I2, big.mark, unit = "%"),
            "\n",
            sep = ""
        )
      ##
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
        ##
        if (is.bin & x$d == 1)
          cat("")
        else if (x$d == 1 | is.bin |
                 is.na(x$Q.heterogeneity) | is.na(x$Q.inconsistency)) {
          Qdata <- cbind(round(Q.overall, digits.Q), df.Q.overall,
                         pval.Q.overall)
          ##
          dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
          ##
          cat("\nTest of ", hi.txt, ":\n", sep = "")
          prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
        }
        else {
          Qs <- c(x$Q, x$Q.heterogeneity, x$Q.inconsistency)
          df.Qs <- c(x$df.Q, x$df.Q.heterogeneity, x$df.Q.inconsistency)
          pval.Qs <- c(x$pval.Q, x$pval.Q.heterogeneity, x$pval.Q.inconsistency)
          pval.Qs <- formatPT(pval.Qs, digits = digits.pval.Q,
                              scientific = scientific.pval)
          cat("\nTests of heterogeneity (within designs) and inconsistency",
              if (options()$width < 77) "\n" else " ",
              "(between designs):\n",
              sep = "")
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
    }
  }
  #
  # Print details of network meta-analysis methods
  #
  if (details.methods) {
    text.details <- catmeth(x, random, text.tau2, digits.tau2, big.mark)
    #
    cat(text.details)
  }
  #
  # Add legend with abbreviated treatment labels
  #
  legendabbr(unique(rownames(TE.common)),
             treats(TE.common, nchar.trts),
             legend)
  
  invisible(NULL)
}
