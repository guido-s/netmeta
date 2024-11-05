#' Print detailed results of network meta-analysis
#' 
#' @description
#' Print method for objects of class \code{summary.netmeta}.
#' 
#' @param x An object of class \code{summary.netmeta}.
#' @param sortvar An optional vector used to sort individual studies
#'   (must be of same length as \code{x$TE}).
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
#' @param details A logical indicating whether further details for
#'   individual studies should be printed.
#' @param nma A logical indicating whether summary results of network
#'   meta-analysis should be printed.
#' @param overall.hetstat A logical indicating whether to print heterogeneity
#'   measures.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param nchar.studlab A numeric defining the minimum number of
#'   characters used to create unique study labels.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard
#'   deviations and standard errors, see \code{print.default}.
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
#' @param truncate An optional vector used to truncate the printout of
#'   results for individual studies (must be a logical vector of
#'   length corresponding to the number of pairwise comparisons
#'   \code{x$TE} or contain numerical values).
#' @param text.truncate A character string printed if study results
#'   were truncated from the printout.
#' @param details.methods A logical specifying whether details on statistical
#'   methods should be printed.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{netmeta}}, \code{\link{summary.netmeta}}
#' 
#' @keywords print
#' 
#' @examples
#' data(smokingcessation)
#' 
#' # Transform data from arm-based format to contrast-based format
#' #
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' 
#' # Conduct random effects network meta-analysis and print detailed
#' # summary
#' #
#' net1 <- netmeta(p1, common = FALSE)
#' summary(net1)
#' 
#' \dontrun{
#' data(Senn2013)
#' 
#' # Conduct common effects network meta-analysis
#' #
#' net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD", random = FALSE, ref = "plac")
#' snet2 <- summary(net2)
#' print(snet2, digits = 3)
#' 
#' # Only show individual study results for multi-arm studies
#' #
#' print(snet2, digits = 3, truncate = multiarm)
#' 
#' # Only show first three individual study results
#' #
#' print(snet2, digits = 3, truncate = 1:3)
#' 
#' # Only show individual study results for Kim2007 and Willms1999
#' #
#' print(snet2, digits = 3, truncate = c("Kim2007", "Willms1999"))
#' 
#' # Only show individual study results for studies starting with the
#' # letter "W"
#' #
#' print(snet2, ref = "plac", digits = 3,
#'   truncate = substring(studlab, 1, 1) == "W")
#' 
#' # Conduct random effects network meta-analysis
#' #
#' net3 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD", common = FALSE, ref = "plac")
#' print(summary(net3), digits = 3)
#' }
#' 
#' @method print summary.netmeta
#' @export


print.summary.netmeta <- function(x,
                                  sortvar,
                                  common = x$x$common,
                                  random = x$x$random,
                                  prediction = x$prediction,
                                  reference.group = x$reference.group,
                                  baseline.reference = x$baseline.reference,
                                  all.treatments = x$all.treatments,
                                  details = TRUE, nma = TRUE,
                                  ##
                                  overall.hetstat = x$overall.hetstat,
                                  backtransf = x$backtransf,
                                  nchar.trts = x$nchar.trts,
                                  nchar.studlab = x$nchar.studlab,
                                  digits = gs("digits"),
                                  digits.se = gs("digits.se"),
                                  digits.pval.Q = max(gs("digits.pval.Q"), 2),
                                  digits.Q = gs("digits.Q"),
                                  digits.tau2 = gs("digits.tau2"),
                                  digits.I2 = gs("digits.I2"),
                                  scientific.pval = gs("scientific.pval"),
                                  big.mark = gs("big.mark"),
                                  truncate,
                                  text.truncate = "*** Output truncated ***",
                                  #
                                  details.methods = gs("details.netmeta"),
                                  legend = gs("legend.netmeta"),
                                  #
                                  warn.deprecated = gs("warn.deprecated"),
                                  ##
                                  ...) {
  
  
  ##
  ##
  ## (1) Check for summary.netmeta object and upgrade object
  ##
  ##
  chkclass(x, "summary.netmeta")
  x <- updateversion(x)
  ##
  k.all <- length(x$x$TE)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chklogical(prediction)
  chklogical(baseline.reference)
  ##
  chknumeric(nchar.trts, min = 1, length = 1)
  if (is.null(nchar.studlab))
    nchar.studlab <- 666
  chknumeric(nchar.studlab, min = 1, length = 1)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.pval.Q, min = 1, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  ##
  chklogical(scientific.pval)
  ##
  chklogical(nma)
  chklogical(details.methods)
  chklogical(legend)
  ##
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  ## Catch 'truncate' from network meta-analysis object:
  ##
  missing.truncate <- missing(truncate)
  if (!missing.truncate) {
    truncate <- catch("truncate", mc, x$x, sfsp)
    ##
    if (is.null(truncate))
      truncate <- catch("truncate", mc, x$x$data, sfsp)
    ##
    if (length(truncate) > k.all)
      stop("Length of argument 'truncate' is too long.",
           call. = FALSE)
    else if (length(truncate) < k.all) {
      if (is.numeric(truncate)) {
        if (any(is.na(truncate)) | max(truncate) > k.all | min(truncate) < 0)
          stop("Numeric values in argument 'truncate' must be between 1 and ",
               k.all, ".",
               call. = FALSE)
        truncate2 <- rep(FALSE, k.all)
        truncate2[truncate] <- TRUE
        truncate <- truncate2
      }
      else if (is.character(truncate)) {
        if (any(!(truncate %in% x$x$studlab)))
          stop("At least one value of argument 'truncate' does not ",
               "match a study label.",
               call. = FALSE)
        truncate2 <- rep(FALSE, k.all)
        truncate2[x$x$studlab %in% truncate] <- TRUE
        truncate <- truncate2
      }
      else
        stop("Argument 'truncate' must contain integers or study labels if ",
             "length differs from number of treatment effects.",
             call. = FALSE)
    }
  }
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
  chklogical(backtransf)
  
  
  ##
  ##
  ## (3) Some additional settings
  ##
  ##
  if (!inherits(x, "summary.netmetabin")) {
    ##
    if (missing(sortvar)) sortvar <- 1:k.all
    ##
    if (length(sortvar) != k.all)
      stop("'x' and 'sortvar' have different length")
    ##
    ci.lab <- paste0(round(100 * x$level, 1), "%-CI")
    ##
    sm <- sm.lab <- x$sm
    #
    if (!backtransf & (is.relative.effect(sm) | sm == "VE"))
      sm.lab <- paste0("log", if (sm == "VE") "VR" else sm)
    ##    
    trts <- x$x$trts
    ##
    treat1 <-
      as.character(factor(x$x$treat1, levels = trts,
                          labels = treats(trts, nchar.trts)))
    treat2 <-
      as.character(factor(x$x$treat2, levels = trts,
                          labels = treats(trts, nchar.trts)))
    
    
    ##
    ##
    ## (4) Print title and details
    ##
    ##
    matitle(x)
    ##  
    if (details) {
      multiarm <- any(x$x$narms > 2)
      cat("Original data",
          ifelse(multiarm & (common | random),
                 paste(" (with adjusted standard errors for",
                       "multi-arm studies)"),
                 ""),
          ":\n\n",
          sep = ""
      )
      ##
      res <- data.frame(treat1,
                        treat2,
                        TE = formatN(x$x$TE, digits, text.NA = "NA",
                                     big.mark = big.mark),
                        seTE = formatN(x$x$seTE, digits.se, text.NA = "NA",
                                       big.mark = big.mark))
      ##
      if (multiarm) {
        if (is.null(x$x$seTE.adj.common))
          seTE.adj <- x$x$seTE.adj
        else
          seTE.adj <- x$x$seTE.adj.common
        ##
        if (common & random & !is.null(x$x$seTE.adj.random)) {
          res$seTE.adj.c <- format(round(seTE.adj, digits.se))
          res$seTE.adj.r <- format(round(x$x$seTE.adj.random, digits.se))
        }
        else if (common)
          res$seTE.adj <- format(round(seTE.adj, digits.se))
        else if (random & !is.null(x$x$seTE.adj.random))
          res$seTE.adj <- format(round(x$x$seTE.adj.random, digits.se))
        ##
        res$narms <- x$x$n.arms
        res$multiarm <- ifelse(x$x$multiarm, "*", "")
      }
      ##
      res <- as.matrix(res)
      dimnames(res)[[1]] <- x$x$studlab
      ##
      if (!missing.truncate) {
        sortvar <- sortvar[truncate]
        res <- res[truncate, , drop = FALSE]
      }
      ##
      res <- res[order(sortvar), , drop = FALSE]
      dimnames(res)[[1]] <- treats(dimnames(res)[[1]], nchar.studlab)
      prmatrix(res, quote = FALSE, right = TRUE)
      if (!missing.truncate)
        cat(text.truncate, "\n")
      cat("\n")
      ##      
      studyarms <- data.frame(narms = x$x$narms, row.names = x$x$studies)
      if (!missing.truncate)
        studyarms <-
          studyarms[rownames(studyarms) %in% rownames(res), , drop = FALSE]
      cat("Number of treatment arms (by study):\n")
      rownames(studyarms) <- treats(rownames(studyarms), nchar.studlab)
      prmatrix(studyarms, quote = FALSE, right = TRUE)
      if (!missing.truncate)
        cat(text.truncate, "\n")
      cat("\n")
    }
    
    
    ##
    ##
    ## (5) Print results for individual studies
    ##
    ##
    TE.f    <- x$comparison.nma.common$TE
    lowTE.f <- x$comparison.nma.common$lower
    uppTE.f <- x$comparison.nma.common$upper
    ##
    if (backtransf) {
      TE.f    <- backtransf(TE.f, sm)
      lowTE.f <- backtransf(lowTE.f, sm)
      uppTE.f <- backtransf(uppTE.f, sm)
      #
      # Switch lower and upper limit for VE if results have been
      # backtransformed
      #
      if (sm == "VE") {
        tmp.l <- lowTE.f
        lowTE.f <- uppTE.f
        uppTE.f <- tmp.l
      }
    }
    #
    TE.r    <- x$comparison.nma.random$TE
    lowTE.r <- x$comparison.nma.random$lower
    uppTE.r <- x$comparison.nma.random$upper
    ##
    if (backtransf) {
      TE.r    <- backtransf(TE.r, sm)
      lowTE.r <- backtransf(lowTE.r, sm)
      uppTE.r <- backtransf(uppTE.r, sm)
      #
      # Switch lower and upper limit for VE if results have been
      # backtransformed
      #
      if (sm == "VE") {
        tmp.l <- lowTE.r
        lowTE.r <- uppTE.r
        uppTE.r <- tmp.l
      }
    }
    ##
    res.f <- cbind(treat1, treat2,
                   formatN(TE.f, digits, text.NA = "NA", big.mark = big.mark),
                   formatCI(formatN(round(lowTE.f, digits), digits, "NA",
                                    big.mark = big.mark),
                            formatN(round(uppTE.f, digits), digits, "NA",
                                    big.mark = big.mark)),
                   if (common)
                     formatN(round(x$x$Q.common, digits.Q), digits.Q, "NA",
                             big.mark = big.mark),
                   if (common & !all(x$x$narms > 2))
                     formatN(round(x$x$leverage.common, 2), 2, ".")
                   )
    dimnames(res.f) <-
      list(treats(x$x$studlab, nchar.studlab),
           c("treat1", "treat2",
             sm.lab, ci.lab,
             if (common) "Q",
             if (common & !all(x$x$narms > 2)) "leverage"))
    ##
    res.r <- cbind(treat1, treat2,
                   formatN(TE.r, digits, text.NA = "NA", big.mark = big.mark),
                   formatCI(formatN(round(lowTE.r, digits), digits, "NA",
                                    big.mark = big.mark),
                            formatN(round(uppTE.r, digits), digits, "NA",
                                    big.mark = big.mark)))
    dimnames(res.r) <-
      list(treats(x$x$studlab, nchar.studlab),
           c("treat1", "treat2", sm.lab, ci.lab))
    ##
    if (common) {
      cat("Results (common effects model):\n\n")
      ##
      if (!missing.truncate)
        res.f <- res.f[truncate, , drop = FALSE]
      ##
      prmatrix(res.f[order(sortvar), , drop = FALSE],
               quote = FALSE, right = TRUE)
      if (!missing.truncate)
        cat(text.truncate, "\n")      
      cat("\n")
    }
    ##
    if (random) {
      cat("Results (random effects model):\n\n")
      ##
      if (!missing.truncate)
        res.r <- res.r[truncate, , drop = FALSE]
      ##
      prmatrix(res.r[order(sortvar), , drop = FALSE],
               quote = FALSE, right = TRUE)
      if (!missing.truncate)
        cat(text.truncate, "\n")      
      cat("\n")
    }
  }
  ##  
  if (reference.group != "" & missing(all.treatments))
    all.treatments <- FALSE
  ##  
  if (reference.group != "")
    reference.group <- setref(reference.group, rownames(x$x$A.matrix))
  ##  
  if (nma) {
    print.netmeta(x$x,
                  common = common, random = random,
                  prediction = prediction,
                  backtransf = backtransf,
                  overall.hetstat = overall.hetstat,
                  reference.group = reference.group,
                  baseline.reference = baseline.reference,
                  all.treatments = all.treatments,
                  header = FALSE, nchar.trts = nchar.trts,
                  ##
                  digits = digits,
                  digits.pval.Q = digits.pval.Q,
                  digits.Q = digits.Q,
                  digits.tau2 = digits.tau2,
                  digits.I2 = digits.I2,
                  scientific.pval = scientific.pval,
                  big.mark = big.mark,
                  #
                  details.methods = details.methods,
                  legend = legend)
  }
  else {
    ##
    ## Add legend with abbreviated treatment labels
    ##
    legendabbr(trts, treats(trts, nchar.trts), legend)
  }
  
  
  invisible(NULL)
}
