#' Print results of network meta-analysis
#' 
#' @description
#' Print method for objects of class \code{netmeta}.
#' 
#' @param x An object of class \code{netmeta}.
#' @param sortvar An optional vector used to sort individual studies
#'   (must be of same length as \code{x$TE}).
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
#' @param details A logical indicating whether further details for
#'   individual studies should be printed.
#' @param ma A logical indicating whether summary results of
#'   meta-analysis should be printed.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
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
#' @param \dots Additional arguments.
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
#' 
#' \dontrun{
#' # Conduct random effects network meta-analysis
#' #
#' net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'                 data = Senn2013, sm = "MD",
#'                 comb.fixed = FALSE)
#' print(net2, ref = "plac", digits = 3)
#' }
#' 
#' @method print netmeta
#' @export
#' @export print.netmeta


print.netmeta <- function(x,
                          sortvar,
                          comb.fixed = x$comb.fixed, comb.random = x$comb.random,
                          prediction = x$prediction,
                          reference.group = x$reference.group,
                          baseline.reference = x$baseline.reference,
                          all.treatments = x$all.treatments,
                          details = TRUE, ma = TRUE,
                          ##
                          backtransf = x$backtransf, nchar.trts = x$nchar.trts,
                          digits = gs("digits"),
                          digits.se = gs("digits.se"),
                          digits.pval.Q = max(gs("digits.pval.Q"), 2),
                          digits.Q = gs("digits.Q"),
                          digits.tau2 = gs("digits.tau2"),
                          digits.I2 = gs("digits.I2"),
                          scientific.pval = gs("scientific.pval"),
                          big.mark = gs("big.mark"),
                          ...
                          ) {
  
  
  meta:::chkclass(x, "netmeta")
  ##  
  x <- upgradenetmeta(x)
  ##
  formatN <- meta:::formatN
  formatCI <- meta:::formatCI
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  
  
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
  chklogical(baseline.reference)
  ##
  chklogical(backtransf)
  chknumeric(nchar.trts, min = 1, single = TRUE)
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.se, min = 0, single = TRUE)
  chknumeric(digits.pval.Q, min = 1, single = TRUE)
  chknumeric(digits.Q, min = 0, single = TRUE)
  chknumeric(digits.tau2, min = 0, single = TRUE)
  chknumeric(digits.I2, min = 0, single = TRUE)
  ##
  chklogical(scientific.pval)
  
  
  ##
  ## Additional arguments
  ##
  fun <- "print.netmeta"
  addargs <- names(list(...))
  ##
  meta:::warnarg("logscale", addargs, fun, otherarg = "backtransf")
  
  
  if (!inherits(x, "netmetabin")) {
    k.all <- length(x$TE)
    ##
    if (missing(sortvar)) sortvar <- 1:k.all
    ##
    if (length(sortvar) != k.all)
      stop("'x' and 'sortvar' have different length")
    ##
    ci.lab <- paste(round(100 * x$level, 1), "%-CI", sep = "")
    
    sm <- x$sm
    
    sm.lab <- sm
    ##
    if (!backtransf & meta:::is.relative.effect(sm))
      sm.lab <- paste("log", sm, sep = "")
    
    
    trts <- x$trts
    trts.abbr <- treats(trts, nchar.trts)
    ##
    treat1 <- as.character(factor(x$treat1, levels = trts, labels = trts.abbr))
    treat2 <- as.character(factor(x$treat2, levels = trts, labels = trts.abbr))
    ##
    if (any(treat1 != x$treat1) | any(treat2 != x$treat2))
      abbr <- c(treat1, treat2)
    else
      abbr <- NULL
    
    
    matitle(x)
    
    
    if (details) {

      multiarm <- any(x$narms > 2)

      cat(paste("Original data",
                ifelse(multiarm,
                       " (with adjusted standard errors for multi-arm studies)",
                       ""),
                ":\n\n", sep = ""))
      
      res <- data.frame(treat1,
                        treat2,
                        TE = formatN(x$TE, digits, text.NA = "NA",
                                     big.mark = big.mark),
                        seTE = formatN(x$seTE, digits.se, text.NA = "NA",
                                       big.mark = big.mark))
      ##
      if (multiarm) {
        res$seTE.adj <- format(round(x$seTE.adj, digits.se))
        res$narms <- x$n.arms
        res$multiarm <- ifelse(x$multiarm, "*", "")
      }
      ##
      res <- as.matrix(res)
      dimnames(res)[[1]] <- x$studlab
      
      prmatrix(res[order(sortvar), , drop = FALSE],
               quote = FALSE, right = TRUE)
      cat("\n")
      
      cat("Number of treatment arms (by study):\n")
      prmatrix(data.frame(narms = x$narms, row.names = x$studies),
               quote = FALSE, right = TRUE)
      cat("\n")
    }
    
    
    tsum <- summary(x, warn = FALSE)
    ##
    TE.f    <- tsum$comparison.nma.fixed$TE
    lowTE.f <- tsum$comparison.nma.fixed$lower
    uppTE.f <- tsum$comparison.nma.fixed$upper
    ##
    if (backtransf & meta:::is.relative.effect(sm)) {
      TE.f    <- exp(TE.f)
      lowTE.f <- exp(lowTE.f)
      uppTE.f <- exp(uppTE.f)
    }
    ##
    ##
    TE.r    <- tsum$comparison.nma.random$TE
    lowTE.r <- tsum$comparison.nma.random$lower
    uppTE.r <- tsum$comparison.nma.random$upper
    ##
    if (backtransf & meta:::is.relative.effect(sm)) {
      TE.r    <- exp(TE.r)
      lowTE.r <- exp(lowTE.r)
      uppTE.r <- exp(uppTE.r)
    }
    
    
    res.f <- cbind(treat1, treat2,
                   formatN(TE.f, digits, text.NA = "NA", big.mark = big.mark),
                   formatCI(formatN(round(lowTE.f, digits), digits, "NA",
                                    big.mark = big.mark),
                            formatN(round(uppTE.f, digits), digits, "NA",
                                    big.mark = big.mark)),
                   if (comb.fixed)
                     formatN(round(x$Q.fixed, digits.Q), digits.Q, "NA",
                             big.mark = big.mark),
                   if (comb.fixed) format(round(x$leverage.fixed, 2)))
    dimnames(res.f) <-
      list(x$studlab, c("treat1", "treat2",
                        sm.lab, ci.lab,
                        if (comb.fixed) "Q",
                        if (comb.fixed) "leverage"))
    
    
    res.r <- cbind(treat1, treat2,
                   formatN(TE.r, digits, text.NA = "NA", big.mark = big.mark),
                   formatCI(formatN(round(lowTE.r, digits), digits, "NA",
                                    big.mark = big.mark),
                            formatN(round(uppTE.r, digits), digits, "NA",
                                    big.mark = big.mark)))
    dimnames(res.r) <-
      list(x$studlab, c("treat1", "treat2", sm.lab, ci.lab))
    
    
    if (comb.fixed) {
      cat("Results (fixed effects model):\n\n")
      
      prmatrix(res.f[order(sortvar), , drop = FALSE],
               quote = FALSE, right = TRUE)
      
      cat("\n")
    }
    
    if (comb.random) {
      cat("Results (random effects model):\n\n")
      
      prmatrix(res.r[order(sortvar), , drop = FALSE],
               quote = FALSE, right = TRUE)
      
      cat("\n")
    }
  }
  else {
    tsum <- summary(x, warn = FALSE)
    abbr <- NULL
  }
  
  
  if (reference.group != "" & missing(all.treatments))
    all.treatments <- FALSE
  
  
  if (reference.group != "")
    reference.group <- setref(reference.group, rownames(x$A.matrix))
  
  
  if (ma)
    print(tsum,
          comb.fixed = comb.fixed, comb.random = comb.random,
          prediction = prediction,
          backtransf = backtransf,
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
          big.mark = big.mark)
  else
    if (!is.null(abbr)) {
      abbr <- unique(abbr)
      full <- unique(c(x$treat1, x$treat2))
      ##
      tmat <- data.frame(abbr, full)
      names(tmat) <- c("Abbreviation", "Treatment name")
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat("Legend:\n")
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(abbr)))
    }
  
  invisible(NULL)
}
