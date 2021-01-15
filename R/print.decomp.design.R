#' Print method for objects of class decomp.design
#' 
#' @description
#' Print method for objects of class \code{decomp.design}.
#' 
#' @param x An object of class \code{decomp.design}.
#' @param digits.Q Minimal number of significant digits for Q
#'   statistics, see \code{print.default}.
#' @param showall A logical indicating whether results should be shown
#'   for all designs or only designs contributing to chi-squared
#'   statistics (default).
#' @param digits.pval.Q Minimal number of significant digits for
#'   p-value of heterogeneity tests, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
#' @param \dots Additional arguments (ignored at the moment).
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}, Ulrike
#'   Krahn \email{ulrike.krahn@@bayer.com}
#' 
#' @seealso \code{\link{decomp.design}}
#' 
#' @keywords print
#' 
#' @examples
#' data(Senn2013)
#' 
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'                 data = Senn2013, sm = "MD")
#' print(decomp.design(net1))
#' 
#' @method print decomp.design
#' @export
#' @export print.decomp.design


print.decomp.design <- function(x,
                                digits.Q = gs("digits.Q"),
                                showall = FALSE,
                                digits.pval.Q = gs("digits.pval.Q"),
                                digits.tau2 = gs("digits.tau2"),
                                scientific.pval = gs("scientific.pval"),
                                big.mark = gs("big.mark"), ...) {
  
  
  meta:::chkclass(x, "decomp.design")
  ##
  formatPT <- meta:::formatPT
  formatN <- meta:::formatN
  chknumeric <- meta:::chknumeric
  ##
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.pval.Q, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  ##
  meta:::chklogical(showall)
  
  
  if (!is.null(attributes(x)$netmetabin)) {
    Qdata <- x$Q.decomp
    ##
    Qdata$Q <- round(Qdata$Q, digits.Q)
    Qdata$pval <- meta:::formatPT(Qdata$pval,
                                  digits = digits.pval.Q,
                                  scientific = scientific.pval)
    ##
    dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
    ##
    cat(paste0("\nTest of inconsistency (between designs):\n"))
    prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
    ##
    return(invisible(NULL))
  }
  
  
  Q.decomp <- x$Q.decomp
  Q.design <- x$Q.het.design
  Q.detach <- x$Q.inc.detach
  ##
  if (!showall) {
    Q.design <- Q.design[Q.design$df > 0, ]
    ##
    sel.designs <- Q.detach$design[Q.detach$df !=
                                   Q.decomp["Between designs", "df"]]
    Q.detach <- Q.detach[Q.detach$design %in% sel.designs, ]
  }
  ##
  Q.inc.random <- x$Q.inc.random
  Q.inc.random.preset <- x$Q.inc.random.preset
  ##
  Q.decomp$Q <- formatN(round(Q.decomp$Q, digits.Q), digits.Q, "NA",
                        big.mark = big.mark)
  Q.decomp$pval <- formatPT(Q.decomp$pval, digits = digits.pval.Q,
                            scientific = scientific.pval)
  ##
  Q.design$Q <- formatN(round(Q.design$Q, digits.Q), digits.Q, "NA",
                        big.mark = big.mark)
  Q.design$pval <- formatPT(Q.design$pval, digits = digits.pval.Q,
                            scientific = scientific.pval)
  ##
  Q.detach$Q <- formatN(round(Q.detach$Q, digits.Q), digits.Q, "NA",
                        big.mark = big.mark)
  Q.detach$pval <- formatPT(Q.detach$pval, digits = digits.pval.Q,
                            scientific = scientific.pval)
  ##
  Q.inc.random$Q <- formatN(round(Q.inc.random$Q, digits.Q), digits.Q, "NA",
                            big.mark = big.mark)
  Q.inc.random$tau2.within <- formatPT(Q.inc.random$tau.within^2,
                                       digits = digits.tau2,
                                       lab.NA = "NA", big.mark = big.mark)
  Q.inc.random$tau.within <- formatPT(Q.inc.random$tau.within,
                                      digits = digits.tau2,
                                      lab.NA = "NA", big.mark = big.mark)
  Q.inc.random$pval <- formatPT(Q.inc.random$pval, digits = digits.pval.Q,
                                scientific = scientific.pval)
  ##
  nam <- names(Q.decomp)
  names(Q.decomp)[nam == "pval"] <- "p-value"
  nam <- names(Q.design)
  names(Q.design)[nam == "design"] <- "Design"
  names(Q.design)[nam == "pval"] <- "p-value"
  nam <- names(Q.detach)
  names(Q.detach)[nam == "design"] <- "Detached design"
  names(Q.detach)[nam == "pval"] <- "p-value"
  nam <- names(Q.inc.random)
  names(Q.inc.random)[nam == "pval"] <- "p-value"
  ##
  Q.design <- as.matrix(Q.design)
  Q.detach <- as.matrix(Q.detach)
  
  cat("Q statistics to assess homogeneity / consistency\n\n")
  print(Q.decomp)
  
  if (nrow(Q.design) > 0) {
    cat("\nDesign-specific decomposition of within-designs Q statistic\n\n")
    dimnames(Q.design) <- list(rep("", dim(Q.design)[[1]]),
                               colnames(Q.design))
    prmatrix(Q.design, quote = FALSE, right = TRUE)
  }
  
  if (nrow(Q.detach) > 0) {
    cat("\nBetween-designs Q statistic after detaching of single designs\n\n")
    dimnames(Q.detach) <- list(rep("", dim(Q.detach)[[1]]),
                               colnames(Q.detach))
    prmatrix(Q.detach, quote = FALSE, right = TRUE)
  }

  cat(paste("\nQ statistic to assess consistency under the assumption of\n",
            "a full design-by-treatment interaction random effects model\n\n",
            sep = ""))
  print(Q.inc.random)

  invisible(NULL)
}
