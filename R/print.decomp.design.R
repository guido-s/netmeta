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
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param sort A logical indicating whether to sort results by
#'   p-values.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param \dots Additional arguments (ignored at the moment).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}, Ulrike
#'   Krahn \email{ulrike.krahn@@bayer.com}
#' 
#' @seealso \code{\link{decomp.design}}
#' 
#' @keywords print
#' 
#' @examples
#' data(Senn2013)
#' 
#' # Only consider first five studies (to reduce runtime of example)
#' #
#' studies <- unique(Senn2013$studlab)
#' Senn2013.5 <- subset(Senn2013, studlab %in% studies[1:5])
#' 
#' # Conduct network meta-analysis with placebo as reference treatment
#' #
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013.5, sm = "MD", reference = "plac")
#' 
#' # Decomposition of Cochran's Q
#' #
#' decomp.design(net1)
#' 
#' @method print decomp.design
#' @export


print.decomp.design <- function(x,
                                digits.Q = gs("digits.Q"),
                                showall = FALSE,
                                digits.pval.Q = gs("digits.pval.Q"),
                                digits.tau2 = gs("digits.tau2"),
                                scientific.pval = gs("scientific.pval"),
                                big.mark = gs("big.mark"),
                                nchar.trts = x$nchar.trts,
                                sort = TRUE,
                                legend = TRUE,
                                ...) {
  
  
  chkclass(x, "decomp.design")
  x <- updateversion(x)
  ##
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.pval.Q, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  ##
  chklogical(showall)
  if (is.null(nchar.trts))
    nchar.trts <- 666
  chknumeric(nchar.trts, min = 1, length = 1)
  chklogical(sort)
  chklogical(legend)
  ##
  if (is.null(x$x$sep.trts))
    sep.trts <- ":"
  else
    sep.trts <- x$x$sep.trts
  
  
  if (!is.null(attributes(x)$netmetabin)) {
    Qdata <- x$Q.decomp
    ##
    Qdata$Q <- round(Qdata$Q, digits.Q)
    Qdata$pval <- formatPT(Qdata$pval,
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
  trts1 <- trts2 <- ""
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
  if (sort)
    Q.design <- Q.design[order(Q.design$pval), , drop = FALSE]
  Q.design$pval <- formatPT(Q.design$pval, digits = digits.pval.Q,
                            scientific = scientific.pval)
  ##
  Q.detach$Q <- formatN(round(Q.detach$Q, digits.Q), digits.Q, "NA",
                        big.mark = big.mark)
  if (sort)
    Q.detach <- Q.detach[order(-Q.detach$pval), , drop = FALSE]
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
    ##
    trts1 <- unique(sort(unlist(compsplit(Q.design[, 1], sep.trts))))
    Q.design[, 1] <- comps(Q.design[, 1], trts1, sep.trts, nchar.trts)
    ##
    dimnames(Q.design) <- list(rep("", dim(Q.design)[[1]]),
                               colnames(Q.design))
    prmatrix(Q.design, quote = FALSE, right = TRUE)
  }
  
  if (nrow(Q.detach) > 0) {
    cat(paste0("\nBetween-designs Q statistic after detaching of ",
               "single designs\n",
               "(influential designs have p-value markedly different from ",
               rmSpace(Q.decomp[3, 3]), ")\n\n"))
    ##
    trts2 <- unique(sort(unlist(compsplit(Q.detach[, 1], sep.trts))))
    Q.detach[, 1] <- comps(Q.detach[, 1], trts2, sep.trts, nchar.trts)
    ##
    dimnames(Q.detach) <- list(rep("", dim(Q.detach)[[1]]),
                               colnames(Q.detach))
    prmatrix(Q.detach, quote = FALSE, right = TRUE)
  }

  cat(paste("\nQ statistic to assess consistency under the assumption of\n",
            "a full design-by-treatment interaction random effects model\n\n",
            sep = ""))
  print(Q.inc.random)
  
  
  ##
  ## Add legend with abbreviated treatment labels
  ##
  trts <- unique(c(trts1, trts2))
  legendabbr(trts, treats(trts, nchar.trts), legend)
  
  
  invisible(NULL)
}
