#' Print method for objects of class netimpact
#' 
#' @description
#' Print method for objects of class \code{netimpact}.
#' 
#' @param x An object of class \code{netimpact}.
#' @param comb.fixed A logical indicating whether results for the
#'   fixed effects (common effects) model should be printed.
#' @param comb.random A logical indicating whether results for the
#'   random effects model should be printed.
#' @param digits Minimal number of significant digits.
#' @param \dots Additional arguments (ignored).
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netimpact}}
#' 
#' @keywords print
#' 
#' @examples
#' data(Franchini2012)
#'
#' # Only consider first four studies (to reduce runtime of example)
#' #
#' studies <- unique(Franchini2012$Study)
#' p1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'                n = list(n1, n2, n3),
#'                mean = list(y1, y2, y3),
#'                sd = list(sd1, sd2, sd3),
#'                data = subset(Franchini2012, Study %in% studies[1:4]),
#'                studlab = Study)
#' 
#' net1 <- netmeta(p1)
#' ni <- netimpact(net1, verbose = TRUE)
#' ni
#' 
#' @method print netimpact
#' @export


print.netimpact <- function(x,
                            comb.fixed = x$x$comb.fixed,
                            comb.random = x$x$comb.random,
                            digits = gs("digits.prop"), ...) {
  
  meta:::chkclass(x, "netimpact")
  
  
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  meta:::chknumeric(digits, min = 0, length = 1)
  
  
  ##
  ## Generate (abbreviated) column names
  ##
  sep.trts <- x$x$sep.trts
  ##
  cn <- colnames(x$impact.fixed)
  mat <- matrix(unlist(strsplit(cn, split = sep.trts)),
                ncol = 2, byrow = TRUE)
  treat1.long <- mat[, 1]
  treat2.long <- mat[, 2]
  ##
  trts <- x$x$trts
  trts.abbr <- treats(trts, x$x$nchar.trts)
  ##
  treat1 <- as.character(factor(treat1.long, levels = trts, labels = trts.abbr))
  treat2 <- as.character(factor(treat2.long, levels = trts, labels = trts.abbr))
  ##
  if (any(treat1 != treat1.long) | any(treat2 != treat2.long))
    abbr <- c(treat1, treat2)
  else
    abbr <- NULL
  
  
  ##
  ## Print results for fixed effects model
  ##
  if (comb.fixed) {
    cat("Fixed effects model: \n\n")
    impact.fixed <- meta:::formatN(x$impact.fixed, digits = digits)
    colnames(impact.fixed) <- paste(treat1, treat2, sep = sep.trts)
    ##
    prmatrix(impact.fixed, quote = FALSE, right = TRUE)
    if (comb.random)
      cat("\n")
  }
  ##
  ## Print results for random effects model
  ##
  if (comb.random) {
    cat("Random effects model: \n\n")
    impact.random <- meta:::formatN(x$impact.random, digits = digits)
    colnames(impact.random) <- paste(treat1, treat2, sep = sep.trts)
    ##
    prmatrix(impact.random, quote = FALSE, right = TRUE)
  }
  
  
  if (comb.fixed || comb.random)
    if (!is.null(abbr)) {
      abbr <- unique(abbr)
      full <- unique(c(treat1.long, treat2.long))
      ##
      tmat <- data.frame(abbr, full)
      names(tmat) <- c("Abbreviation", "Treatment name")
      tmat <- tmat[abbr != full, ]
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat("\nLegend:\n")
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(abbr)))
    }
  
  
  invisible(NULL)
}
