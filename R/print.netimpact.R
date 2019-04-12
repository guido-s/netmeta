#' Print method for objects of class netimpact
#' 
#' @description
#' Print method for objects of class \code{netimpact}.
#' 
#' @param x An object of class \code{netimpact}.
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
#' data(parkinson)
#' 
#' p1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'                n = list(n1, n2, n3),
#'                mean = list(y1, y2, y3),
#'                sd = list(sd1, sd2, sd3),
#'                data = parkinson, studlab = Study)
#' 
#' net1 <- netmeta(p1)
#' ni <- netimpact(net1, verbose = TRUE)
#' ni
#' 
#' @method print netimpact
#' @export
#' @export print.netimpact


print.netimpact <- function(x, digits = gs("digits.prop"), ...) {
  
  meta:::chkclass(x, "netimpact")
  
  
  meta:::chknumeric(digits, min = 0, single = TRUE)
  
  
  impact <- x$impact
  sep.trts <- x$x$sep.trts
  ##
  cn <- colnames(impact)
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
  colnames(impact) <- paste(treat1, treat2, sep = sep.trts)
  
  
  impact <- meta:::formatN(impact, digits = digits)
  ##
  prmatrix(impact, quote = FALSE, right = TRUE)
  
  
  if (!is.null(abbr)) {
    abbr <- unique(abbr)
    full <- unique(c(treat1.long, treat2.long))
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
