#' Print method for objects of class netimpact
#' 
#' @description
#' Print method for objects of class \code{netimpact}.
#' 
#' @param x An object of class \code{netimpact}.
#' @param fixed A logical indicating whether results for the
#'   fixed effects / common effects model should be printed.
#' @param random A logical indicating whether results for the
#'   random effects model should be printed.
#' @param digits Minimal number of significant digits.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
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
#' # Only consider first two studies (to reduce runtime of example)
#' #
#' studies <- unique(Franchini2012$Study)
#' p1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'                n = list(n1, n2, n3),
#'                mean = list(y1, y2, y3),
#'                sd = list(sd1, sd2, sd3),
#'                data = subset(Franchini2012, Study %in% studies[1:2]),
#'                studlab = Study)
#' 
#' net1 <- netmeta(p1)
#' ni <- netimpact(net1, verbose = TRUE)
#' ni
#' 
#' @method print netimpact
#' @export


print.netimpact <- function(x,
                            fixed = x$x$fixed,
                            random = x$x$random,
                            digits = gs("digits.prop"),
                            ##
                            legend = TRUE,
                            warn.deprecated = gs("warn.deprecated"),
                            ##
                            ...) {
  
  ##
  ##
  ## (1) Check for netimpact object and upgrade object
  ##
  ##
  chkclass(x, "netimpact")
  x <- updateversion(x)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chknumeric(digits, min = 0, length = 1)
  chklogical(legend)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  fixed <- deprecated(fixed, missing(fixed), args, "comb.fixed",
                      warn.deprecated)
  chklogical(fixed)
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  
  
  ##
  ##
  ## (3) Generate (abbreviated) column names
  ##
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
  ##
  treat1 <- as.character(factor(treat1.long, levels = trts,
                                labels = treats(trts, x$x$nchar.trts)))
  treat2 <- as.character(factor(treat2.long, levels = trts,
                                labels = treats(trts, x$x$nchar.trts)))
  
  
  ##
  ##
  ## (4) Print results
  ##
  ##
  if (fixed) {
    cat("Fixed effects model: \n\n")
    impact.fixed <- formatN(x$impact.fixed, digits = digits)
    colnames(impact.fixed) <- paste(treat1, treat2, sep = sep.trts)
    ##
    prmatrix(impact.fixed, quote = FALSE, right = TRUE)
    if (random)
      cat("\n")
  }
  ##
  ## Print results for random effects model
  ##
  if (random) {
    cat("Random effects model: \n\n")
    impact.random <- formatN(x$impact.random, digits = digits)
    colnames(impact.random) <- paste(treat1, treat2, sep = sep.trts)
    ##
    prmatrix(impact.random, quote = FALSE, right = TRUE)
  }
  ##
  ## Add legend with abbreviated treatment labels
  ##
  if (fixed | random)
    legendabbr(trts, treats(trts, x$x$nchar.trts), legend)
  
  
  invisible(NULL)
}
