#' Print method for objects of class netbind
#' 
#' @description
#' Print method for objects of class \code{netbind}.
#' 
#' @param x An object of class \code{netbind} or
#'   \code{summary.netbind}.
#' @param fixed A logical indicating whether results for the fixed
#'   effects / common effects model should be printed.
#' @param random A logical indicating whether results for the random
#'   effects model should be printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netbind}}
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
#' # Standard random effects NMA model (with placebo as reference
#' # treatment)
#' #
#' net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'                 data = face, reference.group = "placebo",
#'                 sm = "OR", fixed = FALSE)
#' 
#' # Additive CNMA model with placebo as inactive component and
#' # reference
#' #
#' nc1 <- netcomb(net1, inactive = "placebo")
#' 
#' # Combine results of standard NMA and CNMA
#' #
#' nb1 <- netbind(nc1, net1,
#'                name = c("Additive CNMA", "Standard NMA"),
#'                col.study = c("red", "black"),
#'                col.square = c("red", "black"))
#' 
#' nb1
#' print(nb1, fixed = TRUE)
#' 
#' @method print netbind
#' @export


print.netbind <- function(x,
                          fixed = x$x$fixed,
                          random = x$x$random,
                          ##
                          warn.deprecated = gs("warn.deprecated"),
                          ##
                          ...) {
  
  ##
  ##
  ## (1) Check for netbind object and upgrade object
  ##
  ##
  chkclass(x, "netbind")
  x <- updateversion(x)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
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
  ## (3) Print results
  ##
  ##
  if (fixed) {
    cat("Fixed effects model\n\n")
    print(x$fixed[, c("name", "treat",
                      "TE", "seTE", "lower", "upper", "statistic", "pval")])
    if (random)
      cat("\n")
  }
  ##
  if (random) {
    cat("Random effects model\n\n")
    print(x$random[, c("name", "treat",
                       "TE", "seTE", "lower", "upper", "statistic", "pval")])
  }
  
  invisible(NULL)
}
