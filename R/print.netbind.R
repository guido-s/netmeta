#' Print method for objects of class netbind
#' 
#' @description
#' Print method for objects of class \code{netbind}.
#' 
#' @param x An object of class \code{netbind} or
#'   \code{summary.netbind}.
#' @param common A logical indicating whether results for the common
#'   effects model should be printed.
#' @param random A logical indicating whether results for the random
#'   effects model should be printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
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
#'   data = face, reference.group = "placebo",
#'   sm = "OR", common = FALSE)
#' 
#' # Additive CNMA model with placebo as inactive component and
#' # reference
#' #
#' nc1 <- netcomb(net1, inactive = "placebo")
#' 
#' # Combine results of standard NMA and CNMA
#' #
#' nb1 <- netbind(nc1, net1,
#'   name = c("Additive CNMA", "Standard NMA"),
#'   col.study = c("red", "black"), col.square = c("red", "black"))
#' 
#' nb1
#' print(nb1, common = TRUE)
#' 
#' @method print netbind
#' @export

print.netbind <- function(x,
                          common = x$x$common,
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
  
  
  ##
  ##
  ## (3) Print results
  ##
  ##
  
  nam <- c("name", "treat",
           "TE", "seTE", "lower", "upper", "statistic", "pval")
  #
  if (common) {
    cnames <- nam
    if (!is.null(x$common$method) && length(unique(x$common$method)) > 1)
      cnames <- c(cnames, "method")
    #
    cat("Common effects model\n\n")
    print(x$common[, cnames])
    if (random)
      cat("\n")
  }
  ##
  if (random) {
    rnames <- nam
    if (!is.null(x$random$method) && length(unique(x$random$method)) > 1)
      rnames <- c(rnames, "method")
    #
    cat("Random effects model\n\n")
    print(x$random[, rnames])
  }
  
  invisible(NULL)
}
