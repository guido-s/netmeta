#' Print method for objects of class netbind
#' 
#' @description
#' Print method for objects of class \code{netbind}.
#' 
#' @param x An object of class \code{netbind} or
#'   \code{summary.netbind}.
#' @param comb.fixed A logical indicating whether results for the
#'   fixed effects (common effects) model should be printed.
#' @param comb.random A logical indicating whether results for the
#'   random effects model should be printed.
#' @param \dots Additional arguments (ignored).
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
#'                 sm = "OR", comb.fixed = FALSE)
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
#' print(nb1, comb.fixed = TRUE)
#' 
#' @method print netbind
#' @export
#' @export print.netbind


print.netbind <- function(x,
                          comb.fixed = x$comb.fixed,
                          comb.random = x$comb.random,
                          ...) {
  
  
  meta:::chkclass(x, "netbind")
  
  
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  

  if (comb.fixed) {
    cat("Fixed effects model\n\n")
    print(x$fixed[, c("name", "treat",
                      "TE", "seTE", "lower", "upper", "zval", "pval")])
    if (comb.random)
      cat("\n")
  }
  ##
  if (comb.random) {
    cat("Random effects model\n\n")
    print(x$random[, c("name", "treat",
                       "TE", "seTE", "lower", "upper", "zval", "pval")])
  }
  
  
  invisible(NULL)
}
