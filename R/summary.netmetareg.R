#' Summary method for objects of class netmetareg
#' 
#' @description
#' Summary method for objects of class \code{netmetareg} to print
#' list of studies in subnetworks.
#' 
#' @param object An object of class \code{netmetareg}.
#' @param x An object of class \code{summary.netmetareg}.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard
#'   errors.
#' @param big.mark A character used as thousands separator.
#' @param \dots Additional arguments (ignored).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' @seealso \code{\link{netmetareg}}
#' 
#' @examples
#' \dontrun{
#' data(smokingcessation)
#' # Add variable with (fictitious) risk of bias values
#' # with 1 = "low risk" and 2 = "high risk"
#' #
#' smokingcessation$rob <- rep(1:2, 12)
#' 
#' pw1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' 
#' net1 <- netmeta(pw1, common = FALSE, ref = "A")
#' 
#' # Network meta-regression with continuous covariate and assumption of
#' # independent slopes
#' nr1 <- netmetareg(net1, rob)
#' nr1
#' 
#' summary(nr1)
#' }
#'
#' @method summary netmetareg 
#' @export

summary.netmetareg <- function(object, ...) {
  
  chkclass(object, "netmetareg")
  
  res <- object
  #
  class(res) <- c("summary.netmetareg", class(res))
  #
  res
}


#' @rdname summary.netmetareg
#' @method print summary.netmetareg 
#' @export

print.summary.netmetareg <- function(x,
                                     digits = gs("digits"),
                                     digits.se = gs("digits.se"),
                                     big.mark = gs("big.mark"),
                                     ...) {
  
  chkclass(x, "summary.netmetareg")
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chkchar(big.mark, length = 1)
  #
  dat <- x$full_results
  rownames(dat) <- dat$comparison
  dat$comparison <- NULL
  #
  dat$d <-
    formatN(dat$d, digits = digits, big.mark = big.mark, text.NA = ".")
  dat$se.d <-
    formatN(dat$se.d, digits = digits.se, big.mark = big.mark, text.NA = ".")
  dat$beta <-
    formatN(dat$beta, digits = digits, big.mark = big.mark, text.NA = ".")
  dat$se.beta <-
    formatN(dat$se.beta, digits = digits.se, big.mark = big.mark, text.NA = ".")
  dat$cov <-
    formatN(dat$cov, digits = digits.se, big.mark = big.mark, text.NA = ".")
  #
  # Do not print treatment labels (which are contained in the row names)
  #
  dat$treat1 <- dat$treat2 <- NULL
  #
  prmatrix(dat, quote = FALSE, right = TRUE)
  #
  invisible(NULL)
}
