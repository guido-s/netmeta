#' Create a data frame from an object of class netmeta
#' 
#' @description
#' The \code{as.data.frame} method returns a data frame containing
#' information on individual studies, e.g., estimated treatment effect
#' and its standard error.
#' 
#' @param x An object of class \code{netmeta}.
#' @param row.names \code{NULL} or a character vector giving the row
#'   names for the data frame.
#' @param optional A logical. If \code{TRUE}, setting row names and
#'   converting column names (to syntactic names) is optional.
#' @param details A logical. If \code{TRUE}, additional variables of
#'   less interest are included in data frame.
#' @param \dots Additional arguments.
#' @return A data frame is returned by the function
#'   \code{as.data.frame}.
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' @seealso \code{\link{netmeta}}
#' 
#' @examples
#' data(smokingcessation)
#' 
#' # Transform data from arm-based format to contrast-based format
#' #
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' 
#' # Conduct random effects network meta-analysis and show data frame
#' #
#' net1 <- netmeta(p1, common = FALSE)
#' as.data.frame(net1)
#' 
#' \dontrun{
#' data(Senn2013)
#' 
#' # Conduct network meta-analysis
#' #
#' net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD")
#' 
#' as.data.frame(net2)
#' as.data.frame(net2, details = TRUE)
#' }
#'
#' @method as.data.frame netmeta 
#' @export


as.data.frame.netmeta <- function(x, row.names = NULL,
                                  optional = FALSE,
                                  details = FALSE, ...){
  
  
  chkclass(x, "netmeta")
  ##
  x <- updateversion(x)
  
  
  ## Remove element 'call' from object of class meta to get rid
  ## of an error message in meta-analyses with six studies:
  ## 'Error: evaluation nested too deeply: infinite recursion ...'
  ##
  ## NB: Element 'call' which is of length six contains information
  ##     on the function call.
  ##
  x$call <- NULL
  
  sel1 <- as.vector(lapply(x, length) == length(x$studlab))
  sel2 <- as.vector(unlist(lapply(x, is.vector)))
  sel <- sel1 & sel2
  
  res <- as.data.frame(x[names(x)[sel]], ...)
  ##
  res$studies <- NULL
  res$narms <- NULL
  ##
  res$lower.nma.fixed <- NULL
  res$upper.nma.fixed <- NULL
  res$leverage.fixed <- NULL
  
  if (!details)
    res <- res[, !(names(res) %in% c("treat1.pos", "treat2.pos",
                                     "lower.nma.common", "upper.nma.common",
                                     "lower.nma.random", "upper.nma.random",
                                     "leverage.common"))]
  
  attr(res, "version") <- packageDescription("netmeta")$Version
  
  res
}
