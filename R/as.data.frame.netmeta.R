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
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' @seealso \code{\link{netmeta}}
#' 
#' @examples
#' data(Senn2013)
#' 
#' # Conduct network meta-analysis
#' #
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD")
#' 
#' as.data.frame(net1)
#' as.data.frame(net1, details = TRUE)
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
  
  if (!details)
    res <- res[, !(names(res) %in% c("treat1.pos", "treat2.pos",
                                     "lower.nma.fixed", "upper.nma.fixed",
                                     "lower.nma.random", "upper.nma.random",
                                     "leverage.fixed"))]
  
  attr(res, "version") <- packageDescription("netmeta")$Version
  
  res
}
