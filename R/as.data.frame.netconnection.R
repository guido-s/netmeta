#' Create a data frame from an object of class netconnection
#' 
#' @description
#' The \code{as.data.frame} method returns a data frame containing
#' information on membership of studies / pairwise comparisons to a
#' (sub)network.
#' 
#' @param x An object of class \code{netconnection}.
#' @param \dots Additional arguments (ignored).
#' @return A data frame is returned by the function
#'   \code{as.data.frame}.
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' @seealso \code{\link{netconnection}}
#' 
#' @examples
#' # Artificial example with two subnetworks
#' #
#' t1 <- c("G", "B", "B", "D", "A", "F")
#' t2 <- c("B", "C", "E", "E", "H", "A")
#' #
#' nc2 <- netconnection(t1, t2)
#' print(nc2, details = TRUE)
#' 
#' as.data.frame(nc2)
#'
#' @method as.data.frame netconnection 
#' @export


as.data.frame.netconnection <- function(x, ...){
  
  
  meta:::chkclass(x, "netconnection")
  
  
  ## Drop unnecessary list elements
  ## 
  res <- data.frame(treat1 = x$treat1,
                    treat2 = x$treat2,
                    studlab = x$studlab,
                    design = x$design,
                    subnet = x$subnet)
  ##
  attr(res, "version") <- x$version
  
  res
}
