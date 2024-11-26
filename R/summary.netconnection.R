#' Summary method for objects of class netconnection
#' 
#' @description
#' Summary method for objects of class \code{netconnection} to print
#' list of studies in subnetworks.
#' 
#' @param object An object of class \code{netconnection}.
#' @param x An object of class \code{summary.netconnection}.
#' @param \dots Additional arguments (passed on to
#'   \code{\link{print.netconnection}}.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
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
#' summary(nc2)
#'
#' @method summary netconnection 
#' @export

summary.netconnection <- function(object, ...) {
  
  chkclass(object, "netconnection")
  
  res <- object
  #
  class(res) <- c("summary.netconnection", class(res))
  #
  res
}


#' @rdname summary.netconnection
#' @method print summary.netconnection 
#' @export

print.summary.netconnection <- function(x, ...) {
  
  chkclass(x, "summary.netconnection")
  #
  class(x) <- "netconnection"
  print(x, ...)
  #
  cat("\nStudies in subnetworks\n\n")
  #
  for (i in seq_len(x$n.subnets)) {
    cat(paste0("Subnet ", i, ":\n"))
    cat(paste(unique(x$studlab[x$subnet == i]), collapse = ", "))
    cat("\n")
  }
  #
  invisible(NULL)
}
