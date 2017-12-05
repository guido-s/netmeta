summary.netcomb <- function(object, ...) {
  
  meta:::chkclass(object, "netcomb")

  res <- object
  ##
  class(res) <- c("summary.netcomb", "netcomb")
  
  res
}
