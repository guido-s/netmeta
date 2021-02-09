.onAttach <-
function (libname, pkgname) 
{
  msg <- paste("Loading 'netmeta' package (version ",
               packageDescription("netmeta")$Version,
               ").",
               "\nType 'help(\"netmeta-package\")' for a brief overview.",
               sep = "")
  packageStartupMessage(msg)
}


.special.characters <- c("+", ".", "&", "$", "#", "|", "*", "^")


is.zero <- function(x, n = 10)
  abs(x) < n * .Machine$double.eps


invmat <- function(X) {
  n <- nrow(X)
  J <- matrix(1, nrow = n, ncol = n)
  ##
  res <- solve(X - J / n) + J / n
  ##
  res
}
