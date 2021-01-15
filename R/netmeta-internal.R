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
