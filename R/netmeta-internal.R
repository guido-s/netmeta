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
