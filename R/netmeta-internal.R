.onAttach <-
function (libname, pkgname) 
{
  msg <- paste("Loading 'netmeta' package (version ",
               packageDescription("netmeta")$Version,
               ").", sep = "")
  packageStartupMessage(msg)
}
