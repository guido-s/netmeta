hasse <- function(x, newpage = TRUE) {
  
  meta:::chkclass(x, "netposet")
  
  
  if (!meta:::is.installed.package("hasseDiagram", stop = FALSE))
    stop(paste("Library 'hasseDiagram' missing.",
               "\n  ",
               "Please use the following R commands for installation:",
               "\n  ",
               "source(\"https://bioconductor.org/biocLite.R\")",
               "\n  ",
               "biocLite(\"Rgraphviz\")",
               "\n  ",
               "install.packages(\"hasseDiagram\")",
               sep = ""),
         call. = FALSE)
  
  
  M <- x$M
  n <- length(x$treatments)
  ## Warum nicht M0 aus netposet object ???
  M0 <- matrix(as.logical(M),
               ncol = n, nrow = n,
               dimnames = dimnames(M))
  
  
  hasseDiagram::hasse(M0, parameters = list(newpage = newpage))
  
  
  invisible(NULL)
}
