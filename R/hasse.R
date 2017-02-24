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
  
  M <- x$M.matrix
  ##
  hasseDiagram::hasse(matrix(as.logical(M), dim(M), dimnames = dimnames(M)),
                      parameters = list(newpage = newpage))
  
  
  invisible(NULL)
}
