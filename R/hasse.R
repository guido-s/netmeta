hasse <- function(x,
                  pooled = ifelse(x$comb.random, "random", "fixed"),
                  newpage = TRUE) {
  
  meta:::chkclass(x, "netposet")
  ##
  pooled <- meta:::setchar(pooled, c("fixed", "random"))
  
  
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
  
  if (pooled == "fixed")
    M <- x$M.fixed
  else
    M <- x$M.random
  ##
  hasseDiagram::hasse(matrix(as.logical(M), dim(M), dimnames = dimnames(M)),
                      parameters = list(newpage = newpage))
  
  
  invisible(NULL)
}
