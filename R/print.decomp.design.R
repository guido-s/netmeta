print.decomp.design <- function(x,
                                digits=2, ...){
  
  if (!inherits(x, "decomp.design"))
    stop("Argument 'x' must be an object of class \"decomp.design\"")

  Q.decomp <- x$Q.decomp
  Q.design <- x$Q.design
  Q.detach <- x$Q.detach
  ##
  Q.decomp$Q <- round(Q.decomp$Q, digits)
  Q.decomp$pval <- meta:::format.p(Q.decomp$pval)
  ##
  Q.design$Q <- round(Q.design$Q, digits)
  Q.design$pval <- meta:::format.p(Q.design$pval)
  ##
  Q.detach$Q <- round(Q.detach$Q, digits)
  Q.detach$pval <- meta:::format.p(Q.detach$pval)
  ##
  nam <- names(Q.decomp)
  names(Q.decomp)[nam=="pval"] <- "p.value"
  nam <- names(Q.design)
  names(Q.design)[nam=="design"] <- "Design"
  names(Q.design)[nam=="pval"] <- "p.value"
  nam <- names(Q.detach)
  names(Q.detach)[nam=="design"] <- "Detached design"
  names(Q.detach)[nam=="pval"] <- "p.value"
  ##
  Q.design <- as.matrix(Q.design)
  Q.detach <- as.matrix(Q.detach)
  
  
  cat("Q statistics to assess homogeneity / consistency\n\n")
  print(Q.decomp)
  
  cat("\nDesign-specific decomposition of within-designs Q statistic\n\n")
  dimnames(Q.design) <- list(rep("", dim(Q.design)[[1]]),
                             colnames(Q.design))
  prmatrix(Q.design, quote=FALSE, right=TRUE)
  
  cat("\nBetween-designs Q statistic after detaching of single designs\n\n")
  dimnames(Q.detach) <- list(rep("", dim(Q.detach)[[1]]),
                             colnames(Q.detach))
  prmatrix(Q.detach, quote=FALSE, right=TRUE)
  
  invisible(NULL)
}
