print.decomp.design <- function(x,
                                digits = 2, showall = FALSE,
                                ...) {
  
  
  meta:::chkclass(x, "decomp.design")
  meta:::chklogical(showall)
  
  
  Q.decomp <- x$Q.decomp
  Q.design <- x$Q.het.design
  Q.detach <- x$Q.inc.detach
  ##
  if (!showall) {
    sel.designs <- Q.detach$design[Q.detach$df !=
                                   Q.decomp["Between designs", "df"]]
    Q.design <- Q.design[Q.design$design %in% sel.designs, ]
    Q.detach <- Q.detach[Q.detach$design %in% sel.designs, ]
  }
  ##
  Q.inc.random <- x$Q.inc.random
  Q.inc.random.preset <- x$Q.inc.random.preset
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
  Q.inc.random$Q <- round(Q.inc.random$Q, digits)
  Q.inc.random$tau2.within <- meta:::format.tau(Q.inc.random$tau.within^2)
  Q.inc.random$tau.within  <- meta:::format.tau(Q.inc.random$tau.within)
  Q.inc.random$pval <- meta:::format.p(Q.inc.random$pval)
  ##
  nam <- names(Q.decomp)
  names(Q.decomp)[nam == "pval"] <- "p-value"
  nam <- names(Q.design)
  names(Q.design)[nam == "design"] <- "Design"
  names(Q.design)[nam == "pval"] <- "p-value"
  nam <- names(Q.detach)
  names(Q.detach)[nam == "design"] <- "Detached design"
  names(Q.detach)[nam == "pval"] <- "p-value"
  nam <- names(Q.inc.random)
  names(Q.inc.random)[nam == "pval"] <- "p-value"
  ##
  Q.design <- as.matrix(Q.design)
  Q.detach <- as.matrix(Q.detach)
  
  cat("Q statistics to assess homogeneity / consistency\n\n")
  print(Q.decomp)
  
  if (nrow(Q.design) > 0) {
    cat("\nDesign-specific decomposition of within-designs Q statistic\n\n")
    dimnames(Q.design) <- list(rep("", dim(Q.design)[[1]]),
                               colnames(Q.design))
    prmatrix(Q.design, quote = FALSE, right = TRUE)
  }
  
  if (nrow(Q.detach) > 0) {
    cat("\nBetween-designs Q statistic after detaching of single designs\n\n")
    dimnames(Q.detach) <- list(rep("", dim(Q.detach)[[1]]),
                               colnames(Q.detach))
    prmatrix(Q.detach, quote = FALSE, right = TRUE)
  }

  cat(paste("\nQ statistic to assess consistency under the assumption of\n",
            "a full design-by-treatment interaction random effects model\n\n",
            sep = ""))
  print(Q.inc.random)

  invisible(NULL)
}
