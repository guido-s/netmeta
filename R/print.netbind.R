print.netbind <- function(x,
                          comb.fixed = x$comb.fixed,
                          comb.random = x$comb.random,
                          ...) {
  
  
  meta:::chkclass(x, "netbind")
  
  
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  

  if (comb.fixed) {
    cat("Fixed effects model\n\n")
    print(x$fixed[, c("name", "treat",
                      "TE", "seTE", "lower", "upper", "zval", "pval")])
    if (comb.random)
      cat("\n")
  }
  ##
  if (comb.random) {
    cat("Random effects model\n\n")
    print(x$random[, c("name", "treat",
                       "TE", "seTE", "lower", "upper", "zval", "pval")])
  }
  
  
  invisible(NULL)
}
