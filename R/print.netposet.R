print.netposet <- function(x,
                           pooled = ifelse(x$comb.random, "random", "fixed"),
                           ...) {
  
  pooled <- meta:::setchar(pooled, c("fixed", "random"))
  
  if (pooled == "fixed") {
    if (!(!x$comb.fixed & !x$comb.random))
      cat("Fixed effect model\n")
    print(x$M0.fixed, ...)
    if (x$comb.random)
      cat("\n")
  }
  ##
  else {
    cat("Random effects model\n")
    print(x$M0.random, ...)
  }
  
  
  ## diag(M0) <- NA
  ## M0[is.na(M0)] <- "."
  ## prmatrix(M0, quote = FALSE, right = TRUE)
  
  invisible(NULL)
}
