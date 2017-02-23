print.netposet <- function(x, ...) {
  M0 <- x$M0
  print(M0, ...)
  ## diag(M0) <- NA
  ## M0[is.na(M0)] <- "."
  ## prmatrix(M0, quote = FALSE, right = TRUE)
  ##
  invisible(NULL)
}
