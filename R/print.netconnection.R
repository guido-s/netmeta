print.netconnection <- function(x,
                                digits = max(4, .Options$digits - 3),
                                nchar.trts = x$nchar.trts, ...) {
  
  meta:::chkclass(x, "netconnection")
  ##
  if (is.null(nchar.trts))
    nchar.trts <- 666
  
  
  meta:::chknumeric(digits, single = TRUE)
  meta:::chknumeric(nchar.trts, min = 1, single = TRUE)
  
  
  matitle(x)
  ##
  cat(paste("Number of studies: k = ", x$k, "\n", sep = ""))
  cat(paste("Number of treatments: n = ", x$n, "\n", sep = ""))
  cat(paste("Number of pairwise comparisons: m = ", x$m, "\n", sep = ""))
  ##
  cat("Number of networks: ", x$n.subnets, "\n\n", sep = "")
  
  cat("Distance matrix:\n")
  
  D <- round(x$D.matrix, digits = digits)
  D[is.infinite(D)] <- "."
  ##
  if (x$n.subnets == 1)
    diag(D) <- "."
  ##
  rownames(D) <- treats(rownames(D), nchar.trts)
  colnames(D) <- treats(colnames(D), nchar.trts)
  ##
  prmatrix(D, quote = FALSE, right = TRUE)
  
  
  if (any(rownames(x$D.matrix) != rownames(D))) {
    abbr <- rownames(D)
    full <- rownames(x$D.matrix)
    ##
    tmat <- data.frame(abbr, full)
    names(tmat) <- c("Abbreviation", "Treatment name")
    tmat <- tmat[order(tmat$Abbreviation), ]
    ##
    cat("\nLegend:\n")
    prmatrix(tmat, quote = FALSE, right = TRUE,
             rowlab = rep("", length(abbr)))
  }
  
  
  invisible(NULL)
}
