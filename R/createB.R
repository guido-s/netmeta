createB <- function (pos1, pos2, ncol) {
  if (missing(pos1) | missing(pos2)) {
    ##
    ## Create full edge-vertex incidence matrix
    ##
    nrow <- choose(ncol, 2)
    B <- matrix(0, nrow = nrow, ncol = ncol)
    ##
    i <- 0
    ##
    for (pos1.i in 1:(ncol - 1)) {
      for (pos2.i in (pos1.i + 1):ncol) {
        i <- i + 1
        B[i, pos1.i] <-  1
        B[i, pos2.i] <- -1
      }
    }
  }
  else {
    ##
    ## Create edge-vertex incidence matrix
    ##
    nrow <- length(pos1)
    ncol <- length(unique(c(pos1, pos2)))
    ##
    B <- matrix(0, nrow = nrow, ncol = ncol)
    ##
    for (i in 1:nrow) {
      B[i, pos1[i]] <-  1
      B[i, pos2[i]] <- -1
    }
  }
  
  
  B
}
