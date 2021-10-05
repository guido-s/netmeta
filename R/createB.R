createB <- function (pos1, pos2, ncol, aggr = FALSE) {
  
  
  if (!aggr) {
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
  }
  else {
    nrow <- 0
    ##
    ## Determine number of edges (no. of rows of B)
    ##
    for (i in 1:(ncol - 1)) {
      for (j in (i + 1):ncol) {
        ij.count <- 0
        ## Cycle through every possible edge ij
        ## Search pos1 and pos2 to see if at least one of these
        ## combinations is ij
        for (k in seq_along(pos1)) {
          if (pos1[k] == i & pos2[k] == j) {
            ij.count <- ij.count + 1
          }
          else {
            ij.count <- ij.count
          }
        }
        if (ij.count > 0)
          nrow <- nrow + 1
        else
          nrow <- nrow
      }
    }
    ##
    ## Create aggregate B matrix with dimensions e x n
    ##
    B <- matrix(0, nrow = nrow, ncol = ncol)
    ##
    r <- 0
    ## Cycle through each possible pairwise comparison ij
    for (i in 1:(ncol - 1)) {
      for (j in (i + 1):ncol) {
        ij.count <- 0
        for (k in 1:length(pos1)) {
          ## If there is an edge for that pairwise comparison ...
          if (pos1[k] == i & pos2[k] == j)
            ij.count <- ij.count + 1 # ...then ij.count is no longer = 0 ...
          else
            ij.count <- ij.count
        }
        if (ij.count > 0) {
          ## ...and we add this row to B
          r <- r + 1
          B[r, i] <-  1
          B[r, j] <- -1
        }
      }
    }
  }
  
  B
}
