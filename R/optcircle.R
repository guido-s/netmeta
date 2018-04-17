optcircle <- function(x, start.layout = "eigen") {
  
  n <- x$n
  n1 <- netgraph(x, start.layout = start.layout,
                 iterate = TRUE, figure = FALSE)
  
  
  xpos <- n1$nodes$xpos
  ypos <- n1$nodes$ypos
  ##
  slope <- ypos / xpos
  
  
  angle <- rep(0, n)
  ##
  for (i in 1:n) {
    if (xpos[i] >= 0 & ypos[i] >= 0)
      angle[i] <- atan(slope[i])
    else if (xpos[i] >= 0 & ypos[i] <  0)
      angle[i] <- atan(slope[i]) + 2 * pi 
    else if (xpos[i] < 0)
      angle[i] <- atan(slope[i]) + pi 
  }
  ##
  seq <- order(angle)
  ##
  ## Permutation matrix P
  ##
  P <- diag(n)
  ##
  for (i in 1:n) {
    k <- seq[i]
    if (k != i) {
      P[i, i] <- P[k, k] <- 0
      P[i, k] <- 1
    }
  }
  
  
  ##
  ## Distance matrix D of n-circle, here used to penalize distances
  ##
  D <- t(P) %*% netconnection(LETTERS[1:n], LETTERS[c(2:n, 1)])$D.matrix %*% P
  ##
  A.matrix <- sign(x$A.matrix) * D # A.matrix is the adjacency matrix, edges penalized
  ##
  distsum <- sum(A.matrix) / 2 # Weighted distance sum
  
  
  res <- list(xpos = xpos, ypos = ypos, seq = seq,
              A.matrix = A.matrix, distsum = distsum,
              start.layout = start.layout, x = x)
  ##
  res
}
