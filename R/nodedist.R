nodedist <- function(A){
  ## Distance algorithm by Mueller, Knop, Szymanski and Trinajstic
  ## using triangle inequality
  ##
  ## Calculates the distance matrix D of a given adjacency matrix A
  
  n <- dim(A)[1]

  ## Starting value for D is A, with 0 replaced by n
  
  D <- A
  ##
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (D[i,j]==0) {
        D[i,j] <- n
        D[j,i] <- n
      }
    }
  }
  ##
  for (d in 1:(n-1)) {
    for (i in 1:n) {
      for (j in 1:n) {
        if (D[i,j]==d) {
          for (k in 1:n) {
            akj <- D[k,i] + d   # = D[k,i] + D[i,j]
            D[k,j] <- min(D[k,j],akj)
          }
        }
      }
    } 
  }
  D
}
