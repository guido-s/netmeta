stress <- function(x,
                   A.matrix=x$A.matrix,
                   N.matrix=sign(A.matrix),
                   ##
                   start.layout="eigen",
                   iterate=TRUE,
                   eig1=2, eig2=3,
                   tol=0.0001,
                   maxit=500,
                   ##
                   allfigures=FALSE,
                   ...){
  
  
  ##
  ## Function stress for optimising 2D representation of a graph,
  ## given its adjacency matrix A.matrix
  ##
  
  
  ##
  ## Calculate Laplacian matrix
  ##
  As.matrix <- sign(A.matrix)
  n <- dim(As.matrix)[1]                          # Dimension
  e <- rep(1, n)                                  # Vector of ones
  L <- diag(as.vector(As.matrix%*%e)) - As.matrix # Laplacian
  Lt <- solve(L-e%*%t(e)/n) + e%*%t(e)/n          # Its pseudoinverse
  
  
  ##
  ## Calculate eigenvectors
  ##
  if (start.layout=="eigen"){
    Eig <- eigen(L, symmetric=TRUE)
    eigvectors <- Eig$vectors
  }
  else if (start.layout=="prcomp"){
    Eig <- prcomp(L, scale=TRUE)
    eigvectors <- Eig$rotation
  }
  
  
  ##
  ## Weight matrix and its pseudoinverse from complete graph of n
  ## vertices
  ##
  K <- 1 - diag(rep(1,n))               # K = Complete graph of n vertices
  W <- diag(rep(n-1,n))-K               # Weight matrix = Laplacian of K
  Wt <- solve(W-e%*%t(e)/n)+e%*%t(e)/n  # Its pseudoinverse
  
  
  ##
  ## Starting layout
  ##
  if (start.layout=="circle")
    coord <- cbind(cos(2*pi*(1:n)/n), sin(2*pi*(1:n)/n))/sqrt(n)
  else if (start.layout=="eigen" | start.layout=="prcomp")
    coord <- cbind(eigvectors[,n-eig1+1], eigvectors[,n-eig2+1])
  else if (start.layout=="random")
    coord <- cbind(rnorm(n, 0, 1/sqrt(n-1)), rnorm(n, 0, 1/sqrt(n-1)))
  
  
  D.matrix <- nodedist(N.matrix)
  
  
  if (allfigures)
    netgraph(x,
             xpos=coord[,1], ypos=coord[,2],
             iterate=iterate,
             allfigures=FALSE, ...)
  
  
  ## Loop
  ##
  convergence <- 10L
  ##
  if (!iterate){
    coord.new <- coord
    convergence <- NA
  }
  else{
    n.iter <- 0
    ##
    while (convergence!=0L & convergence!=1L){
      n.iter <- n.iter+1
      ##
      Lwd <- matrix(nrow=n, ncol=n)
      ##Lwd
      for (i in 1:(n-1)){
        for (j in (i+1):n){
          Lwd[i,j] <- Lwd[j,i] <- W[i,j]*D.matrix[i,j]/sqrt((coord[i,1] - coord[j,1])^2+
                                                            (coord[i,2] - coord[j,2])^2)
          ##
          if (is.na(Lwd[i,j]))
            Lwd[i,j] <- Lwd[j,i] <- 0
          else if (Lwd[i,j]==Inf)
            Lwd[i,j] <- Lwd[j,i] <- 1
          else if (Lwd[i,j]==-Inf)
            Lwd[i,j] <- Lwd[j,i] <- -1
        }
      }
      ##
      for (i in 1:n)
        Lwd[i,i] <- -sum(Lwd[i,], na.rm=TRUE)
      ##
      RS <- Lwd%*%coord
      coord.new <- Wt%*%RS
      coord.new <- cbind(coord.new[,1]/sqrt(sum(coord.new[,1]^2)),
                         coord.new[,2]/sqrt(sum(coord.new[,2]^2)))
      ##
      if (allfigures)
        netgraph(x,
                 xpos=coord.new[,1], ypos=coord.new[,2],
                 iterate=iterate,
                 allfigures=FALSE,
                 ...)
      ##
      maxdiff <- max(coord-coord.new)
      ##
      if (maxdiff < tol)
        convergence <- 0L
      else if (maxit == n.iter)
        convergence <- 1L
      ##
      coord <- coord.new
    }
  }
  
  
  res <- data.frame(labels=rownames(As.matrix),
                    x=coord.new[,1], y=coord.new[,2])
  ##
  attr(res, "convergence") <- convergence
  if (!is.na(convergence) & convergence != 10L){
    attr(res, "n.iter") <- n.iter
    attr(res, "tol") <- tol
    attr(res, "maxdiff") <- maxdiff
    }
  ##
  ##if (!is.na(convergence) & convergence == 1L)
  ## print(attributes(res))
  ##
  res
}
