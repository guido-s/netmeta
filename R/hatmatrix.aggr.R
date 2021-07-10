hatmatrix.aggr <- function(x, model, type = "full") {
  
  model <- meta:::setchar(model, c("fixed", "random"))
  type <- meta:::setchar(type, c("full", "long", "short"))
  
  ## Create aggregate B matrix
  if (!is.null(x$B.matrix.aggr))
    B <- x$B.matrix.aggr
  else
    B <- createB(x$treat1.pos, x$treat2.pos, x$n, aggr = TRUE)
  ## Number of edges (direct comparisons) = no. of rows of B
  e <- nrow(B)
  ##
  ## Aggregate weights
  ##
  ## Create a matrix of the aggregate weights e x e
  ##
  W <- matrix(0, nrow = e, ncol = e)
  for (i in 1:e){
    ## idxA element of row i of B that corresponds to baseline treatment
    ## idxB element of row i of B that corresponds to the other treatment
    for (j in seq_len(x$n))
      if (B[i, j] == 1)
        idxA <- j
      else if (B[i, j] == -1)
        idxB <- j
    ## Find all rows of x$B.matrix (dimension m x n) that have the 
    ## same elements as row i of B and calculate the inverse variance 
    ## of the weighted mean
    WAB <- 0.0
    for (k in seq_len(x$m)) {
        if (x$B.matrix[k, idxA] == 1 & x$B.matrix[k, idxB] == -1) {
          if (model == "fixed")
            WAB <- WAB + 1.0 / (x$seTE.adj.fixed[k])^2
          else if (model == "random")
            WAB <- WAB + 1.0 / (x$seTE.adj.random[k])^2
        }
    }
    ##
    W[i, i] = WAB
  }
  ##
  ## Laplacian of the aggregate model
  ##
  L <- t(B) %*% W %*% B
  ## Pseudo-Inverse of L 
  L.plus <- invmat(L)
  L.plus[is.zero(L.plus)] <- 0
    ##
  ## Aggregate Hat matrix
  ##
  ## type = "short"
  ## e x e hat matrix
  H.short <- B %*% L.plus %*% t(B) %*% W
  ##
  ## type = "long"
  ## n (n - 1) / 2 x e hat matrix
  ## Get edge incidence matrix for a fully connected network with n
  ## nodes
  B.full <- createB(ncol = x$n)
  H.long <- B.full %*% L.plus %*% t(B) %*% W
  ## type = "full"
  ## n (n - 1) / 2 x n (n - 1) / 2 hat matrix
  ## Insert columns of zeroes for missing edges
  c <- (x$n) * (x$n - 1.0) / 2.0
  H.full <- matrix(0, nrow = c, ncol = c)
  ##
  k <- 0 # counts H.full rows and columns
  l <- 1 # counts B rows and H.long columns
  for (i in 1:(x$n - 1)) {
    for (j in (i + 1):x$n){
      k <- k + 1
      ## is there a direct comparison of ij?
      if (l <= e){
        ## once l is > e (no. of edges) we have reached the end of B
        if (B[l, i] == 1 & B[l, j] == -1){
          for (m in 1:c)
            H.full[m, k] <- H.long[m, l]
          l <- l+1
        }
      }
    }
  }
  ##    
  if (type == "short")
    H <- H.short
  else if (type == "long")
    H <- H.long
  else if (type == "full")
    H <- H.full
  ##
  H
}
