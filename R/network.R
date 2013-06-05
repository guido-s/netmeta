network <- function(TE, seTE,
                    treat1, treat2,
                    treat1.pos, treat2.pos,
                    narms, studlab,
                    data=NULL,
                    sm="",
                    level=0.95, level.comb=0.95,
                    comb.fixed, comb.random){
  
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  ## Catch TE, treat1, treat2, treat1.pos, treat2.pos, seTE, studlab,
  ## narms from data:
  ##
  mf <- match.call()
  mf$level <- mf$level.comb <- NULL
  mf$data <- NULL
  mf[[1]] <- as.name("data.frame")
  mf <- eval(mf, data)

  TE <- mf$TE
  treat1 <- mf$treat1
  treat2 <- mf$treat2
  treat1.pos <- mf$treat1.pos
  treat2.pos <- mf$treat2.pos
  seTE <- mf$seTE
  ##
  if (length(mf$studlab)!=0){
    if (is.factor(mf$studlab))
      studlab <- as.character(mf$studlab)
  }
  else
    studlab <- seq(along=TE)
  ##
  if (length(mf$narms)==0)
    stop("Argument 'narms' has to be provided.")
  else
    narms <- mf$narms
  
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  
  ##
  ## Check for levels of confidence interval
  ##
  if (!is.numeric(level) | length(level)!=1)
    stop("parameter 'level' must be a numeric of length 1")
  if (level <= 0 | level >= 1)
    stop("parameter 'level': no valid level for confidence interval")
  
  ##
  ## Check for levels of confidence interval
  ##
  if (!is.numeric(level.comb) | length(level.comb)!=1)
    stop("parameter 'level.comb' must be a numeric of length 1")
  if (level.comb <= 0 | level.comb >= 1)
    stop("parameter 'level.comb': no valid level.comb for confidence interval")
  
  w.fixed <- 1/seTE^2
  
  m <- length(TE)                        ## Number of pairwise comparisons (edges)
  n <- length(unique(c(treat1, treat2))) ## Number of treatments (vertices)
  W <- diag(w.fixed)                     ## Weighted degree diagonal matrix
  df1 <- 2*sum(1/narms)                  ## Sum of degrees of freedom per study
  
  ##
  ## B is the edge-vertex incidence matrix (m x n)
  ##
  B <- matrix(0, nrow=m, ncol=n)
  for (i in 1:m){
    if (treat1[i]!=treat2[i]) {
      B[i, treat1.pos[i]] <-  1
      B[i, treat2.pos[i]] <- -1
    }
  }
  ##
  ## M is the unweighted Laplacian, D its diagonal,
  ## A is the adjacency matrix
  ##
  M <- t(B)%*%B      ## unweighted Laplacian matrix
  D <- diag(diag(M)) ## diagonal matrix
  A <- D - M         ## adjacency matrix (n x n)
  ##
  ## L is the weighted Laplacian (Kirchhoff) matrix (n x n)
  ## Lplus is its Moore-Penrose pseudoinverse
  ##
  L <- t(B)%*%W%*%B 
  Lplus <- solve(L-1/n)+1/n
  ##
  ## R resistance distance (variance) matrix (n x n)
  ##
  R <- matrix(0, nrow=n, ncol=n) 
  for (i in 1:n) {
    for (j in 1:n) {
      R[i,j] <- Lplus[i,i] + Lplus[j,j] - 2*Lplus[i,j]
    }
  }
  ##
  ## V is the vector of effective variances
  ##
  V <- vector(length=m, mode="numeric")
  for (i in 1:m){
    V[i] <- R[treat1.pos[i], treat2.pos[i]]
  }
  ##
  ## G is the matrix B%*%Lplus%*%t(B) GW is the projection matrix (also
  ## called "hat matrix")
  ##
  ## Interpretation:
  ## (i)    diag(G) = V               The effective variances
  ## (ii)   diag(GW) = V%*%W = V*w    The leverages
  ## (iii)  sum(diag(GW)) = n-1       Rank of projection
  ## (iv)   mean(diag(GW)) = (n-1)/m  Mean leverage = average efficiency
  ##
  G <- B%*%Lplus%*%t(B)
  GW <- G%*%W
  ##
  ## Resulting effects and variances at numbered edges
  ##
  v <- as.vector(GW%*%TE)
  ci.v <- meta:::ci(v, sqrt(V), level=level)
  ##
  ## Resulting effects, all edges, as a n x n matrix:
  ##
  all <- matrix(NA, nrow=n, ncol=n)
  for (i in 1:m){
    all[treat1.pos[i], treat2.pos[i]] <- v[i]
  }
  for(i in 1:n){
    for(j in 1:n){
      for(k in 1:n){
        if(is.na(all[i,k])){
          all[i,k] <- all[i,j]-all[k,j]
        }
        if(is.na(all[k,j])) { 
          all[k,j] <- all[i,j]-all[i,k]
        }
      }
    }
  }
  ##
  ## Test of total heterogeneity/inconsistency:
  ##
  Q <- as.vector(t(TE-v)%*%W%*%(TE-v))
  df <- df1-(n-1)
  ##
  ## Heterogeneity variance and random effects weights
  ##
  tau2 <- max(0, (Q-df)/sum(w.fixed*(1-V*w.fixed)))
  tau <- sqrt(tau2)
  w.random <- 1/(1/w.fixed+tau2)
  ##
  ## Decomposition of total Q into parts from pairwise meta-analyses
  ## and residual inconsistency
  ##
  Q.matrix <- matrix(0, nrow=n, ncol=n)
  n.pairwise <- 0
  data <- data.frame(treat1.pos, treat2.pos,
                     TE, w.fixed,
                     stringsAsFactors=FALSE)
  ##
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (A[i,j] > 1) {
        n.pairwise <- n.pairwise + 1
        sub <- data[(data$treat1.pos==i & data$treat2.pos==j) |
                    (data$treat1.pos==j & data$treat2.pos==i),]
        l <- nrow(sub)
        Q.matrix[i,j] <- metagen(TE=sub$TE, seTE=1/sqrt(sub$w.fixed))$Q
      }
    }
  }
  q <- t1 <- t2 <- dfs <- vector(length=n.pairwise, mode="numeric")
  p <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (A[i,j] > 0) {
        p <- p+1
        t1[p] <- i
        t2[p] <- j
        q[p] <- Q.matrix[i,j]
        dfs[p] <- A[i,j]-1
      }
    }
  }
  ##
  Q.heterogeneity <- sum(Q.matrix)
  Q.inconsistency <- Q - Q.heterogeneity
  ##
  names.treat <- sort(unique(c(treat1, treat2)))
  ##
  Q.decomp <- data.frame(treat1=names.treat[t1],
                         treat2=names.treat[t2],
                         Q=q,
                         df=dfs,
                         pval.Q=1-pchisq(q, dfs))
  
  TE.fixed <- all
  seTE.fixed <- sqrt(R)
  ##
  ci.fixed <- meta:::ci(all, sqrt(R), level=level.comb)
  ##
  lower.fixed <- ci.fixed$lower
  upper.fixed <- ci.fixed$upper
  zval.fixed <- ci.fixed$z
  pval.fixed <- ci.fixed$p
  ##
  rownames(TE.fixed) <- colnames(TE.fixed) <- names.treat
  rownames(seTE.fixed) <- colnames(seTE.fixed) <- names.treat
  rownames(lower.fixed) <- colnames(lower.fixed) <- names.treat
  rownames(upper.fixed) <- colnames(upper.fixed) <- names.treat
  rownames(zval.fixed) <- colnames(zval.fixed) <- names.treat
  rownames(pval.fixed) <- colnames(pval.fixed) <- names.treat

  A.matrix <- A
  L.matrix <- L
  Lplus.matrix <- Lplus
  ##
  rownames(A.matrix) <- colnames(A.matrix) <- names.treat
  rownames(L.matrix) <- colnames(L.matrix) <- names.treat
  rownames(Lplus.matrix) <- colnames(Lplus.matrix) <- names.treat
  rownames(Q.matrix) <- colnames(Q.matrix) <- names.treat
  
  G.matrix <- G
  H.matrix <- GW
  ##
  rownames(G.matrix) <- colnames(G.matrix) <- studlab
  rownames(H.matrix) <- colnames(H.matrix) <- studlab
  
  
  ## Contribution of individual studies to Q
  ##
  Q.fixed <- w.fixed*(TE-v)^2
  
  
  res <- list(
              studlab=studlab,
              treat1=treat1, treat2=treat2,
              TE=TE, seTE=seTE,
              TE.nma.fixed=v,
              seTE.nma.fixed=sqrt(V),
              lower.nma.fixed=ci.v$lower,
              upper.nma.fixed=ci.v$upper,
              TE.nma.random=NA,
              seTE.nma.random=NA,
              lower.nma.random=NA,
              upper.nma.random=NA,
              leverage.fixed=diag(GW),
              leverage.random=NA,
              w.fixed=w.fixed,
              Q.fixed=Q.fixed,
              w.random=w.random,
              ##Q.random=Q.random,
              treat1.pos=treat1.pos,
              treat2.pos=treat2.pos,
              ##
              TE.fixed=TE.fixed,
              seTE.fixed=seTE.fixed,
              lower.fixed=lower.fixed,
              upper.fixed=upper.fixed,
              zval.fixed=zval.fixed,
              pval.fixed=pval.fixed,
              ##
              TE.random=NA,
              seTE.random=NA,
              lower.random=NA,
              upper.random=NA,
              zval.random=NA,
              pval.random=NA,
              ##
              k=length(unique(studlab)),
              m=length(TE),
              Q=Q,
              df=df,
              pval.Q=1-pchisq(Q, df),
              I2=max(0, 100*(Q-df)/Q),
              tau=tau,
              Q.heterogeneity=Q.heterogeneity,
              Q.inconsistency=Q.inconsistency,
              ##
              sm=sm,
              level=level,
              level.comb=level.comb,
              comb.fixed=comb.fixed,
              comb.random=comb.random,
              ##
              A.matrix=A.matrix,
              L.matrix=L.matrix,
              Lplus.matrix=Lplus.matrix,
              Q.matrix=Q.matrix,
              G.matrix=G.matrix,
              H.matrix=H.matrix,
              ##
              Q.decomp = Q.decomp
              )
  
  res$version <- packageDescription("netmeta")$Version
  
  res
}
