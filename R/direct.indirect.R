direct.indirect <- function(x, tol.direct = 0.0005) {
  
  chkclass(x, "netmeta")
  chknumeric(tol.direct, min = 0, length = 1)
    
  ##
  ## Direct and indirect estimates from common effects model
  ##
  sel.df <- abs(x$P.common) < tol.direct
  ##
  TE.direct.common <- setNA_subset(x$TE.direct.common, sel.df)
  seTE.direct.common <- setNA_subset(x$seTE.direct.common, sel.df)
  ##
  sel.if <- abs(x$P.common - 1) < tol.direct
  ##
  TE.indirect.common <- setNA_subset(x$TE.indirect.common, sel.if)
  seTE.indirect.common <- setNA_subset(x$seTE.indirect.common, sel.if)
  ## Set indirect estimate to network estimate if k = 0
  TE.indirect.common[x$A.matrix == 0] <- x$TE.common[x$A.matrix == 0]
  seTE.indirect.common[x$A.matrix == 0] <- x$seTE.common[x$A.matrix == 0]
  diag(TE.indirect.common) <- diag(seTE.indirect.common) <- NA
  
  ##
  ## Direct and indirect estimates from random effects model
  ##
  sel.dr <- abs(x$P.random) < tol.direct
  ##
  TE.direct.random <- setNA_subset(x$TE.direct.random, sel.df)
  seTE.direct.random <- setNA_subset(x$seTE.direct.random, sel.df)
  ##
  sel.ir <- abs(x$P.random - 1) < tol.direct
  ##
  TE.indirect.random <- setNA_subset(x$TE.indirect.random, sel.ir)
  seTE.indirect.random <- setNA_subset(x$seTE.indirect.random, sel.ir)
  ## Set indirect estimate to network estimate if k = 0
  TE.indirect.random[x$A.matrix == 0] <- x$TE.random[x$A.matrix == 0]
  seTE.indirect.random[x$A.matrix == 0] <- x$seTE.random[x$A.matrix == 0]
  diag(TE.indirect.random) <- diag(seTE.indirect.random) <- NA
  
  ##
  ## Calculate confidence limits
  ##
  ci.nf <- ci(x$TE.common, x$seTE.common, x$level.ma)
  ci.df <- ci(TE.direct.common, seTE.direct.common, x$level.ma)
  ci.if <- ci(TE.indirect.common, seTE.indirect.common, x$level.ma)
  ##
  ci.nr <- ci(x$TE.random, x$seTE.random, x$level.ma)
  ci.dr <- ci(TE.direct.random, seTE.direct.random, x$level.ma)
  ci.ir <- ci(TE.indirect.random, seTE.indirect.random, x$level.ma)
  ##
  ci.nf$z <- ci.df$z <- ci.if$z <- 
    ci.nr$z <- ci.dr$z <- ci.ir$z <- NULL
  
  res <- list(nma.common = ci.nf, nma.random = ci.nr,
              direct.common = ci.df, direct.random = ci.dr,
              indirect.common = ci.if, indirect.random = ci.ir,
              x = x)
  ##
  res
}
