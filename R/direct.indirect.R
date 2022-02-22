direct.indirect <- function(x, tol.direct = 0.0005) {
  
  
  chkclass(x, "netmeta")
  chknumeric(tol.direct, min = 0, length = 1)
  
  
  setNA <- function(x, select) {
    x[select] <- NA
    x
  }
  
  
  ##
  ## Direct and indirect estimates from fixed effects model
  ##
  sel.df <- abs(x$P.fixed) < tol.direct
  ##
  TE.direct.fixed <- setNA(x$TE.direct.fixed, sel.df)
  seTE.direct.fixed <- setNA(x$seTE.direct.fixed, sel.df)
  ##
  sel.if <- abs(x$P.fixed - 1) < tol.direct
  ##
  TE.indirect.fixed <- setNA(x$TE.indirect.fixed, sel.if)
  seTE.indirect.fixed <- setNA(x$seTE.indirect.fixed, sel.if)
  ## Set indirect estimate to network estimate if k = 0
  TE.indirect.fixed[x$A.matrix == 0] <- x$TE.fixed[x$A.matrix == 0]
  seTE.indirect.fixed[x$A.matrix == 0] <- x$seTE.fixed[x$A.matrix == 0]
  diag(TE.indirect.fixed) <- diag(seTE.indirect.fixed) <- NA
  
  
  ##
  ## Direct and indirect estimates from random effects model
  ##
  sel.dr <- abs(x$P.random) < tol.direct
  ##
  TE.direct.random <- setNA(x$TE.direct.random, sel.df)
  seTE.direct.random <- setNA(x$seTE.direct.random, sel.df)
  ##
  sel.ir <- abs(x$P.random - 1) < tol.direct
  ##
  TE.indirect.random <- setNA(x$TE.indirect.random, sel.ir)
  seTE.indirect.random <- setNA(x$seTE.indirect.random, sel.ir)
  ## Set indirect estimate to network estimate if k = 0
  TE.indirect.random[x$A.matrix == 0] <- x$TE.random[x$A.matrix == 0]
  seTE.indirect.random[x$A.matrix == 0] <- x$seTE.random[x$A.matrix == 0]
  diag(TE.indirect.random) <- diag(seTE.indirect.random) <- NA
  
  
  ##
  ## Calculate confidence limits
  ##
  ci.nf <- ci(x$TE.fixed, x$seTE.fixed, x$level.ma)
  ci.df <- ci(TE.direct.fixed, seTE.direct.fixed, x$level.ma)
  ci.if <- ci(TE.indirect.fixed, seTE.indirect.fixed, x$level.ma)
  ##
  ci.nr <- ci(x$TE.random, x$seTE.random, x$level.ma)
  ci.dr <- ci(TE.direct.random, seTE.direct.random, x$level.ma)
  ci.ir <- ci(TE.indirect.random, seTE.indirect.random, x$level.ma)
  ##
  ci.nf$z <- ci.df$z <- ci.if$z <- 
    ci.nr$z <- ci.dr$z <- ci.ir$z <- NULL
  
  
  res <- list(nma.fixed = ci.nf, nma.random = ci.nr,
              direct.fixed = ci.df, direct.random = ci.dr,
              indirect.fixed = ci.if, indirect.random = ci.ir,
              x = x)
  ##
  res
}
