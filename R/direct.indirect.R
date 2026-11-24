direct.indirect <- function(x, tol.direct = 0.0005) {
  
  chkclass(x, "netmeta")
  chknumeric(tol.direct, min = 0, length = 1)
  #
  method.random.ci <- replaceNULL(x$method.random.ci, "classic")

  ##
  ## Direct and indirect estimates from common effects model
  ##
  sel.direct.c <- abs(x$P.common) < tol.direct
  ##
  TE.direct.common <- setNA_subset(x$TE.direct.common, sel.direct.c)
  seTE.direct.common <- setNA_subset(x$seTE.direct.common, sel.direct.c)
  ##
  sel.indirect.c <- abs(x$P.common - 1) < tol.direct
  ##
  TE.indirect.common <- setNA_subset(x$TE.indirect.common, sel.indirect.c)
  seTE.indirect.common <- setNA_subset(x$seTE.indirect.common, sel.indirect.c)
  ## Set indirect estimate to network estimate if k = 0
  TE.indirect.common[x$A.matrix == 0] <- x$TE.common[x$A.matrix == 0]
  seTE.indirect.common[x$A.matrix == 0] <- x$seTE.common[x$A.matrix == 0]
  diag(TE.indirect.common) <- diag(seTE.indirect.common) <- NA
  
  ##
  ## Direct and indirect estimates from random effects model
  ##
  sel.direct.r <- abs(x$P.random) < tol.direct
  ##
  TE.direct.random <- setNA_subset(x$TE.direct.random, sel.direct.r)
  seTE.direct.random <- setNA_subset(x$seTE.direct.random, sel.direct.r)
  ##
  sel.indirect.r <- abs(x$P.random - 1) < tol.direct
  ##
  TE.indirect.random <- setNA_subset(x$TE.indirect.random, sel.indirect.r)
  seTE.indirect.random <- setNA_subset(x$seTE.indirect.random, sel.indirect.r)
  ## Set indirect estimate to network estimate if k = 0
  TE.indirect.random[x$A.matrix == 0] <- x$TE.random[x$A.matrix == 0]
  seTE.indirect.random[x$A.matrix == 0] <- x$seTE.random[x$A.matrix == 0]
  diag(TE.indirect.random) <- diag(seTE.indirect.random) <- NA
  
  #
  # Calculate confidence limits
  #
  
  # Common effects model
  #
  ci.direct.c <- ci(TE.direct.common, seTE.direct.common, x$level.ma)
  ci.indirect.c <- ci(TE.indirect.common, seTE.indirect.common, x$level.ma)
  ci.mixed.c <- ci(x$TE.common, x$seTE.common, x$level.ma)
  #
  # Random effects model
  #
  if (method.random.ci == "classic") {
    ci.mixed.r <- ci(x$TE.random, x$seTE.random, x$level.ma)
    ci.direct.r <- ci(TE.direct.random, seTE.direct.random, x$level.ma)
    ci.indirect.r <- ci(TE.indirect.random, seTE.indirect.random, x$level.ma)
  }
  else if (method.random.ci == "t-dist") {
    ci.mixed.r <-
      ci(x$TE.random, x$seTE.random, x$level.ma, x$df.Q)
    ci.direct.r <-
      ci(TE.direct.random, seTE.direct.random, x$level.ma, x$df.Q)
    ci.indirect.r <-
      ci(TE.indirect.random, seTE.indirect.random, x$level.ma, x$df.Q)
  }
  else
    stop("Method to calculate confidence intervals is not defined: ",
         "'method.random.ci = ", method.random.ci, "'.",
         call. = FALSE)
  
  res <- list(nma.common = ci.mixed.c, nma.random = ci.mixed.r,
              direct.common = ci.direct.c, direct.random = ci.direct.r,
              indirect.common = ci.indirect.c, indirect.random = ci.indirect.r,
              x = x)
  ##
  res
}
