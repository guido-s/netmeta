summary.netmeta <- function(object,
                            all.treatments=object$all.treatments,
                            reference.group=object$reference.group,
                            level=object$level,
                            level.comb=object$level.comb,
                            comb.fixed=object$comb.fixed,
                            comb.random=object$comb.random,
                            warn=object$warn,
                            ...){
  
  
  if (!inherits(object, "netmeta"))
    stop("Argument 'object' must be an object of class \"netmeta\"")
  
  if (length(warn)==0){
    warn <- TRUE
  }
  
  k <- object$k
  m <- object$m
  Q <- object$Q
  
  ci.lab <- paste(round(100*level.comb, 1), "%-CI", sep="")
  ##
  ci.comp <- ci(object$TE, object$seTE, level)
  ci.comp.nma.fixed <- ci(object$TE.nma.fixed,
                          object$seTE.nma.fixed, level)
  ci.f <- ci(object$TE.fixed , object$seTE.fixed , level.comb)
  ci.r <- ci(object$TE.random, object$seTE.random, level.comb)
  
  ci.comp$studlab <- object$studlab
  ci.comp$treat1 <- object$treat1
  ci.comp$treat2 <- object$treat2
  ##
  ci.comp.nma.fixed$studlab <- object$studlab
  ci.comp.nma.fixed$treat1 <- object$treat1
  ci.comp.nma.fixed$treat2 <- object$treat2
  ci.comp.nma.fixed$leverage <- object$leverage.fixed
  
  
  res <- list(
              comparison=ci.comp,
              comparison.nma.fixed=ci.comp.nma.fixed,
              fixed=ci.f, random=ci.r,
              studies=object$studies,
              narms=object$narms,
              k=k, m=m, Q=Q, df=object$df,
              tau=object$tau, I2=object$I2,
              sm=object$sm,
              call=match.call(),
              ci.lab=ci.lab,
              comb.fixed=comb.fixed,
              comb.random=comb.random,
              level=level,
              level.comb=level.comb
              )
  
  
  if (reference.group!="" & missing(all.treatments))
    all.treatments <- FALSE
  
  res$all.treatments <- all.treatments
  res$reference.group <- reference.group
  ##
  res$title   <- object$title
  ##
  class(res) <- "summary.netmeta"
  
  res$version <- packageDescription("netmeta")$Version
  
  res
}
