summary.netmeta <- function(object,
                            comb.fixed = object$comb.fixed,
                            comb.random = object$comb.random,
                            prediction = object$prediction,
                            reference.group = object$reference.group,
                            baseline.reference = object$baseline.reference,
                            all.treatments = object$all.treatments,
                            warn = object$warn,
                            ...) {
  
  
  if (!inherits(object, "netmeta"))
    stop("Argument 'object' must be an object of class \"netmeta\"")
  
  
  object <- upgradenetmeta(object)
  
  
  if (length(warn) == 0)
    warn <- TRUE

  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  meta:::chklogical(prediction)
  meta:::chklogical(baseline.reference)
  ##
  cl <- "netmeta()"
  addargs <- names(list(...))
  ##
  fun <- "summary.netmeta"
  ##
  meta:::warnarg("level", addargs, fun, cl)
  meta:::warnarg("level.comb", addargs, fun, cl)
  
  
  ci.lab <- paste(round(100 * object$level.comb, 1), "%-CI", sep = "")
  ##
  ci.comp <- ci(object$TE, object$seTE, object$level)
  ci.comp$studlab <- object$studlab
  ci.comp$treat1 <- object$treat1
  ci.comp$treat2 <- object$treat2
  ##
  ci.comp.nma.fixed <- meta::ci(object$TE.nma.fixed,
                                object$seTE.nma.fixed, object$level)
  ci.comp.nma.fixed$studlab <- object$studlab
  ci.comp.nma.fixed$treat1 <- object$treat1
  ci.comp.nma.fixed$treat2 <- object$treat2
  ci.comp.nma.fixed$leverage <- object$leverage.fixed
  ##
  ci.comp.nma.random <- meta::ci(object$TE.nma.random,
                                 object$seTE.nma.random, object$level)
  ci.comp.nma.random$studlab <- object$studlab
  ci.comp.nma.random$treat1 <- object$treat1
  ci.comp.nma.random$treat2 <- object$treat2
  ##
  ci.f <- list(TE = object$TE.fixed,
               seTE = object$seTE.fixed,
               lower = object$lower.fixed,
               upper = object$upper.fixed,
               z = object$zval.fixed,
               p = object$pval.fixed,
               level = object$level.comb)
  ##
  ci.r <- list(TE = object$TE.random,
               seTE = object$seTE.random,
               lower = object$lower.random,
               upper = object$upper.random,
               z = object$zval.random,
               p = object$pval.random,
               level = object$level.comb)
  ##
  ci.p <- list(TE = NA,
               seTE = object$seTE.predict,
               lower = object$lower.predict,
               upper = object$upper.predict,
               z = NA,
               p = NA,
               level = object$level.predict,
               df = object$df.Q - 1)
  
  
  res <- list(comparison = ci.comp,
              comparison.nma.fixed = ci.comp.nma.fixed,
              comparison.nma.random = ci.comp.nma.random,
              fixed = ci.f, random = ci.r,
              predict = ci.p,
              studies = object$studies,
              narms = object$narms,
              k = object$k, m = object$m, n = object$n, d = object$d,
              Q = object$Q, df.Q = object$df.Q, pval.Q = object$pval.Q,
              Q.heterogeneity = object$Q.heterogeneity,
              df.Q.heterogeneity = object$df.Q.heterogeneity,
              pval.Q.heterogeneity = object$pval.Q.heterogeneity,
              Q.inconsistency = object$Q.inconsistency,
              df.Q.inconsistency = object$df.Q.inconsistency,
              pval.Q.inconsistency = object$pval.Q.inconsistency,
              tau = object$tau, I2 = object$I2,
              tau.preset = object$tau.preset,
              sm = object$sm,
              ci.lab = ci.lab,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              prediction = prediction,
              seq = object$seq
              )
  
  if (reference.group != "" & missing(all.treatments))
    all.treatments <- FALSE
  
  if (reference.group != "")
    reference.group <- setref(reference.group, rownames(object$A.matrix))
  
  res$reference.group <- reference.group
  res$baseline.reference <- baseline.reference
  res$all.treatments <- all.treatments
  ##
  res$nchar.trts <- object$nchar.trts
  ##
  res$backtransf <- object$backtransf
  ##
  res$title <- object$title
  ##
  res$call <- match.call()
  res$version <- packageDescription("netmeta")$Version
  ##
  class(res) <- "summary.netmeta"
  
  res
}
