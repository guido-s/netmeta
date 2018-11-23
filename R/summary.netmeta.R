summary.netmeta <- function(object,
                            comb.fixed = object$comb.fixed,
                            comb.random = object$comb.random,
                            prediction = object$prediction,
                            reference.group = object$reference.group,
                            baseline.reference = object$baseline.reference,
                            all.treatments = object$all.treatments,
                            warn = object$warn,
                            ...) {
  
  ##
  ##
  ## (1) Check for netmeta object and upgrade older meta objects
  ##
  ##
  meta:::chkclass(object, "netmeta")
  ##  
  object <- upgradenetmeta(object)
  ##  
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
  
  
  ##
  ##
  ## (3) Summarise results for individual studies and network
  ##     meta-analyses
  ##
  ##
  keepvars <- c("TE", "seTE", "lower", "upper", "z", "p")
  ##
  ci.comp <- data.frame(studlab = object$studlab,
                        treat1 = object$treat1, treat2 = object$treat2,
                        ci(object$TE, object$seTE, object$level)[keepvars],
                        stringsAsFactors = FALSE)
  ##
  ci.nma.fixed <- data.frame(studlab = object$studlab,
                             treat1 = object$treat1,
                             treat2 = object$treat2,
                             TE = object$TE.nma.fixed,
                             seTE = object$seTE.nma.fixed,
                             lower = object$lower.nma.fixed,
                             upper = object$upper.nma.fixed,
                             z = object$zval.nma.fixed,
                             p = object$pval.nma.fixed,
                             leverage = object$leverage.fixed,
                             stringsAsFactors = FALSE)
  ##
  ci.nma.random <- data.frame(studlab = object$studlab,
                             treat1 = object$treat1,
                             treat2 = object$treat2,
                             TE = object$TE.nma.random,
                             seTE = object$seTE.nma.random,
                             lower = object$lower.nma.random,
                             upper = object$upper.nma.random,
                             z = object$zval.nma.random,
                             p = object$pval.nma.random,
                             stringsAsFactors = FALSE)
  ##
  ci.f <- list(TE = object$TE.fixed,
               seTE = object$seTE.fixed,
               lower = object$lower.fixed,
               upper = object$upper.fixed,
               z = object$zval.fixed,
               p = object$pval.fixed)
  ##
  ci.r <- list(TE = object$TE.random,
               seTE = object$seTE.random,
               lower = object$lower.random,
               upper = object$upper.random,
               z = object$zval.random,
               p = object$pval.random)
  ##
  ci.p <- list(seTE = object$seTE.predict,
               lower = object$lower.predict,
               upper = object$upper.predict)
  
  
  ##
  ##
  ## (4) Create summary.netmeta object
  ##
  ##
  res <- list(comparison = ci.comp,
              comparison.nma.fixed = ci.nma.fixed,
              comparison.nma.random = ci.nma.random,
              fixed = ci.f, random = ci.r,
              predict = ci.p,
              ##
              studies = object$studies,
              narms = object$narms,
              ##
              k = object$k, m = object$m, n = object$n, d = object$d,
              ##
              Q = object$Q,
              df.Q = object$df.Q,
              pval.Q = object$pval.Q,
              I2 = object$I2,
              tau = object$tau,
              ##
              Q.heterogeneity = object$Q.heterogeneity,
              df.Q.heterogeneity = object$df.Q.heterogeneity,
              pval.Q.heterogeneity = object$pval.Q.heterogeneity,
              ##
              Q.inconsistency = object$Q.inconsistency,
              df.Q.inconsistency = object$df.Q.inconsistency,
              pval.Q.inconsistency = object$pval.Q.inconsistency,
              ##
              Q.decomp = object$Q.decomp,
              ##
              sm = object$sm,
              method = object$method,
              level = object$level,
              level.comb = object$level.comb,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              ##
              prediction = prediction,
              level.predict = object$level.predict,
              ##
              incr = object$incr,
              allincr = object$allincr,
              addincr = object$addincr,
              allstudies = object$allstudies,
              cc.pooled = object$cc.pooled,
              ##
              ci.lab = paste0(round(100 * object$level.comb, 1),"%-CI"),
              ##
              reference.group = NA,
              baseline.reference = NA,
              all.treatments = NA,
              seq = object$seq,
              ##
              tau.preset = object$tau.preset,
              ##
              sep.trts = object$sep.trts,
              nchar.trts = object$nchar.trts,
              ##
              backtransf = object$backtransf,
              ##
              title = object$title,
              ##
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  if (reference.group != "" & missing(all.treatments))
    all.treatments <- FALSE
  ##
  if (reference.group != "")
    reference.group <- setref(reference.group, rownames(object$A.matrix))
  ##
  res$reference.group <- reference.group
  res$baseline.reference <- baseline.reference
  res$all.treatments <- all.treatments
  ##
  class(res) <- "summary.netmeta"
  
  res
}
