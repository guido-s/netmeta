summary.netcomb <- function(object,
                            comb.fixed = object$comb.fixed,
                            comb.random = object$comb.random,
                            ...) {
  
  ##
  ##
  ## (1) Check for netcomb object
  ##
  ##
  meta:::chkclass(object, "netcomb")
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  
  
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
  ci.cnma.fixed <- data.frame(studlab = object$studlab,
                              treat1 = object$treat1,
                              treat2 = object$treat2,
                              TE = object$TE.cnma.fixed,
                              seTE = object$seTE.cnma.fixed,
                              lower = object$lower.cnma.fixed,
                              upper = object$upper.cnma.fixed,
                              z = object$zval.cnma.fixed,
                              p = object$pval.cnma.fixed,
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
  ci.cnma.random <- data.frame(studlab = object$studlab,
                               treat1 = object$treat1,
                               treat2 = object$treat2,
                               TE = object$TE.cnma.random,
                               seTE = object$seTE.cnma.random,
                               lower = object$lower.cnma.random,
                               upper = object$upper.cnma.random,
                               z = object$zval.cnma.random,
                               p = object$pval.cnma.random,
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
  ci.comp.f <- data.frame(comps = object$comps,
                          TE = object$Comp.fixed,
                          seTE = object$seComp.fixed,
                          lower = object$lower.Comp.fixed,
                          upper = object$upper.Comp.fixed,
                          z = object$zval.Comp.fixed,
                          p = object$pval.Comp.fixed,
                          stringsAsFactors = FALSE)
  ##
  ci.comp.r <- data.frame(comps = object$comps,
                          TE = object$Comp.random,
                          seTE = object$seComp.random,
                          lower = object$lower.Comp.random,
                          upper = object$upper.Comp.random,
                          z = object$zval.Comp.random,
                          p = object$pval.Comp.random,
                          stringsAsFactors = FALSE)
  ##
  ci.comb.f <- data.frame(trts = object$trts,
                          TE = object$Comb.fixed,
                          seTE = object$seComb.fixed,
                          lower = object$lower.Comb.fixed,
                          upper = object$upper.Comb.fixed,
                          z = object$zval.Comb.fixed,
                          p = object$pval.Comb.fixed,
                          stringsAsFactors = FALSE)
  ##
  ci.comb.r <- data.frame(trts = object$trts,
                          TE = object$Comb.random,
                          seTE = object$seComb.random,
                          lower = object$lower.Comb.random,
                          upper = object$upper.Comb.random,
                          z = object$zval.Comb.random,
                          p = object$pval.Comb.random,
                          stringsAsFactors = FALSE)
  
  
  ##
  ##
  ## (4) Create summary.netmeta object
  ##
  ##
  res <- list(k = object$k,
              m = object$m,
              n = object$n,
              d = object$d,
              c = object$c,
              ##
              trts = object$trts,
              k.trts = object$k.trts,
              n.trts = object$n.trts,
              events.trts = object$events.trts,
              ##
              studies = object$studies,
              narms = object$narms,
              ##
              designs = object$designs,
              ##
              comps = object$comps,
              k.comps = object$k.comps,
              n.comps = object$n.comps,
              events.comps = object$events.comps,
              ##
              comparison = ci.comp,
              comparison.nma.fixed = ci.nma.fixed,
              comparison.nma.random = ci.nma.random,
              comparison.cnma.fixed = ci.cnma.fixed,
              comparison.cnma.random = ci.cnma.random,
              ##
              components.fixed = ci.comp.f,
              components.random = ci.comp.r,
              ##
              combinations.fixed = ci.comb.f,
              combinations.random = ci.comb.r,
              ##
              fixed = ci.f, random = ci.r,
              ##
              Q.additive = object$Q.additive,
              df.Q.additive = object$df.Q.additive,
              pval.Q.additive = object$pval.Q.additive,
              tau = object$tau,
              I2 = object$I2,
              ##
              Q.standard = object$Q.standard,
              df.Q.standard = object$df.Q.standard,
              pval.Q.standard = object$pval.Q.standard,
              ##
              Q.diff = object$Q.diff,
              df.Q.diff = object$df.Q.diff,
              pval.Q.diff = object$pval.Q.diff, 
              ##
              sm = object$sm,
              method = object$method,
              level = object$level,
              level.comb = object$level.comb,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
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
  class(res) <- "summary.netcomb"
  
  res
}
