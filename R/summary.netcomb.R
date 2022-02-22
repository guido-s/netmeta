#' Summary method for objects of class netcomb
#' 
#' @description
#' Summary method for objects of class \code{netcomb}.
#' 
#' @param object An object of class \code{netcomb}.
#' @param fixed A logical indicating whether results for the fixed
#'   effects / common effects model should be printed.
#' @param random A logical indicating whether results for the random
#'   effects model should be printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#'
#' @return
#' A list is returned with the same elements as a
#' \code{\link{netcomb}} object.
#'
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netcomb}}, \code{\link{discomb}}
#' 
#' @keywords print
#' 
#' @examples
#' data(Linde2016)
#' 
#' # Only consider studies including Face-to-face PST (to reduce
#' # runtime of example)
#' #
#' face <- subset(Linde2016, id %in% c(16, 24, 49, 118))
#' 
#' # Conduct random effects network meta-analysis
#' #
#' net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'   data = face, reference.group = "placebo",
#'   sm = "OR", fixed = FALSE)
#' 
#' # Additive model for treatment components
#' #
#' nc1 <- netcomb(net1)
#' summary(nc1)
#' print(summary(nc1), digits = 2, digits.stat = 3)
#' 
#' \dontrun{
#' # Conduct random effects network meta-analysis
#' #
#' net2 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'   data = Linde2016, reference.group = "placebo",
#'   sm = "OR", fixed = FALSE)
#' 
#' # Additive model for treatment components
#' #
#' nc2 <- netcomb(net2)
#' summary(nc2)
#' print(summary(nc2), digits = 2, digits.stat = 3)
#' }
#' 
#' @method summary netcomb
#' @export


summary.netcomb <- function(object,
                            fixed = object$fixed,
                            random = object$random,
                            warn.deprecated = gs("warn.deprecated"),
                            ...) {
  
  ##
  ##
  ## (1) Check for netcomb object
  ##
  ##
  chkclass(object, "netcomb")
  object <- updateversion(object)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  missing.fixed <- missing(fixed)
  fixed <- deprecated(fixed, missing.fixed, args, "comb.fixed",
                      warn.deprecated)
  chklogical(fixed)
  ##
  random <-
    deprecated(random, missing(random), args, "comb.random", warn.deprecated)
  chklogical(random)
  
  
  ##
  ##
  ## (3) Summarise results for individual studies and network
  ##     meta-analyses
  ##
  ##
  keepvars <- c("TE", "seTE", "lower", "upper", "statistic", "p")
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
                             statistic = object$statistic.nma.fixed,
                             p = object$pval.nma.fixed,
                             stringsAsFactors = FALSE)
  ##
  ci.cnma.fixed <- data.frame(studlab = object$studlab,
                              treat1 = object$treat1,
                              treat2 = object$treat2,
                              TE = object$TE.cnma.fixed,
                              seTE = object$seTE.cnma.fixed,
                              lower = object$lower.cnma.fixed,
                              upper = object$upper.cnma.fixed,
                              statistic = object$statistic.cnma.fixed,
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
                              statistic = object$statistic.nma.random,
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
                               statistic = object$statistic.cnma.random,
                               p = object$pval.cnma.random,
                               stringsAsFactors = FALSE)
  ##
  ci.f <- list(TE = object$TE.fixed,
               seTE = object$seTE.fixed,
               lower = object$lower.fixed,
               upper = object$upper.fixed,
               statistic = object$statistic.fixed,
               p = object$pval.fixed)
  ##
  ci.r <- list(TE = object$TE.random,
               seTE = object$seTE.random,
               lower = object$lower.random,
               upper = object$upper.random,
               statistic = object$statistic.random,
               p = object$pval.random)
  ##
  ci.comp.f <- data.frame(TE = object$Comp.fixed,
                          seTE = object$seComp.fixed,
                          lower = object$lower.Comp.fixed,
                          upper = object$upper.Comp.fixed,
                          statistic = object$statistic.Comp.fixed,
                          p = object$pval.Comp.fixed,
                          stringsAsFactors = FALSE)
  rownames(ci.comp.f) <- object$comps
  ##
  ci.comp.r <- data.frame(TE = object$Comp.random,
                          seTE = object$seComp.random,
                          lower = object$lower.Comp.random,
                          upper = object$upper.Comp.random,
                          statistic = object$statistic.Comp.random,
                          p = object$pval.Comp.random,
                          stringsAsFactors = FALSE)
  rownames(ci.comp.r) <- object$comps
  ##
  ci.comb.f <- data.frame(TE = object$Comb.fixed,
                          seTE = object$seComb.fixed,
                          lower = object$lower.Comb.fixed,
                          upper = object$upper.Comb.fixed,
                          statistic = object$statistic.Comb.fixed,
                          p = object$pval.Comb.fixed,
                          stringsAsFactors = FALSE)
  rownames(ci.comb.f) <- object$trts
  ##
  ci.comb.r <- data.frame(TE = object$Comb.random,
                          seTE = object$seComb.random,
                          lower = object$lower.Comb.random,
                          upper = object$upper.Comb.random,
                          statistic = object$statistic.Comb.random,
                          p = object$pval.Comb.random,
                          stringsAsFactors = FALSE)
  rownames(ci.comb.r) <- object$trts
  
  
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
              s = object$s,
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
              lower.I2 = object$lower.I2, upper.I2 = object$upper.I2,
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
              level.ma = object$level.ma,
              fixed = fixed,
              random = random,
              ##
              ci.lab = paste0(round(100 * object$level.ma, 1),"%-CI"),
              ##
              reference.group = object$reference.group,
              baseline.reference = object$baseline.reference,
              all.treatments = object$all.treatments,
              seq = object$seq,
              ##
              tau.preset = object$tau.preset,
              ##
              sep.comps = object$sep.comps,
              nchar.comps = replaceNULL(object$nchar.comps, 666),
              ##
              inactive = object$inactive,
              sep.comps = object$sep.comps,
              ##
              backtransf = object$backtransf,
              ##
              title = object$title,
              ##
              x = object,
              ##
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- c(if (inherits(object, "discomb")) "summary.discomb",
                  "summary.netcomb")
  
  res
}
