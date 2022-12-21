#' Summary method for objects of class netcomb
#' 
#' @description
#' Summary method for objects of class \code{netcomb}.
#' 
#' @param object An object of class \code{netcomb}.
#' @param common A logical indicating whether results for the common
#'   effects model should be printed.
#' @param random A logical indicating whether results for the random
#'   effects model should be printed.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots.
#' @param nchar.comps A numeric defining the minimum number of
#'   characters used to create unique component names.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#'
#' @return
#' A list is returned with the same elements as a
#' \code{\link{netcomb}} object.
#'
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
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
#'   sm = "OR", common = FALSE)
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
#'   sm = "OR", common = FALSE)
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
                            common = object$common,
                            random = object$random,
                            backtransf = object$backtransf,
                            nchar.comps = object$nchar.comps,
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
  chklogical(backtransf)
  nchar.comps <- replaceNULL(nchar.comps, 666)
  chknumeric(nchar.comps, min = 1, length = 1)
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
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
  ci.nma.common <- data.frame(studlab = object$studlab,
                              treat1 = object$treat1,
                              treat2 = object$treat2,
                              TE = object$TE.nma.common,
                              seTE = object$seTE.nma.common,
                              lower = object$lower.nma.common,
                              upper = object$upper.nma.common,
                              statistic = object$statistic.nma.common,
                              p = object$pval.nma.common,
                              stringsAsFactors = FALSE)
  ##
  ci.cnma.common <- data.frame(studlab = object$studlab,
                               treat1 = object$treat1,
                               treat2 = object$treat2,
                               TE = object$TE.cnma.common,
                               seTE = object$seTE.cnma.common,
                               lower = object$lower.cnma.common,
                               upper = object$upper.cnma.common,
                               statistic = object$statistic.cnma.common,
                               p = object$pval.cnma.common,
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
  ci.f <- list(TE = object$TE.common,
               seTE = object$seTE.common,
               lower = object$lower.common,
               upper = object$upper.common,
               statistic = object$statistic.common,
               p = object$pval.common)
  ##
  ci.r <- list(TE = object$TE.random,
               seTE = object$seTE.random,
               lower = object$lower.random,
               upper = object$upper.random,
               statistic = object$statistic.random,
               p = object$pval.random)
  ##
  ci.comp.f <- data.frame(TE = object$Comp.common,
                          seTE = object$seComp.common,
                          lower = object$lower.Comp.common,
                          upper = object$upper.Comp.common,
                          statistic = object$statistic.Comp.common,
                          p = object$pval.Comp.common,
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
  ci.comb.f <- data.frame(TE = object$Comb.common,
                          seTE = object$seComb.common,
                          lower = object$lower.Comb.common,
                          upper = object$upper.Comb.common,
                          statistic = object$statistic.Comb.common,
                          p = object$pval.Comb.common,
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
              comparison.nma.common = ci.nma.common,
              comparison.nma.random = ci.nma.random,
              comparison.cnma.common = ci.cnma.common,
              comparison.cnma.random = ci.cnma.random,
              ##
              components.common = ci.comp.f,
              components.random = ci.comp.r,
              ##
              combinations.common = ci.comb.f,
              combinations.random = ci.comb.r,
              ##
              common = ci.f, random = ci.r,
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
              nchar.comps = nchar.comps,
              ##
              inactive = object$inactive,
              sep.comps = object$sep.comps,
              ##
              backtransf = backtransf,
              ##
              title = object$title,
              ##
              x = object,
              ##
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  res$x$common <- common
  res$x$random <- random
  ##
  class(res) <- c(if (inherits(object, "discomb")) "summary.discomb",
                  "summary.netcomb")
  
  res
}
