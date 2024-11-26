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
#' @param overall.hetstat A logical indicating whether to print heterogeneity
#'   measures.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots.
#' @param nchar.comps A numeric defining the minimum number of
#'   characters used to create unique component names.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @return
#'
#' A list of class "summary.netcomb" is returned with the following elements:
#' 
#' \item{k}{Total number of studies.}
#' \item{m}{Total number of pairwise comparisons.}
#' \item{n}{Total number of treatments.}
#' \item{d}{Total number of designs (corresponding to the unique set
#'   of treatments compared within studies).}
#' \item{c}{Total number of components.}
#' \item{s}{Total number of subnetworks (for \code{\link{discomb}}).}
#' 
#' \item{trts}{Treatments included in network meta-analysis.}
#' \item{k.trts}{Number of studies evaluating a treatment.}
#' \item{n.trts}{Number of observations receiving a treatment (if
#'   available).}
#' \item{events.trts}{Number of events observed for a treatment (if
#'   available).}
#' 
#' \item{studies}{Study labels coerced into a factor with its levels
#'   sorted alphabetically.}
#' \item{narms}{Number of arms for each study.}
#' 
#' \item{designs}{Vector with unique designs present in the network. A
#'   design corresponds to the set of treatments compared within a
#'   study.}
#' 
#' \item{comps}{Components included in network meta-analysis.}
#' \item{k.comps}{Number of studies evaluating a component.}
#' \item{n.comps}{Number of observations receiving a component (if
#'   available).}
#' \item{events.comps}{Number of events observed for a component (if
#'   available).}
#' 
#' \item{comparison}{Results for pairwise comparisons (data frame with
#'   columns studlab, treat1, treat2, TE, seTE, lower, upper, z, p).}
#' \item{comparison.nma.common}{Results for pairwise comparisons based
#'   on common effects NMA model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p, leverage).}
#' \item{comparison.nma.random}{Results for pairwise comparisons based
#'   on random effects NMA model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p).}
#' \item{comparison.cnma.common}{Results for pairwise comparisons based
#'   on common effects CNMA model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p, leverage).}
#' \item{comparison.cnma.random}{Results for pairwise comparisons based
#'   on random effects CNMA model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p).}
#' 
#' \item{components.common}{Results for components based
#'   on common effects CNMA model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p, leverage).}
#' \item{components.random}{Results for components based
#'   on random effects CNMA model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p).}
#' 
#' \item{combinations.common}{Results for available combinations based
#'   on common effects CNMA model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p, leverage).}
#' \item{combinations.random}{Results for available combinations based
#'   on random effects CNMA model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p).}
#' 
#' \item{common}{Results for common effects model (a list with elements
#'   TE, seTE, lower, upper, z, p).}
#' \item{random}{Results for random effects model (a list with
#'   elements TE, seTE, lower, upper, z, p).}
#' \item{predict}{Prediction intervals (a list with elements seTE,
#'   lower, upper).}
#' 
#' \item{Q.additive}{Overall heterogeneity / inconsistency statistic.}
#' \item{df.Q.additive}{Degrees of freedom for test of heterogeneity /
#'   inconsistency.}
#' \item{pval.Q.additive}{P-value for test of heterogeneity / inconsistency.}
#' \item{I2, lower.I2, upper.I2}{I-squared, lower and upper confidence
#'   limits.}
#' \item{tau}{Square-root of between-study variance.}
#' 
#' \item{Q.additive}{Overall heterogeneity / inconsistency statistic (CNMA model).}
#' \item{df.Q.additive}{Degrees of freedom for test of heterogeneity /
#'   inconsistency (CNMA model).}
#' \item{pval.Q.additive}{P-value for test of heterogeneity / inconsistency
#'   (CNMA model).}
#' 
#' \item{Q.standard}{Overall heterogeneity statistic (NMA model).}
#' \item{df.Q.heterogeneity}{Degrees of freedom for test of overall
#'   heterogeneity (NMA model).}
#' \item{pval.Q.heterogeneity}{P-value for test of overall
#'   heterogeneity (NMA model).}
#' 
#' \item{Q.diff}{Q statistic for difference between CNMA and NMA model.}
#' \item{df.Q.diff, pval.Q.diff}{Corresponding degrees of freedom and p-value.}
#'
#' \item{sm}{A character string indicating underlying summary
#'   measure.}
#' \item{method}{A character string indicating which method is to be
#'   used for pooling of studies.}
#' \item{level}{The level used to calculate confidence intervals for
#'   individual studies.}
#' \item{level.ma}{The level used to calculate confidence intervals
#'   for pooled estimates.}
#' 
#' \item{ci.lab}{Label for confidence interval.}
#' 
#' \item{reference.group, baseline.reference}{As defined above.}
#' \item{all.treatments}{As defined above.}
#' \item{seq}{A character specifying the sequence of treatments.}
#' 
#' \item{tau.preset}{An optional value for the square-root of the
#'   between-study variance \eqn{\tau^2}.}
#' 
#' \item{sep.comps}{A character used in comparison names as separator
#'   between component labels.}
#' \item{nchar.comps}{A numeric defining the minimum number of
#'   characters used to create unique component names.}
#' 
#' \item{overall.hetstat, backtransf}{As defined above.}
#'  
#' \item{title}{Title of meta-analysis / systematic review.}
#' \item{x}{\code{\link{netmeta}} object (if available).}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{netmeta} used to create object.}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netcomb}}, \code{\link{discomb}}
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
                            overall.hetstat = object$overall.hetstat,
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
  chklogical(overall.hetstat)
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
              I2 = object$I2,
              lower.I2 = object$lower.I2, upper.I2 = object$upper.I2,
              tau = object$tau,
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
              #
              incr = object$incr,
              method.incr = object$method.incr,
              allstudies = object$allstudies,
              cc.pooled = object$cc.pooled,
              #
              overall.hetstat = overall.hetstat,
              backtransf = backtransf,
              #
              title = object$title,
              ##
              x = object,
              ##
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  #
  res$x$common <- common
  res$x$random <- random
  res$x$overall.hetstat <- overall.hetstat
  res$x$backtransf <- backtransf
  res$x$nchar.comps <- nchar.comps
  #
  class(res) <- c(if (inherits(object, "discomb")) "summary.discomb",
                  "summary.netcomb")
  
  res
}
