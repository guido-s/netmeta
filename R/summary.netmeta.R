#' Summary method for objects of class netmeta
#' 
#' @description
#' Summary method for objects of class \code{netmeta}.
#' 
#' @param object An object of class \code{netmeta}.
#' @param common A logical indicating whether results for the common
#'   effects model should be printed.
#' @param random A logical indicating whether results for the random
#'   effects model should be printed.
#' @param prediction A logical indicating whether prediction intervals
#'   should be printed.
#' @param reference.group Reference treatment.
#' @param baseline.reference A logical indicating whether results
#'   should be expressed as comparisons of other treatments versus the
#'   reference treatment (default) or vice versa. This argument is
#'   only considered if \code{reference.group} has been specified.
#' @param all.treatments A logical or \code{"NULL"}. If \code{TRUE},
#'   matrices with all treatment effects, and confidence limits will
#'   be printed.
#' @param overall.hetstat A logical indicating whether to print heterogeneity
#'   measures.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @return
#'
#' A list of class "summary.netmeta" is returned with the following elements:
#' 
#' \item{k}{Total number of studies.}
#' \item{m}{Total number of pairwise comparisons.}
#' \item{n}{Total number of treatments.}
#' \item{d}{Total number of designs (corresponding to the unique set
#'   of treatments compared within studies).}
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
#' \item{comparisons}{Vector with unique direct comparisons present in the
#'   network.}
#' 
#' \item{comparison}{Results for pairwise comparisons (data frame with
#'   columns studlab, treat1, treat2, TE, seTE, lower, upper, z, p).}
#' \item{comparison.nma.common}{Results for pairwise comparisons based
#'   on common effects model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p, leverage).}
#' \item{comparison.nma.random}{Results for pairwise comparisons based
#'   on random effects model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p).}
#' 
#' \item{common}{Results for common effects model (a list with elements
#'   TE, seTE, lower, upper, z, p).}
#' \item{random}{Results for random effects model (a list with
#'   elements TE, seTE, lower, upper, z, p).}
#' \item{predict}{Prediction intervals (a list with elements seTE,
#'   lower, upper).}
#' 
#' \item{Q}{Overall heterogeneity / inconsistency statistic.}
#' \item{df.Q}{Degrees of freedom for test of heterogeneity /
#'   inconsistency.}
#' \item{pval.Q}{P-value for test of heterogeneity / inconsistency.}
#' \item{I2, lower.I2, upper.I2}{I-squared, lower and upper confidence
#'   limits.}
#' \item{tau}{Square-root of between-study variance.}
#' 
#' \item{Q.heterogeneity}{Overall heterogeneity statistic.}
#' \item{df.Q.heterogeneity}{Degrees of freedom for test of overall
#'   heterogeneity.}
#' \item{pval.Q.heterogeneity}{P-value for test of overall
#'   heterogeneity.}
#' 
#' \item{Q.inconsistency}{Overall inconsistency statistic.}
#' \item{df.Q.inconsistency}{Degrees of freedom for test of overall
#'   inconsistency.}
#' \item{pval.Q.inconsistency}{P-value for test of overall
#'   inconsistency.}
#'
#' \item{Q.decomp}{Data frame with columns 'treat1', 'treat2', 'Q',
#'   'df' and 'pval.Q', providing heterogeneity statistics for each
#'   pairwise meta-analysis of direct comparisons.}
#' 
#' \item{sm}{A character string indicating underlying summary
#'   measure.}
#' \item{method}{A character string indicating which method is to be
#'   used for pooling of studies.}
#' \item{level}{The level used to calculate confidence intervals for
#'   individual studies.}
#' \item{level.ma}{The level used to calculate confidence intervals
#'   for pooled estimates.}
#' \item{level.predict}{The level used to calculate prediction
#'   intervals for a new study.}
#' 
#' \item{ci.lab}{Label for confidence interval.}
#' 
#' \item{incr}{Numerical value added to cell frequencies
#'   (if applicable).}
#' \item{method.incr}{A character string indicating which continuity
#'   correction method was used (if applicable).}
#' \item{allstudies}{A logical indicating whether studies with zero
#'   events or non-events in all treatment arms should be included in
#'   an inverse variance meta-analysis (if applicable).}
#' \item{cc.pooled}{A logical indicating whether \code{incr} should be
#'   used as a continuity correction (if applicable).}
#' 
#' \item{reference.group, baseline.reference}{As defined above.}
#' \item{all.treatments}{As defined above.}
#' \item{seq}{A character specifying the sequence of treatments.}
#' 
#' \item{tau.preset}{An optional value for the square-root of the
#'   between-study variance \eqn{\tau^2}.}
#' 
#' \item{sep.trts}{A character used in comparison names as separator
#'   between treatment labels.}
#' \item{nchar.trts}{A numeric defining the minimum number of
#'   characters used to create unique treatment names.}
#' 
#' \item{prediction, overall.hetstat, backtransf}{As defined above.}
#'  
#' \item{title}{Title of meta-analysis / systematic review.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{netmeta} used to create object.}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}
#' 
#' @keywords summary
#' 
#' @examples
#' data(smokingcessation)
#' 
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' net1 <- netmeta(p1)
#'
#' summary(net1)
#' 
#' \dontrun{
#' data(Senn2013)
#' 
#' # Conduct common effects network meta-analysis
#' #
#' net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD", random = FALSE)
#' print(net2, ref = "plac", digits = 3)
#' summary(net2)
#'
#' # Conduct random effects network meta-analysis
#' #
#' net3 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD", common = FALSE)
#' print(net3, ref = "plac", digits = 3)
#' summary(net3)
#' }
#' 
#' @method summary netmeta
#' @export


summary.netmeta <- function(object,
                            common = object$common,
                            random = object$random,
                            prediction = object$prediction,
                            reference.group = object$reference.group,
                            baseline.reference = object$baseline.reference,
                            all.treatments = object$all.treatments,
                            overall.hetstat = object$overall.hetstat,
                            backtransf = object$backtransf,
                            nchar.trts = object$nchar.trts,
                            warn.deprecated = gs("warn.deprecated"),
                            ...) {
  
  ##
  ##
  ## (1) Check for netmeta object and upgrade object
  ##
  ##
  chkclass(object, "netmeta")
  object <- updateversion(object)
  ##
  is.bin <- inherits(object, "netmetabin")
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chklogical(prediction)
  chklogical(baseline.reference)
  chklogical(overall.hetstat)
  chklogical(backtransf)
  chknumeric(nchar.trts, min = 1, length = 1)
  ##
  fun <- "summary.netmeta"
  ##
  ## Check for deprecated arguments in '...'
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
  ci.nma.common <-
    data.frame(studlab = object$studlab,
               treat1 = object$treat1,
               treat2 = object$treat2,
               TE = if (!is.bin) object$TE.nma.common else NA,
               seTE = if (!is.bin) object$seTE.nma.common else NA,
               lower = if (!is.bin) object$lower.nma.common else NA,
               upper = if (!is.bin) object$upper.nma.common else NA,
               statistic = if (!is.bin) object$statistic.nma.common else NA,
               p = if (!is.bin) object$pval.nma.common else NA,
               leverage = if (!is.bin) object$leverage.common else NA,
               stringsAsFactors = FALSE)
  ##
  ci.nma.random <-
    data.frame(studlab = object$studlab,
               treat1 = object$treat1,
               treat2 = object$treat2,
               TE = if (!is.bin) object$TE.nma.random else NA,
               seTE = if (!is.bin) object$seTE.nma.random else NA,
               lower = if (!is.bin) object$lower.nma.random else NA,
               upper = if (!is.bin) object$upper.nma.random else NA,
               statistic = if (!is.bin) object$statistic.nma.random else NA,
               p = if (!is.bin) object$pval.nma.random else NA,
               stringsAsFactors = FALSE)
  ##
  ci.c <- list(TE = object$TE.common,
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
  ci.p <- list(seTE = object$seTE.predict,
               lower = object$lower.predict,
               upper = object$upper.predict)
  
  
  ##
  ##
  ## (4) Create summary.netmeta object
  ##
  ##
  object$common <- common
  object$random <- random
  ##
  res <- list(k = object$k,
              m = object$m,
              n = object$n,
              d = object$d,
              #
              trts = object$trts,
              k.trts = object$k.trts,
              n.trts = object$n.trts,
              events.trts = object$events.trts,
              #
              studies = object$studies,
              narms = object$narms,
              #
              designs = object$designs,
              comparisons = object$comparisons,
              #
              comparison = ci.comp,
              comparison.nma.common = ci.nma.common,
              comparison.nma.random = ci.nma.random,
              #
              common = ci.c,
              random = ci.r,
              predict = ci.p,
              #
              Q = object$Q,
              df.Q = object$df.Q,
              pval.Q = object$pval.Q,
              I2 = object$I2,
              lower.I2 = object$lower.I2, upper.I2 = object$upper.I2,
              tau = object$tau,
              #
              Q.heterogeneity = object$Q.heterogeneity,
              df.Q.heterogeneity = object$df.Q.heterogeneity,
              pval.Q.heterogeneity = object$pval.Q.heterogeneity,
              #
              Q.inconsistency = object$Q.inconsistency,
              df.Q.inconsistency = object$df.Q.inconsistency,
              pval.Q.inconsistency = object$pval.Q.inconsistency,
              #
              Q.decomp = object$Q.decomp,
              #
              sm = object$sm,
              method = object$method,
              level = object$level,
              level.ma = object$level.ma,
              level.predict = object$level.predict,
              #
              ci.lab = paste0(round(100 * object$level.ma, 1),"%-CI"),
              #
              incr = object$incr,
              method.incr = object$method.incr,
              allstudies = object$allstudies,
              cc.pooled = object$cc.pooled,
              #
              reference.group = NA,
              baseline.reference = NA,
              all.treatments = NA,
              seq = object$seq,
              #
              tau.preset = object$tau.preset,
              #
              sep.trts = object$sep.trts,
              nchar.trts = nchar.trts,
              #
              prediction = prediction,
              overall.hetstat = overall.hetstat,
              backtransf = backtransf,
              #
              title = object$title,
              #
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
  res$x <- object
  #
  res$x$common <- common
  res$x$random <- random
  res$x$prediction <- prediction
  res$x$reference.group <- reference.group
  res$x$baseline.reference <- baseline.reference
  res$x$all.treatments <- all.treatments
  res$x$overall.hetstat <- overall.hetstat
  res$x$backtransf <- backtransf
  res$x$nchar.trts <- nchar.trts
  ##
  ## Backward compatibility
  ##
  res$fixed <- res$common
  ##
  if (is.bin)
    class(res) <- c("summary.netmeta", "summary.netmetabin")
  else
    class(res) <- "summary.netmeta"
  
  res
}
