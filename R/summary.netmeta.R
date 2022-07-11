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
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @return
#'
#' A list is returned with the following elements:
#' \item{comparison}{Results for pairwise comparisons (data frame with
#'   columns studlab, treat1, treat2, TE, seTE, lower, upper, z, p).}
#' \item{comparison.nma.common}{Results for pairwise comparisons based
#'   on common effects model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p, leverage).}
#' \item{comparison.nma.random}{Results for pairwise comparisons based
#'   on random effects model (data frame with columns studlab, treat1,
#'   treat2, TE, seTE, lower, upper, z, p).}
#' \item{common}{Results for common effects model (a list with elements
#'   TE, seTE, lower, upper, z, p).}
#' \item{random}{Results for random effects model (a list with
#'   elements TE, seTE, lower, upper, z, p).}
#' \item{predict}{Prediction intervals (a list with elements seTE,
#'   lower, upper).}
#' \item{studies}{Study labels coerced into a factor with its levels
#'   sorted alphabetically.}
#' \item{narms}{Number of arms for each study.}
#' \item{k}{Total number of studies.}
#' \item{m}{Total number of pairwise comparisons.}
#' \item{n}{Total number of treatments.}
#' \item{d}{Total number of designs (corresponding to the unique set
#'   of treatments compared within studies).}
#' \item{Q}{Overall heterogeneity / inconsistency statistic.}
#' \item{df.Q}{Degrees of freedom for test of heterogeneity /
#'   inconsistency.}
#' \item{pval.Q}{P-value for test of heterogeneity / inconsistency.}
#' \item{I2, lower.I2, upper.I2}{I-squared, lower and upper confidence
#'   limits.}
#' \item{tau}{Square-root of between-study variance.}
#' \item{Q.heterogeneity}{Overall heterogeneity statistic.}
#' \item{df.Q.heterogeneity}{Degrees of freedom for test of overall
#'   heterogeneity.}
#' \item{pval.Q.heterogeneity}{P-value for test of overall
#'   heterogeneity.}
#' \item{Q.inconsistency}{Overall inconsistency statistic.}
#' \item{df.Q.inconsistency}{Degrees of freedom for test of overall
#'   inconsistency.}
#' \item{pval.Q.inconsistency}{P-value for test of overall
#'   inconsistency.}
#' \item{sm}{A character string indicating underlying summary
#'   measure.}
#' \item{method}{A character string indicating which method is to be
#'   used for pooling of studies.}
#' \item{level}{The level used to calculate confidence intervals for
#'   individual studies.}
#' \item{level.ma}{The level used to calculate confidence intervals
#'   for pooled estimates.}
#' \item{common, random}{As defined above.}
#' \item{prediction, level.predict}{As defined above.}
#' \item{reference.group, baseline.reference}{As defined above.}
#' \item{all.treatments, backtransf}{As defined above.}
#' \item{ci.lab}{Label for confidence interval.}
#' \item{seq}{A character specifying the sequence of treatments.}
#' \item{tau.preset}{An optional value for the square-root of the
#'   between-study variance \eqn{\tau^2}.}
#' \item{sep.trts}{A character used in comparison names as separator
#'   between treatment labels.}
#' \item{nchar.trts}{A numeric defining the minimum number of
#'   characters used to create unique treatment names.}
#' \item{title}{Title of meta-analysis / systematic review.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}
#' 
#' @keywords summary
#' 
#' @examples
#' data(Senn2013)
#' 
#' # Conduct common effects network meta-analysis
#' #
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD", random = FALSE)
#' print(net1, ref = "plac", digits = 3)
#' summary(net1)
#'
#' \dontrun{
#' # Conduct random effects network meta-analysis
#' #
#' net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD", common = FALSE)
#' print(net2, ref = "plac", digits = 3)
#' summary(net2)
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
  res <- list(comparison = ci.comp,
              comparison.nma.common = ci.nma.common,
              comparison.nma.random = ci.nma.random,
              common = ci.c,
              random = ci.r,
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
              lower.I2 = object$lower.I2, upper.I2 = object$upper.I2,
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
              level.ma = object$level.ma,
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
              ci.lab = paste0(round(100 * object$level.ma, 1),"%-CI"),
              ##
              reference.group = NA,
              baseline.reference = NA,
              all.treatments = NA,
              seq = object$seq,
              ##
              tau.preset = object$tau.preset,
              ##
              n.trts = object$n.trts,
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
  res$x <- object
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
