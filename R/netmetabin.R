#' Network meta-analysis of binary outcome data
#' 
#' @description
#' Provides three models for the network meta-analysis of binary data
#' (Mantel-Haenszel method, based on the non-central hypergeometric
#' distribution, and the inverse variance method).
#' 
#' @param event1 Number of events (first treatment).
#' @param n1 Number of observations (first treatment).
#' @param event2 Number of events (second treatment).
#' @param n2 Number of observations (second treatment)
#' @param treat1 Label/Number for first treatment.
#' @param treat2 Label/Number for second treatment.
#' @param studlab An optional - but important! - vector with study
#'   labels (see Details).
#' @param data An optional data frame containing the study
#'   information.
#' @param subset An optional vector specifying a subset of studies to
#'   be used.
#' @param sm A character string indicating underlying summary measure,
#'   i.e., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"}.
#' @param method A character string indicating which method is to be
#'   used for pooling of studies. One of \code{"Inverse"},
#'   \code{"MH"}, \code{"NCH"}, or \code{"LRP"}, can be abbreviated.
#' @param cc.pooled A logical indicating whether \code{incr} should be
#'   used as a continuity correction, when calculating the network
#'   meta-analysis estimates.
#' @param incr A numerical value which is added to each cell count,
#'   i.e., to the numbers of events and non-events, of all treatment
#'   arms in studies with zero events or non-events in any of the
#'   treatment arms ("continuity correction").
#' @param method.incr A character string indicating which continuity
#'   correction method should be used (\code{"only0"},
#'   \code{"if0all"}, or \code{"all"}), see \code{\link[meta]{metabin}}.
#' @param allstudies A logical indicating whether studies with zero
#'   events or non-events in all treatment arms should be included in
#'   an inverse variance meta-analysis (applies only if \code{method =
#'   "Inverse"} and \code{sm} is equal to either \code{"RR"} or
#'   \code{"OR"}).
#' @param level The level used to calculate confidence intervals for
#'   individual studies.
#' @param level.ma The level used to calculate confidence intervals
#'   for network estimates.
#' @param common A logical indicating whether a common effects network
#'   meta-analysis should be conducted.
#' @param random A logical indicating whether a random effects network
#'   meta-analysis should be conducted.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed (only considered if \code{method =
#'   "Inverse"}).
#' @param level.predict The level used to calculate prediction
#'   interval for a new study (only considered if \code{method =
#'   "Inverse"}).
#' @param reference.group Reference treatment.
#' @param baseline.reference A logical indicating whether results
#'   should be expressed as comparisons of other treatments versus the
#'   reference treatment (default) or vice versa. This argument is
#'   only considered if \code{reference.group} has been specified.
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"desirable"}) or
#'   harmful (\code{"undesirable"}) effect (passed on to
#'   \code{\link{netrank}}, can be abbreviated.
#' @param all.treatments A logical or \code{"NULL"}. If \code{TRUE},
#'   matrices with all treatment effects, and confidence limits will
#'   be printed.
#' @param seq A character or numerical vector specifying the sequence
#'   of treatments in printouts.
#' @param tau.preset An optional value for manually setting the
#'   square-root of the between-study variance \eqn{\tau^2} (only
#'   considered if \code{method = "Inverse"}).
#' @param tol.multiarm A numeric for the tolerance for consistency of
#'   treatment estimates in multi-arm studies which are consistent by
#'   design (only considered if \code{method = "Inverse"}).
#' @param tol.multiarm.se A numeric for the tolerance for consistency
#'   of standard errors in multi-arm studies which are consistent by
#'   design (only considered if the argument is not \code{NULL} and
#'   \code{method = "Inverse"}).
#' @param details.chkmultiarm A logical indicating whether treatment
#'   estimates and / or variances of multi-arm studies with
#'   inconsistent results or negative multi-arm variances should be
#'   printed (only considered if \code{method = "Inverse"}).
#' @param details.chkdata A logical indicating whether number of
#'   events and participants of studies with inconsistent data should
#'   be printed.
#' @param sep.trts A character used in comparison names as separator
#'   between treatment labels.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
#' @param func.inverse R function used to calculate the pseudoinverse
#'   of the Laplacian matrix L (see \code{\link{netmeta}}).
#' @param overall.hetstat A logical indicating whether to print heterogeneity
#'   measures.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param title Title of meta-analysis / systematic review.
#' @param keepdata A logical indicating whether original data(set)
#'   should be kept in netmeta object.
#' @param addincr Deprecated argument (replaced by 'method.incr');
#'   see \code{\link[meta]{metabin}}.
#' @param allincr Deprecated argument (replaced by 'method.incr');
#'   see \code{\link[meta]{metabin}}.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if studies are excluded from meta-analysis due to zero
#'   standard errors).
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @details
#' This function implements four models for the network meta-analysis
#' of binary data:
#' \itemize{
#' \item The Mantel-Haenszel network meta-analysis model, as described
#'   in Efthimiou et al. (2019) (\code{method = "MH"});
#' \item a network meta-analysis model using the non-central
#'   hypergeometric distribution with the Breslow approximation, as
#'   described in Stijnen et al. (2010) (\code{method = "NCH"});
#' \item a logistic regression with penalised likelihood
#'   (\code{method = "LRP"});
#' \item the inverse variance method for network meta-analysis
#'   (\code{method = "Inverse"}), also provided by
#'   \code{\link{netmeta}}.
#' }
#' 
#' Comparisons belonging to multi-arm studies are identified by
#' identical study labels (argument \code{studlab}). It is therefore
#' important to use identical study labels for all comparisons
#' belonging to the same multi-arm study.
#' 
#' Data entry for this function is in \emph{contrast-based} format,
#' that is, each line of the data corresponds to a single pairwise
#' comparison between two treatments (arguments \code{treat1},
#' \code{treat2}, \code{event1}, \code{n1}, \code{event2}, and
#' \code{n2}). If data are provided in \emph{arm-based} format, that
#' is, number of events and participants are given for each treatment
#' arm separately, function \code{\link[meta]{pairwise}} can be used to
#' transform the data to \emph{contrast-based} format (see help page
#' of function \code{\link[meta]{pairwise}}).
#' 
#' Note, all pairwise comparisons must be provided for a multi-arm
#' study. Consider a multi-arm study of \emph{p} treatments with known
#' variances. For this study, the number of events and observations
#' must be provided for each treatment, for each of \emph{p}(\emph{p}
#' - 1) / 2 possible comparisons in separate lines in the data. For
#' instance, a three-arm study contributes three pairwise comparisons,
#' a four-arm study even six pairwise comparisons. Function
#' \code{\link[meta]{pairwise}} automatically calculates all pairwise
#' comparisons for multi-arm studies.
#' 
#' For \code{method = "Inverse"}, both common and random effects
#' models are calculated regardless of values choosen for arguments
#' \code{common} and \code{random}. Accordingly, the network estimates
#' for the random effects model can be extracted from component
#' \code{TE.random} of an object of class \code{"netmeta"} even if
#' argument \code{random = FALSE}.  However, all functions in R
#' package \bold{netmeta} will adequately consider the values for
#' \code{common} and \code{random}. E.g. function
#' \code{\link{print.summary.netmeta}} will not print results for the
#' random effects model if \code{random = FALSE}.
#' 
#' For the random-effects model, the direct treatment estimates are
#' based on the common between-study variance \eqn{\tau^2} from the
#' network meta-analysis.
#'
#' For \code{method = "MH"} and \code{method = "NCH"}, only a common
#' effects model is available.
#' 
#' By default, treatment names are not abbreviated in
#' printouts. However, in order to get more concise printouts,
#' argument \code{nchar.trts} can be used to define the minimum number
#' of characters for abbreviated treatment names (see
#' \code{\link{abbreviate}}, argument \code{minlength}). R function
#' \code{\link{treats}} is utilised internally to create abbreviated
#' treatment names.
#' 
#' Names of treatment comparisons are created by concatenating
#' treatment labels of pairwise comparisons using \code{sep.trts} as
#' separator (see \code{\link{paste}}). These comparison names are
#' used in the covariance matrices \code{Cov.common} and
#' \code{Cov.random} and in some R functions, e.g,
#' \code{\link{decomp.design}}. By default, a colon is used as the
#' separator. If any treatment label contains a colon the following
#' characters are used as separator (in consecutive order):
#' \code{"-"}, \code{"_"}, \code{"/"}, \code{"+"}, \code{"."},
#' \code{"|"}, and \code{"*"}. If all of these characters are used in
#' treatment labels, a corresponding error message is printed asking
#' the user to specify a different separator.
#'
#' @return
#' An object of class \code{netmetabin} and \code{netmeta} with
#' corresponding \code{print}, \code{summary}, \code{forest}, and
#' \code{netrank} functions. The object is a list containing the
#' following components:
#' \item{studlab, treat1, treat2}{As defined above.}
#' \item{n1, n2, event1, event2}{As defined above.}
#' \item{TE}{Estimate of treatment effect, i.e. difference between
#'   first and second treatment (e.g. log odds ratio).}
#' \item{seTE}{Standard error of treatment estimate.}
#' \item{k}{Total number of studies.}
#' \item{m}{Total number of pairwise comparisons.}
#' \item{n}{Total number of treatments.}
#' \item{d}{Total number of designs (corresponding to the unique set
#'   of treatments compared within studies).}
#' \item{trts}{Treatments included in network meta-analysis.}
#' \item{k.trts}{Number of studies evaluating a treatment.}
#' \item{n.trts}{Number of observations receiving a treatment.}
#' \item{events.trts}{Number of events observed for a treatment.}
#' \item{studies}{Study labels coerced into a factor with its levels
#'   sorted alphabetically.}
#' \item{narms}{Number of arms for each study.}
#' \item{designs}{Unique list of designs present in the network. A
#'   design corresponds to the set of treatments compared within a
#'   study.}
#' \item{TE.common, seTE.common}{\emph{n}x\emph{n} matrix with estimated
#'   overall treatment effects and standard errors for common effects
#'   model.}
#' \item{lower.common, upper.common}{\emph{n}x\emph{n} matrices with
#'   lower and upper confidence interval limits for common effects
#'   model.}
#' \item{statistic.common, pval.common}{\emph{n}x\emph{n} matrices with
#'   z-value and p-value for test of overall treatment effect under
#'   common effects model.}
#' \item{TE.random, seTE.random}{\emph{n}x\emph{n} matrix with
#'   estimated overall treatment effects and standard errors for
#'   random effects model (only available if \code{method =
#'   "Inverse"}).}
#' \item{lower.random, upper.random}{\emph{n}x\emph{n} matrices with
#'   lower and upper confidence interval limits for random effects
#'   model (only available if \code{method = "Inverse"}).}
#' \item{statistic.random, pval.random}{\emph{n}x\emph{n} matrices
#'   with z-value and p-value for test of overall treatment effect
#'   under random effects model (only available if \code{method =
#'   "Inverse"}).}
#' \item{TE.direct.common, seTE.direct.common}{\emph{n}x\emph{n} matrix
#'   with estimated treatment effects and standard errors from direct
#'   evidence under common effects model.}
#' \item{lower.direct.common, upper.direct.common}{\emph{n}x\emph{n}
#'   matrices with lower and upper confidence interval limits from
#'   direct evidence under common effects model.}
#' \item{statistic.direct.common, pval.direct.common}{\emph{n}x\emph{n}
#'   matrices with z-value and p-value for test of overall treatment
#'   effect from direct evidence under common effects model.}
#' \item{TE.direct.random, seTE.direct.random}{\emph{n}x\emph{n}
#'   matrix with estimated treatment effects and standard errors from
#'   direct evidence under random effects model (only available if
#'   \code{method = "Inverse"}).}
#' \item{lower.direct.random, upper.direct.random}{\emph{n}x\emph{n}
#'   matrices with lower and upper confidence interval limits from
#'   direct evidence under random effects model (only available if
#'   \code{method = "Inverse"}).}
#' \item{statistic.direct.random,
#'   pval.direct.random}{\emph{n}x\emph{n} matrices with z-value and
#'   p-value for test of overall treatment effect from direct evidence
#'   under random effects model (only available if \code{method =
#'   "Inverse"}).}
#' \item{Q}{Overall heterogeneity / inconsistency statistic. (only
#'   available if \code{method = "Inverse"})}
#' \item{df.Q}{Degrees of freedom for test of heterogeneity /
#'   inconsistency.}
#' \item{pval.Q}{P-value for test of heterogeneity / inconsistency.}
#' \item{I2, lower.I2, upper.I2}{I-squared, lower and upper confidence
#'   limits (only available if \code{method = "Inverse"}).}
#' \item{tau}{Square-root of between-study variance (only available if
#'   \code{method = "Inverse"}).}
#' \item{Q.heterogeneity}{Overall heterogeneity statistic. (only
#'   available if \code{method = "Inverse"})}
#' \item{df.Q.heterogeneity}{Degrees of freedom for test of overall
#'   heterogeneity.}
#' \item{pval.Q.heterogeneity}{P-value for test of overall
#'   heterogeneity.}
#' \item{Q.inconsistency}{Overall inconsistency statistic.}
#' \item{df.Q.inconsistency}{Degrees of freedom for test of overall
#'   inconsistency.}
#' \item{pval.Q.inconsistency}{P-value for test of overall
#'   inconsistency.}
#' \item{A.matrix}{Adjacency matrix (\emph{n}x\emph{n}).}
#' \item{H.matrix.common}{Hat matrix (\emph{m}x\emph{m})}
#' \item{n.matrix}{\emph{n}x\emph{n} matrix with number of
#'   observations in direct comparisons.}
#' \item{events.matrix}{\emph{n}x\emph{n} matrix with number of events
#'   in direct comparisons.}
#' \item{sm, method, level, level.ma}{As defined above.}
#' \item{incr, method.incr, allstudies, cc.pooled}{As defined
#'   above.}
#' \item{addincr, allincr}{As defined above.}
#' \item{common, random}{As defined above.}
#' \item{prediction, level.predict}{As defined above.}
#' \item{reference.group, baseline.reference, small.values, all.treatments}{As
#'   defined above.}
#' \item{seq, tau.preset, tol.multiarm, tol.multiarm.se}{As defined
#'   above.}
#' \item{details.chkmultiarm, details.chkdata}{As defined above.}
#' \item{sep.trts, nchar.trts, overall.hetstat}{As defined above.}
#' \item{backtransf, title, warn, warn.deprecated}{As defined above.}
#' \item{data}{Dataset (in contrast-based format).}
#' \item{data.design}{List with data in arm-based format (each list
#'   element corresponds to a single design).}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Orestis Efthimiou \email{oremiou@@gmail.com}, Guido
#'   Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link[meta]{pairwise}}, \code{\link{netmeta}},
#'   \code{\link[metadat]{dat.gurusamy2011}},
#'   \code{\link[metadat]{dat.dong2013}}
#' 
#' @references
#' Efthimiou O, Rücker G, Schwarzer G, Higgins J, Egger M, Salanti G
#' (2019):
#' A Mantel-Haenszel model for network meta-analysis of rare events.
#' \emph{Statistics in Medicine},
#' \bold{38}, 2992--3012
#' 
#' Senn S, Gavini F, Magrez D, Scheen A (2013):
#' Issues in performing a network meta-analysis.
#' \emph{Statistical Methods in Medical Research},
#' \bold{22}, 169--89
#' 
#' Stijnen T, Hamza TH, Ozdemir P (2010):
#' Random effects meta-analysis of event outcome in the framework of
#' the generalized linear mixed model with applications in sparse
#' data.
#' \emph{Statistics in Medicine},
#' \bold{29}, 3046--67

#' @examples
#' # Only consider first ten studies (to reduce runtime of example)
#' #
#' first10 <- subset(dat.dong2013, id <= 10)
#' 
#' # Transform data from long arm-based format to contrast-based
#' # format. Argument 'sm' has to be used for odds ratio as summary
#' # measure; by default the risk ratio is used in the metabin
#' # function called internally.
#' #
#' p1 <- pairwise(treatment, death, randomized, studlab = id,
#'   data = first10, sm = "OR")
#' 
#' # Conduct Mantel-Haenszel network meta-analysis (without continuity
#' # correction)
#' #
#' nb1 <- netmetabin(p1, ref = "plac")
#' nb1
#' 
#' # Obtain the league table
#' #
#' netleague(nb1)
#' 
#' \dontrun{
#' # Conduct Mantel-Haenszel network meta-analysis for the whole
#' # dataset
#' #
#' p2 <- pairwise(treatment, death, randomized, studlab = id,
#'   data = dat.dong2013, sm = "OR")
#' netmetabin(p2, ref = "plac")
#'   
#' # Conduct network meta-analysis using the non-central
#' # hypergeometric model (without continuity correction)
#' #
#' netmetabin(p2, ref = "plac", method = "NCH")
#' 
#' # Conduct Mantel-Haenszel network meta-analysis (with continuity
#' # correction of 0.5; include all studies)
#' #
#' netmetabin(p2, ref = "plac", cc.pooled = TRUE)
#' 
#' p3 <- pairwise(treatment, death, n, studlab = study,
#'   data = dat.gurusamy2011, sm = "OR")
#' 
#' # Conduct Mantel-Haenszel network meta-analysis (without continuity
#' # correction)
#' #
#' netmetabin(p3, ref = "cont")
#' }
#' 
#' @export netmetabin

netmetabin <- function(event1, n1, event2, n2,
                       treat1, treat2, studlab,
                       data = NULL, subset = NULL,
                       sm,
                       method = "MH",
                       cc.pooled = FALSE,
                       #
                       incr, method.incr, allstudies,
                       #
                       level = gs("level"),
                       level.ma = gs("level.ma"),
                       common = gs("common"),
                       random = method %in% c("Inverse", "LRP") &
                         (gs("random") | !is.null(tau.preset)),
                       ##
                       prediction = gs("prediction"),
                       level.predict = gs("level.predict"),
                       ##
                       reference.group = "",
                       baseline.reference = gs("baseline.reference"),
                       small.values = gs("small.values"),
                       all.treatments = gs("all.treatments"),
                       seq = gs("seq"),
                       #
                       tau.preset = NULL,
                       ##
                       tol.multiarm = 0.001,
                       tol.multiarm.se = NULL,
                       details.chkmultiarm = FALSE,
                       details.chkdata = TRUE,
                       ##
                       sep.trts = ":",
                       nchar.trts = 666,
                       ##
                       func.inverse = invmat,
                       ##
                       overall.hetstat = gs("overall.hetstat"),
                       backtransf = gs("backtransf"),
                       ##
                       title = gs("title"),
                       keepdata = gs("keepdata"),
                       #
                       addincr, allincr,
                       #
                       warn = gs("warn"),
                       warn.deprecated = gs("warn.deprecated"),
                       ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  
  missing.common <- missing(common)
  missing.random <- missing(random)
  #
  missing.level.ma <- missing(level.ma)
  missing.sep.trts <- missing(sep.trts)
  #
  missing.incr <- missing(incr)
  missing.method.incr <- missing(method.incr)
  missing.allstudies <- missing(allstudies)
  #
  missing.addincr <- missing(addincr)
  missing.allincr <- missing(allincr)
  #
  modtext <-
    paste0("must be equal to 'Inverse' (classic network meta-analysis), ",
           "'MH' (Mantel-Haenszel, the default), ",
           "'NCH' (common-effects non-central hypergeometric), or",
           "'LRP' (penalized logistic regression).")
  method <- setchar(method, c("Inverse", "MH", "NCH", "LRP"), modtext)
  is.mh.nch <- !(method %in% c("Inverse", "LRP"))
  is.lrp <- method == "LRP"
  #
  chklogical(cc.pooled)
  ##
  chklevel(level)
  chklevel(level.predict)
  ##
  chklogical(prediction)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  level.ma <- deprecated(level.ma, missing.level.ma, args, "level.comb",
                         warn.deprecated)
  chklevel(level.ma)
  ##
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  random <-
    deprecated(random, missing.random, args, "comb.random", warn.deprecated)
  chklogical(random)
  ##
  if (is.mh.nch & !common) {
    warning("Argument 'common' set to TRUE for Mantel-Haenszel ",
            "method and non-central hypergeometric distribution.",
            call. = FALSE)
    common <- TRUE
  }
  ##
  if (is.mh.nch & random) {
    warning("Argument 'random' set to FALSE for Mantel-Haenszel ",
            "method and non-central hypergeometric distribution.",
            call. = FALSE)
    random <- FALSE
  }
  ##
  if (is.mh.nch & prediction) {
    warning("Argument 'prediction' set to FALSE for Mantel-Haenszel ",
            "method and non-central hypergeometric distribution.",
            call. = FALSE)
    prediction <- FALSE
  }
  #
  chkchar(reference.group, length = 1)
  chklogical(baseline.reference)
  small.values <- setsv(small.values)
  #
  if (!is.null(all.treatments))
    chklogical(all.treatments)
  else {
    if (reference.group == "")
      all.treatments <- TRUE
    else
      all.treatments <- FALSE
  }
  #
  if (!is.null(tau.preset))
    chknumeric(tau.preset, min = 0, length = 1)
  ##
  chknumeric(tol.multiarm, min = 0, length = 1)
  if (!is.null(tol.multiarm.se))
    chknumeric(tol.multiarm.se, min = 0, length = 1)
  chklogical(details.chkmultiarm)
  chklogical(details.chkdata)
  ##
  chkchar(sep.trts)
  chknumeric(nchar.trts, min = 1, length = 1)
  #
  overall.hetstat <- replaceNULL(overall.hetstat, TRUE)
  chklogical(overall.hetstat)
  chklogical(backtransf)
  #
  chkchar(title)
  chklogical(keepdata)
  chklogical(warn)
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  if (nulldata)
    data <- sfsp
  ##
  ## Catch 'event1', 'event2', 'n1', 'n2', 'treat1', 'treat2', and
  ## 'studlab' from data:
  ##
  event1 <- catch("event1", mc, data, sfsp)
  ##
  if (is.data.frame(event1)  &&
      (!is.null(attr(event1, "pairwise")) ||
       inherits(event1, "pairwise"))) {
    is.pairwise <- TRUE
    ##
    if (missing.incr)
      incr <- replaceNULL(attr(event1, "incr"), gs("incr"))
    if (missing.method.incr)
      method.incr <- replaceNULL(attr(event1, "method.incr"), gs("method.incr"))
    if (missing(allstudies))
      allstudies <- replaceNULL(attr(event1, "allstudies"), gs("allstudies"))
    #
    missing.incr <- FALSE
    missing.method.incr <- FALSE
    missing.allstudies <- FALSE
    #
    if (missing(sm) & method == "Inverse")
      sm <- replaceNULL(attr(event1, "sm"), gs("smbin"))
    else if (method != "Inverse") {
      if (!missing(sm) && tolower(sm) != "or")
        warning("Argument 'sm' set to 'OR'.", call. = FALSE)
      sm <- "OR"
    }
    ##
    n1 <- event1$n1
    event2 <- event1$event2
    n2 <- event1$n2
    treat1 <- event1$treat1
    treat2 <- event1$treat2
    studlab <- event1$studlab
    ##
    pairdata <- event1
    data <- event1
    ##
    event1 <- event1$event1
    ##
    chknull(event1, text = "Variable")
    chknull(n1, text = "Variable")
    chknull(event2, text = "Variable")
    chknull(n2, text = "Variable")
    chknull(treat1, text = "Variable")
    chknull(treat2, text = "Variable")
    chknull(studlab, text = "Variable")
  }
  else {
    is.pairwise <- FALSE
    ##
    if (missing.incr)
      incr <- gs("incr")
    if (missing(method.incr))
      method.incr <- gs("method.incr")
    if (missing(allstudies))
      allstudies <- gs("allstudies")
    #
    if (missing(sm) & method == "Inverse") {
      if (!nulldata && !is.null(attr(data, "sm")))
        sm <- attr(data, "sm")
      else
        sm <- ""
    }
    else if (method != "Inverse") {
      if (!missing(sm) && tolower(sm) != "or")
        warning("Argument 'sm' set to 'OR'.", call. = FALSE)
      sm <- "OR"
    }
    ##
    event2 <- catch("event2", mc, data, sfsp)
    ##
    n1 <- catch("n1", mc, data, sfsp)
    n2 <- catch("n2", mc, data, sfsp)
    ##
    treat1 <- catch("treat1", mc, data, sfsp)
    treat2 <- catch("treat2", mc, data, sfsp)
    ##
    studlab <- catch("studlab", mc, data, sfsp)
    #
    if (!missing.incr)
      incr <- catch("incr", mc, data, sfsp)
  }
  #
  addincr <-
    deprecated2(method.incr, missing.method.incr, addincr, missing.addincr,
                warn.deprecated)
  allincr <-
    deprecated2(method.incr, missing.method.incr, allincr, missing.allincr,
                warn.deprecated)
  #
  if (missing.method.incr) {
    method.incr <- gs("method.incr")
    ##
    if (is.logical(addincr) && addincr)
      method.incr <- "all"
    else if (is.logical(allincr) && allincr)
      method.incr <- "if0all"
  }
  #
  addincr <- allincr <- FALSE
  if (!(sm == "ASD" | method %in% c("Peto", "GLMM"))) {
    if (method.incr == "all")
      addincr <- TRUE
    else if (method.incr == "if0all")
      allincr <- TRUE
  }
  #
  chknull(event1)
  chknull(n1)
  chknull(event2)
  chknull(n2)
  chknull(treat1)
  chknull(treat2)
  chknull(studlab)
  #
  chknumeric(event1)
  chknumeric(n1)
  chknumeric(event2)
  chknumeric(n2)
  ##
  k.Comp <- length(event1)
  ##
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  if (is.factor(studlab))
    studlab <- as.character(studlab)
  ##
  ## Remove leading and trailing whitespace
  ##
  treat1 <- rmSpace(rmSpace(treat1, end = TRUE))
  treat2 <- rmSpace(rmSpace(treat2, end = TRUE))
  ##
  ## Keep original order of studies
  ##
  .order <- seq_along(studlab)
  ##
  subset <- catch("subset", mc, data, sfsp)
  missing.subset <- is.null(subset)
  
  
  ##
  ##
  ## (2b) Store complete dataset in list object data
  ##
  ##
  if (nulldata & !is.pairwise)
    data <- data.frame(.event1 = event1)
  else if (nulldata & is.pairwise) {
    data <- pairdata
    data$.order <- .order
    data$.event1 <- event1
  }
  else
    data$.event1 <- event1
  ##
  data$.n1 <- n1
  data$.event2 <- event2
  data$.n2 <- n2
  data$.treat1 <- treat1
  data$.treat2 <- treat2
  data$.studlab <- studlab
  data$.order <- .order
  ##
  if (!missing.subset) {
    if (length(subset) == dim(data)[1])
      data$.subset <- subset
    else {
      data$.subset <- FALSE
      data$.subset[subset] <- TRUE
    }
  }
  ##
  chklogical(incr)
  chklogical(allstudies)
  ##
  m.data <- metabin(event1, n1, event2, n2,
                    sm = sm, method = "Inverse",
                    incr = incr, method.incr = method.incr,
                    allstudies = allstudies,
                    method.tau = "DL", method.tau.ci = "",
                    warn = FALSE,
                    warn.deprecated = FALSE)
  ##
  data$.TE <- m.data$TE
  data$.seTE <- m.data$seTE
  
  
  ##
  ##
  ## (3) Use subset for analysis
  ##
  ##
  if (!missing.subset) {
    if ((is.logical(subset) & (sum(subset) > k.Comp)) ||
        (length(subset) > k.Comp))
      stop("Length of subset is larger than number of studies.",
           call. = FALSE)
    ##
    studlab <- studlab[subset]
    treat1 <- treat1[subset]
    treat2 <- treat2[subset]
    ##
    event1 <- event1[subset]
    event2 <- event2[subset]
    n1 <- n1[subset]
    n2 <- n2[subset]
    ##
    .order <- .order[subset]
  }
  ##
  ## Check for correct treatment order within comparison
  ##
  wo <- treat1 > treat2
  ##
  if (any(wo)) {
    tevent1 <- event1
    event1[wo] <- event2[wo]
    event2[wo] <- tevent1[wo]
    ##
    tn1 <- n1
    n1[wo] <- n2[wo]
    n2[wo] <- tn1[wo]
    ##
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
  }
  ##
  labels <- sort(unique(c(treat1, treat2)))
  ##
  if (compmatch(labels, sep.trts)) {
    if (!missing.sep.trts)
      warning("Separator '", sep.trts, "' used in at least ",
              "one treatment label. Try to use predefined separators: ",
              "':', '-', '_', '/', '+', '.', '|', '*'.",
              call. = FALSE)
    ##
    if (!compmatch(labels, ":"))
      sep.trts <- ":"
    else if (!compmatch(labels, "-"))
      sep.trts <- "-"
    else if (!compmatch(labels, "_"))
      sep.trts <- "_"
    else if (!compmatch(labels, "/"))
      sep.trts <- "/"
    else if (!compmatch(labels, "+"))
      sep.trts <- "+"
    else if (!compmatch(labels, "."))
      sep.trts <- "-"
    else if (!compmatch(labels, "|"))
      sep.trts <- "|"
    else if (!compmatch(labels, "*"))
      sep.trts <- "*"
    else
      stop("All predefined separators (':', '-', '_', '/', '+', '.', '|', '*') ",
           "are used in at least one treatment label.\n   ",
           "Please specify a different character that should be used ",
           " as separator (argument 'sep.trts').",
           call. = FALSE)
  }
  ##
  if (reference.group != "")
    reference.group <- setref(reference.group, labels)
  ##
  rm(labels)


  ##
  ##
  ## (4) Additional checks
  ##
  ##
  if (any(treat1 == treat2))
    stop("Treatments must be different (arguments 'treat1' and 'treat2').",
         call. = FALSE)
  ##
  if (length(studlab) != 0)
    studlab <- as.character(studlab)
  else {
    if (warn)
      warning("No information given for argument 'studlab'. ",
              "Assuming that comparisons are from independent studies.",
              call. = FALSE)
    studlab <- seq(along = event1)
  }
  ##
  ## Check for correct number of comparisons
  ##
  tabnarms <- table(studlab)
  sel.narms <- !is.wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  ##
  if (sum(sel.narms) == 1)
    stop(paste0("Study '", names(tabnarms)[sel.narms],
                "' has a wrong number of comparisons.",
                "\n  Please provide data for all treatment comparisons",
                " (two-arm: 1; three-arm: 3; four-arm: 6, ...)."),
         call. = FALSE)
  if (sum(sel.narms) > 1)
    stop(paste0("The following studies have a wrong number of comparisons: ",
                paste(paste0("'", names(tabnarms)[sel.narms], "'"),
                      collapse = " - "),
                "\n  Please provide data for all treatment comparisons",
                " (two-arm: 1; three-arm: 3; four-arm: 6, ...)."),
         call. = FALSE)
  ##
  ## Check number of subgraphs
  ##
  n.subnets <- netconnection(treat1, treat2, studlab)$n.subnets
  ##
  if (n.subnets > 1)
    stop(paste0("Network consists of ", n.subnets, " separate sub-networks.\n",
                "  Use R function 'netconnection' to identify sub-networks."),
         call. = FALSE)
  #
  chknumeric(incr, min = 0)
  ##
  if (is.mh.nch & length(incr) > 1) {
    warning("Argument 'incr' must be a single value for ",
            "Mantel-Haenszel and common-effects non-central ",
            "hypergeometric method. Set to zero (default).",
            call. = FALSE)
    incr <- 0
  }
  ##
  if (all(incr == 0) & is.mh.nch & cc.pooled == TRUE)
    cc.pooled <- FALSE
  ##
  lengthunique <- function(x) length(unique(x))


  ##
  ##
  ## (5) Create analysis dataset
  ##
  ##
  dat.long <- rbind(data.frame(studlab, treat = treat1, event = event1, n = n1,
                               .order,
                               stringsAsFactors = FALSE),
                    data.frame(studlab, treat = treat2, event = event2, n = n2,
                               .order,
                               stringsAsFactors = FALSE))
  ##
  dat.wide <- data.frame(studlab = studlab, treat1 = treat1, treat2 = treat2,
                         event1 = event1, n1 = n1, event2 = event2, n2 = n2,
                         .order,
                         stringsAsFactors = FALSE)
  dat.wide <- dat.wide[order(dat.wide$studlab,
                             dat.wide$treat1, dat.wide$treat2), ]
  ##
  ## Check for correct number of treatments
  ##
  d.long <- unique(dat.long[, c("studlab", "treat", "event", "n")])
  d.long <- d.long[order(d.long$studlab, d.long$treat), , drop = FALSE]
  ##
  tabnarms <- table(d.long$studlab, d.long$treat)
  larger.one <- function(x) any(x > 1)
  sel.study <- apply(tabnarms, 1, larger.one)
  ##
  if (sum(sel.study) == 1) {
    if (details.chkdata) {
      cat("Inconsistent number of events and participants:\n")
      d.l <- d.long[d.long$studlab == rownames(tabnarms)[sel.study], ]
      prmatrix(d.l,
               quote = FALSE, right = TRUE,
               rowlab = rep("", nrow(d.l)))
    }
    stop(paste0("Study '", rownames(tabnarms)[sel.study],
                "' has inconsistent data. Please check the original data."),
         call. = FALSE)
  }
  if (sum(sel.study) > 1) {
    if (details.chkdata) {
      cat("Inconsistent number of events and participants:\n")
      d.l <- d.long[d.long$studlab %in% rownames(tabnarms)[sel.study], ]
      prmatrix(d.l,
               quote = FALSE, right = TRUE,
               rowlab = rep("", nrow(d.l)))
    }
    ##
    stop(paste0("The following studies have inconsistent data: ",
                paste(paste0("'", rownames(tabnarms)[sel.study], "'"),
                      collapse = " - "),
                "\n  Please check the original data."),
         call. = FALSE)
  }
  ##
  get.designs <- function(x) {
    net1 <-
      netmeta(pairwise(studlab = x$studlab, treat = x$treat,
                       event = x$event, n = x$n,
                       sm = "RD"),
              warn = FALSE)
    ##
    if (net1$n == 2)
      res <- data.frame(studlab = net1$studlab, design = net1$designs)
    else
      res <- nma_krahn(net1)$studies[, c("studlab", "design")]
    ##
    res$design <- as.character(res$design)
    res <- unique(res)
    ##
    res
  }
  ##
  data$.drop <- rep(FALSE, nrow(data))
  ##
  ## Add variable 'non.event'
  ##
  dat.long$non.event <- dat.long$n - dat.long$event
  ##
  dat.wide$non.event1 <- dat.wide$n1 - dat.wide$event1
  dat.wide$non.event2 <- dat.wide$n2 - dat.wide$event2
  ##
  ## Remove duplicated rows from dataset (due to multi-arm studies)
  ##
  dupl <- duplicated(dat.long[, c("studlab", "treat", "event", "n")])
  ##
  if (any(dupl))
    dat.long <- dat.long[!dupl, ]
  ##
  rm(dupl)
  ##
  all.designs <- get.designs(dat.long)
  ##
  ## Add variable 'design' with study design
  ##
  tdat.design <- get.designs(dat.long)
  ##
  dat.long <- merge(dat.long, tdat.design, by = "studlab",
                    all.x = TRUE)
  dat.wide <- merge(dat.wide, tdat.design, by = "studlab",
                    all.x = TRUE)
  ##
  dat.long <- dat.long[order(dat.long$.order), ]
  dat.wide <- dat.wide[order(dat.wide$.order), ]
  ##
  names(tdat.design) <- c("studlab", ".design")
  ##
  data <- merge(data, tdat.design,
                by.x = ".studlab", by.y = "studlab",
                all.x = TRUE)
  data <- data[order(data$.order), ]
  ##
  rm(tdat.design)
  
  
  ##
  ##
  ## (6) Stage I: setting up the data (Efthimiou et al., 2019)
  ##
  ##
  ## Step i. Remove all-zero or all-event studies (only MH and NCH methods)
  ##
  if (is.mh.nch) {
    events.study <- with(dat.long, tapply(event, studlab, sum))
    nonevents.study <- with(dat.long, tapply(n - event, studlab, sum))
    ##
    if (any(events.study == 0) | any(nonevents.study == 0)) {
      zeroevents <- events.study == 0
      allevents <- nonevents.study == 0 & events.study != 0
      keep <- !(zeroevents | allevents)
      ##
      if (warn) {
        if (sum(zeroevents) == 1)
          warning("Study '", names(events.study)[zeroevents],
                  "' without any events excluded from network meta-analysis.",
                  call. = FALSE)
        else if (sum(zeroevents) > 1)
          warning("Studies without any events excluded ",
                  "from network meta-analysis: ",
                  paste(paste0("'", names(events.study)[zeroevents], "'"),
                        collapse = " - "),
                  call. = FALSE)
        ##
        if (sum(allevents) == 1)
          warning("Study '", names(nonevents.study)[allevents],
                  "' with all events excluded from network meta-analysis.",
                  call. = FALSE)
        else if (sum(allevents) > 1)
          warning("Studies with all events excluded ",
                  "from network meta-analysis: ",
                  paste(paste0("'", names(nonevents.study)[allevents], "'"),
                        collapse = " - "),
                  call. = FALSE)
      }
      ##
      dat.long <- dat.long[dat.long$studlab %in% names(events.study)[keep], ,
                           drop = FALSE]
      dat.wide <- dat.wide[dat.wide$studlab %in% names(events.study)[keep], ,
                           drop = FALSE]
      ##
      data$.drop <- data$.drop | data$studlab %in% names(events.study)[!keep]
      ##
      rm(zeroevents, allevents, keep)
    }
    ##
    rm(events.study)
  }
  ##
  ## Add variable 'study' with study numbers
  ##
  dat.long$study <- as.numeric(factor(dat.long$studlab,
                                      levels = unique(dat.long$studlab)))
  dat.wide$study <- as.numeric(factor(dat.wide$studlab,
                                      levels = unique(dat.wide$studlab)))
  ##
  ## Sort by study and treatment
  ##
  dat.long <- dat.long[order(dat.long$studlab, dat.long$treat), ]
  dat.wide <- dat.wide[order(dat.wide$studlab,
                              dat.wide$treat1, dat.wide$treat2), ]
  ##
  dat.long$incr <- 0
  dat.wide$incr <- 0
  data$.incr <- 0
  ##
  ## Step iii. Drop treatment arms without or all events from individual
  ##           designs (argument 'cc.pooled' is FALSE) or add
  ##           increment if argument 'cc.pooled' is TRUE and argument
  ##           'incr' is larger than zero
  ##
  if (is.mh.nch) {
    ##
    events.arm <- with(dat.long, tapply(event, list(design, treat), sum))
    nonevents.arm <- with(dat.long, tapply(n - event, list(design, treat), sum))
    ##
    if (any(events.arm == 0, na.rm = TRUE) |
        any(nonevents.arm == 0, na.rm = TRUE)) {
      ##
      ## Identify entries without events
      ##
      zeroevents <- events.arm == 0
      zerocells <- as.data.frame(which(zeroevents, arr.ind = TRUE),
                                 stringsAsFactors = FALSE)
      ##
      zerocells$design <- rownames(zeroevents)[zerocells$row]
      zerocells$treat <- colnames(zeroevents)[zerocells$col]
      ##
      zero.long <- rep(0, nrow(dat.long))
      zero.wide <- rep(0, nrow(dat.wide))
      zero.data <- rep(0, nrow(data))
      ##
      for (i in seq_along(zerocells$design)) {
        if (!cc.pooled)
          zero.long <- zero.long + (dat.long$design == zerocells$design[i] &
                                    dat.long$treat == zerocells$treat[i])
        else
          zero.long <- zero.long + (dat.long$design == zerocells$design[i])
        ##
        zero.wide <- zero.wide + (dat.wide$design == zerocells$design[i] &
                                  (dat.wide$treat1 == zerocells$treat[i] |
                                   dat.wide$treat2 == zerocells$treat[i]))
        ##
        zero.data <- zero.data + (data$.design == zerocells$design[i] &
                                  (data$.treat1 == zerocells$treat[i] |
                                   data$.treat2 == zerocells$treat[i]))
      }
      ##
      zero.long <- zero.long > 0
      zero.wide <- zero.wide > 0
      zero.data <- zero.data > 0
      ##
      ## Identify entries with all events
      ##
      allevents <- nonevents.arm == 0
      allcells <- as.data.frame(which(allevents, arr.ind = TRUE),
                                 stringsAsFactors = FALSE)
      ##
      allcells$design <- rownames(allevents)[allcells$row]
      allcells$treat <- colnames(allevents)[allcells$col]
      ##
      all.long <- rep(0, nrow(dat.long))
      all.wide <- rep(0, nrow(dat.wide))
      all.data <- rep(0, nrow(data))
      ##
      for (i in seq_along(allcells$design)) {
        if (!cc.pooled)
          all.long <- all.long + (dat.long$design == allcells$design[i] &
                                  dat.long$treat == allcells$treat[i])
        else
          all.long <- all.long + (dat.long$design == allcells$design[i])
        ##
        all.wide <- all.wide + (dat.wide$design == allcells$design[i] &
                                (dat.wide$treat1 == allcells$treat[i] |
                                 dat.wide$treat2 == allcells$treat[i]))
        ##
        all.data <- all.data + (data$.design == allcells$design[i] &
                                (data$.treat1 == allcells$treat[i] |
                                 data$.treat2 == allcells$treat[i]))
      }
      ##
      all.long <- all.long > 0
      all.wide <- all.wide > 0
      all.data <- all.data > 0
      ##
      ## Identify entries with sparse binary data (i.e., with zero or
      ## all events)
      ##
      sparse.long <- zero.long | all.long
      sparse.wide <- zero.wide | all.wide
      sparse.data <- zero.data | all.data
      ##
      if (!cc.pooled) {
        if (warn) {
          if (sum(zeroevents, na.rm = TRUE) == 1)
            warning("Treatment arm '", zerocells$treat,
                    "' without events in design '",
                    zerocells$design, "' excluded from network meta-analysis.",
                    call. = FALSE)
          else if (sum(zeroevents, na.rm = TRUE) > 1)
            warning("Treatment arms without events in a design excluded ",
                    "from network meta-analysis:\n    ",
                    paste0("'",
                           paste0(paste0(zerocells$treat, " in "),
                                  zerocells$design),
                           "'",
                           collapse = " - "),
                    call. = FALSE)
          ##
          if (sum(allevents, na.rm = TRUE) == 1)
            warning("Treatment arm '", allcells$treat,
                    "' with all events in design '",
                    allcells$design, "' excluded from network meta-analysis.",
                    call. = FALSE)
          else if (sum(allevents, na.rm = TRUE) > 1)
            warning("Treatment arms with all events in a design excluded ",
                    "from network meta-analysis:\n    ",
                    paste0("'",
                           paste0(paste0(allcells$treat, " in "),
                                  allcells$design),
                           "'",
                           collapse = " - "),
                    call. = FALSE)
        }
        ##
        dat.long <- dat.long[!sparse.long, , drop = FALSE]
        dat.wide <- dat.wide[!sparse.wide, , drop = FALSE]
        data$.drop <- data$.drop | sparse.data
      }
      else {
        dat.long$incr[sparse.long] <- incr
        ##
        dat.long$event[sparse.long] <- dat.long$event[sparse.long] + incr
        dat.long$non.event[sparse.long] <- dat.long$non.event[sparse.long] + incr
        dat.long$n[sparse.long] <- dat.long$n[sparse.long] + 2 * incr
        ##
        dat.wide$incr[sparse.wide] <- incr
        ##
        data$.incr[sparse.data] <- incr
      }
      ##
      rm(zeroevents, zerocells, zero.long, zero.wide, zero.data,
         allevents, allcells, all.long, all.wide, all.data,
         sparse.long, sparse.wide, sparse.data)
    }
    ##
    rm(events.arm, nonevents.arm)
  }
  ##
  ## Step iv. Remove designs with single treatment arm from dataset
  ##
  if (is.mh.nch) {
    d.single <- with(dat.long, tapply(treat, design, lengthunique))
    ##
    if (any(d.single == 1, na.rm = TRUE)) {
      single <- !is.na(d.single) & d.single == 1
      design.single <- names(d.single)[single]
      ##
      if (warn)
        if (sum(single) == 1)
          warning("Design '", design.single,
                  "' with single treatment arm excluded ",
                  "from network meta-analysis.",
                  call. = FALSE)
        else
          warning("Designs with single treatment arm excluded ",
                  "from network meta-analysis: ",
                  paste(paste0("'", design.single, "'"),
                        collapse = " - "),
                  call. = FALSE)
      ##
      dat.long <- dat.long[!(dat.long$design %in% design.single), , drop = FALSE]
      dat.wide <- dat.wide[!(dat.wide$design %in% design.single), , drop = FALSE]
      ##
      data$.drop <- data$.drop | data$.design %in% design.single
      ##
      rm(single, design.single)
    }
    ##
    rm(d.single)
  }
  ##
  dat.long <- dat.long[, c("studlab", "treat",
                           "event", "non.event", "n", "incr",
                           "design", "study", ".order")]
  dat.wide <- dat.wide[, c("studlab", "treat1", "treat2",
                           "event1", "event2", "non.event1", "non.event2",
                           "n1", "n2", "incr",
                           "design", "study", ".order")]
  ##
  dat.long$studlab <- as.character(dat.long$studlab)
  dat.long$design <- as.character(dat.long$design)
  dat.long$treat <- as.character(dat.long$treat)
  ##
  dat.wide$studlab <- as.character(dat.wide$studlab)
  dat.wide$design <- as.character(dat.wide$design)
  dat.wide$treat1 <- as.character(dat.wide$treat1)
  dat.wide$treat2 <- as.character(dat.wide$treat2)
  ##
  o <- order(dat.wide$.order)
  studlab <- dat.wide$studlab[o]
  treat1 <- dat.wide$treat1[o]
  treat2 <- dat.wide$treat2[o]
  ##
  event1 <- dat.wide$event1[o]
  event2 <- dat.wide$event2[o]
  n1 <- dat.wide$n1[o]
  n2 <- dat.wide$n2[o]
  ##
  trts <- sort(unique(c(treat1, treat2)))
  n.treat <- length(trts)
  n.studies <- length(unique(studlab))
  m <- length(studlab)
  ##
  treat1.pos <- match(treat1, trts)
  treat2.pos <- match(treat2, trts)
  ##
  ## Calculate adjacency matrix
  ##
  B <- createB(treat1.pos, treat2.pos, ncol = n.treat)
  ##
  M <- t(B) %*% B    # unweighted Laplacian matrix
  D <- diag(diag(M)) # diagonal matrix
  A <- D - M         # adjacency matrix (n x n)
  ##
  rownames(A) <- colnames(A) <- trts
  ##
  rm(B, M, D)
  ##
  ## Empty matrices for results
  ##
  NAmatrix <- matrix(NA, nrow = n.treat, ncol = n.treat)
  rownames(NAmatrix) <- colnames(NAmatrix) <- trts
  ##
  TE.common <- seTE.common <- NAmatrix
  TE.direct.common <- seTE.direct.common <- NAmatrix
  Q.direct <- tau2.direct <- I2.direct <- NAmatrix
  #
  if (is.lrp) {
    TE.direct.random <- seTE.direct.random <- phi.direct <- NAmatrix
  }
  
  
  ##
  ##
  ## (7) Conduct classic network meta-analysis using inverse variance
  ##     method
  ##
  ##
  dat.iv <- dat.long[order(dat.long$.order), ]
  ##
  if (method == "Inverse") {
    if (missing(warn))
      warn.iv <- gs("warn")
    else
      warn.iv <- warn
    incr.iv <- incr
    allstudies.iv <- allstudies
  }
  else {
    warn.iv <- FALSE
    incr.iv <- incr * cc.pooled
    allstudies.iv <- cc.pooled
  }
  ##
  p.iv <- pairwise(studlab = dat.iv$studlab,
                   treat = dat.iv$treat,
                   event = dat.iv$event,
                   n = dat.iv$n,
                   data = dat.iv,
                   sm = sm,
                   incr = incr.iv,
                   method.incr = method.incr,
                   allstudies = allstudies.iv)
  ##
  net.iv <- netmeta(p.iv,
                    level = level, level.ma = level.ma,
                    common = common, random = random,
                    prediction = prediction, level.predict = level.predict,
                    reference.group = reference.group,
                    baseline.reference = baseline.reference,
                    all.treatments = all.treatments,
                    seq = seq,
                    tau.preset = tau.preset,
                    tol.multiarm = tol.multiarm,
                    tol.multiarm.se = tol.multiarm.se,
                    details.chkmultiarm = details.chkmultiarm,
                    sep.trts = sep.trts,
                    nchar.trts = nchar.trts,
                    overall.hetstat = overall.hetstat,
                    backtransf = backtransf,
                    title = title,
                    keepdata = keepdata,
                    warn = warn.iv)
  ##
  if (method == "Inverse")
    return(net.iv)
  else {
    net.iv$Cov.common[!is.na(net.iv$Cov.common)] <- NA
    net.iv$Cov.random[!is.na(net.iv$Cov.random)] <- NA
  }
  
  
  ##
  ##
  ## (8) Stage 2: Direct meta-analyses per design (MH and NCH methods)
  ##
  ##
  n.d <- tapply(dat.long$treat,   dat.long$design, lengthunique)
  k.d <- tapply(dat.long$studlab, dat.long$design, lengthunique)
  ##
  designs <- names(k.d)
  d <- length(designs)
  seq.d <- seq_len(d)
  ##
  dat.design <- vector("list", d)
  names(dat.design) <- designs
  ##
  for (i in seq.d)
    dat.design[[i]] <- dat.long[dat.long$design == designs[i], , drop = FALSE]
  ##
  if (method == "MH") {
    ##
    ## MH method
    ##
    treat.per.design <- vector("list", d)
    studies.in.design <- vector("list", d)
    ##
    c.xy <- vector("list", d)
    d.xy <- vector("list", d)
    C.xy <- vector("list", d)
    L.xy <- vector("list", d)
    ##
    U.xyy <- vector("list", d)
    U.xyz <- vector("list", d)
    ##
    t.pl <- vector("list", d)
    ##
    ps1 <- ps2 <- ps3 <- ps4 <- vector("list", d)
    ##
    L.bar.xy <- vector("list", d)
    ##
    U.plus.xx <- vector("list", d)
    ##
    U.new <- vector("list", d)
    ##
    U.plus.xy <- vector("list", d)
    ##
    Var.Lbar <- vector("list", d)
    ##
    CoVar.Lbar <- vector("list", d)
    ##
    length.y <- sum(n.d - 1)
    ##
    y <- c(rep(0, length.y))
    ##
    V1 <- matrix(rep(0, length.y * length.y), nrow = length.y)
    V2 <- matrix(rep(0, length.y * length.y), nrow = length.y)
    ##
    X <- matrix(0, nrow = length.y, ncol = n.treat - 1)
    ##
    counter1 <- counter2 <- counter3 <- 0
    ##
    list1 <- matrix(0, nrow = length.y, ncol = 2)
    ##
    N.j <- rep(0, d)
    ##
    if (d > 1)
      for (i in 2:d)
        N.j[i] <- N.j[i - 1] + n.d[i - 1] - 1
    ##
    ## Loop for designs
    ##
    for (i in seq.d) {
      treat.per.design[[i]] <- unique(dat.design[[i]]$treat)
      studies.in.design[[i]] <- unique(dat.design[[i]]$study)
      ##
      ## List of studies in each design
      ##
      ## Recode the studies in each design
      ##
      for (j in seq_along(dat.design[[i]]$study))
        for (k in seq_len(k.d[i]))
          if (dat.design[[i]]$study[j] == studies.in.design[[i]][k])
            dat.design[[i]]$id.s[j] <- k
      ##
      ## Recode the treatments in each design
      ##
      for (j in seq_along(dat.design[[i]]$study))
        for (k in seq_len(n.d[i]))
          if (dat.design[[i]]$treat[j] == treat.per.design[[i]][k])
            dat.design[[i]]$id.t[j] <- k
      ##
      ## Calculate total number of patients per studies
      ##
      for (j in seq_along(dat.design[[i]]$studlab))
        dat.design[[i]]$n.study[j] <-
          with(dat.design[[i]], sum(n[which(study == study[j])]))
      ##
      ## Define c.xy and d.xy
      ##
      c.xy[[i]] <- array(rep(0, k.d[i] * n.d[i] * n.d[i]),
                         dim = c(n.d[i], n.d[i], k.d[i]))
      ##
      d.xy[[i]] <- array(rep(0, k.d[i] * n.d[i] * n.d[i]),
                         dim = c(n.d[i], n.d[i], k.d[i]))
      ##
      for (st in seq_len(k.d[i]))
        for (t1 in seq_len(n.d[i]))
          for (t2 in seq_len(n.d[i])) {
            c.xy[[i]][t1, t2, st] <-
              with(dat.design[[i]],
                   sum(    event[which(id.s == st & id.t == t1)]) *
                   sum(non.event[which(id.s == st & id.t == t2)]) /
                   n.study[id.s == st][1])
            ##
            d.xy[[i]][t1, t2, st] <-
              with(dat.design[[i]],
              (sum(    event[which(id.s == st & id.t == t1)]) +
               sum(non.event[which(id.s == st & id.t == t2)])) /
              n.study[id.s == st][1])
          }
      ##
      ## Define C.xy
      ##
      C.xy[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      ##
      for (j in seq_len(n.d[i]))
        for (k in seq_len(n.d[i]))
          C.xy[[i]][j, k] <- sum(c.xy[[i]][j, k, ])
      ##
      ## Define L.xy
      ##
      L.xy[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      for (j in seq_len(n.d[i]))
        for (k in seq_len(n.d[i]))
          L.xy[[i]][j, k] <- log(C.xy[[i]][j, k] / C.xy[[i]][k, j])
      ##
      ## Calculate the variance of L.xy (dimension n.d[i] x n.d[i])
      ##
      U.xyy[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      for (j in seq_len(n.d[i]))
        for (k in seq_len(n.d[i]))
          U.xyy[[i]][j, k] <-
            sum(c.xy[[i]][j, k, ] * d.xy[[i]][j, k, ]) /
            (2 * C.xy[[i]][j, k]^2) +
            sum(c.xy[[i]][j, k, ] * d.xy[[i]][k, j, ] + c.xy[[i]][k, j, ] *
                d.xy[[i]][j, k, ]) / (2 * C.xy[[i]][j, k] * C.xy[[i]][k, j]) +
            sum(c.xy[[i]][k, j, ] * d.xy[[i]][k, j, ]) / (2 * C.xy[[i]][k, j]^2)
      ##
      ## Calculate the covariance matrix U.xyz
      ##
      t.pl[[i]] = array(rep(0, k.d[i] * n.d[i] * n.d[i]),
                        dim = c(n.d[i], n.d[i], k.d[i]))
      ##
      for (j in seq_len(k.d[i]))
        t.pl[[i]][ , , j] <- with(dat.design[[i]], n.study[id.s == j][1])
      ##
      U.xyz[[i]] <- array(rep(0, n.d[i] * n.d[i] * n.d[i]),
                          dim = c(n.d[i], n.d[i], n.d[i]))
      ##
      ## Per study ...
      ##
      for (t1 in seq_len(n.d[i])) {
        for (t2 in seq_len(n.d[i])) {
          for (t3 in seq_len(n.d[i])) {
            for (st in seq_len(k.d[i])) {
              sel1 <- with(dat.design[[i]], which(id.s == st & id.t == t1))
              sel2 <- with(dat.design[[i]], which(id.s == st & id.t == t2))
              sel3 <- with(dat.design[[i]], which(id.s == st & id.t == t3))
              ##
              ps1[[i]][st] <-
                with(dat.design[[i]],
                     event[sel1] /
                     (t.pl[[i]][1, 1, st])^2 *
                     non.event[sel2] * non.event[sel3] *
                     (t1 != t2) * (t1 != t3) * (t2 != t3))
              ##
              ps2[[i]][st] <-
                with(dat.design[[i]],
                     n[sel1] / (t.pl[[i]][1, 1, st])^2 *
                     non.event[sel2] * event[sel3] *
                     (t1 != t2) * (t1 != t3) * (t2 != t3))
              ##
              ps3[[i]][st] <-
                with(dat.design[[i]],
                     n[sel1] / (t.pl[[i]][1, 1, st])^2 *
                     event[sel2] * non.event[sel3] *
                     (t1 != t2) * (t1 != t3) * (t2 != t3))
              ##
              ps4[[i]][st] <-
                with(dat.design[[i]],
                     non.event[sel1] / (t.pl[[i]][1, 1, st])^2 *
                     event[sel2] * event[sel3] *
                     (t1 != t2) * (t1 != t3) * (t2 != t3))
            }
            ##
            U.xyz[[i]][t1, t2, t3] <-
              sum(ps1[[i]][]) / (3 * C.xy[[i]][t1, t2] * C.xy[[i]][t1, t3]) +
              sum(ps2[[i]][]) / (3 * C.xy[[i]][t1, t2] * C.xy[[i]][t3, t1]) +
              sum(ps3[[i]][]) / (3 * C.xy[[i]][t2, t1] * C.xy[[i]][t1, t3]) +
              sum(ps4[[i]][]) / (3 * C.xy[[i]][t2, t1] * C.xy[[i]][t3, t1])
          }
        }
      }
      ##
      ## Calculate L.bar.xy
      ##
      L.bar.xy[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          L.bar.xy[[i]][t1, t2] <-
            (sum(L.xy[[i]][t1, ]) - sum(L.xy[[i]][t2, ])) / n.d[i]
      ##
      ## Calculate U.plus.xx
      ##
      for (t1 in seq_len(n.d[i]))
        U.xyy[[i]][t1, t1] <- 0
      ##
      U.plus.xx[[i]] <- rep(0, n.d[i])
      for (t1 in seq_len(n.d[i]))
        U.plus.xx[[i]][t1] <-
          sum(U.xyy[[i]][t1, seq_len(n.d[i])]) + sum(U.xyz[[i]][t1, , ])
      ##
      ## Calculate U.plus.xy
      ##
      U.new[[i]] <- array(rep(0, n.d[i] * n.d[i] * n.d[i]),
                          dim = c(n.d[i], n.d[i], n.d[i]))
      ##
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          for (t3 in seq_len(n.d[i]))
            U.new[[i]][t1, t2, t3] <-
              (t1 != t2) * (t1 != t3) * (t2 != t3) *
              U.xyz[[i]][t1, t2, t3] +
              (t1 != t2) * (t2 == t3) * U.xyy[[i]][t1, t2]
      ##
      U.plus.xy[[i]] = matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      ##
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          U.plus.xy[[i]][t1, t2] <-
            sum(U.new[[i]][seq_len(n.d[i]), t1, t2]) -
            sum(U.new[[i]][t1, t2, seq_len(n.d[i])]) -
            sum(U.new[[i]][t2, t1, seq_len(n.d[i])]) +
            U.new[[i]][t1, t2, t2]
      ##
      ## Variance of L.bar.xy
      ##
      Var.Lbar[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          Var.Lbar[[i]][t1, t2] <-
            (U.plus.xx[[i]][t1] - 2 * U.plus.xy[[i]][t1, t2] +
             U.plus.xx[[i]][t2]) / n.d[i]^2
      ##
      ## Covariance of L.bar.xy. Only a subset of covariances are
      ## calculate here, i.e. cov(L.bar_(1, t1), L.bar_(1, t2))
      ##
      CoVar.Lbar[[i]] <- matrix(rep(0, n.d[i]^2), nrow = n.d[i])
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          CoVar.Lbar[[i]][t1, t2] <-
            (U.plus.xx[[i]][1] - U.plus.xy[[i]][1, t2] -
             U.plus.xy[[i]][t1, 1] +
             U.plus.xy[[i]][t1, t2]) / n.d[i]^2
      ##
      ## Define matrix X
      ##
      for (j in seq_along(dat.design[[i]]$studlab))
        for (k in seq_along(trts))
          if (dat.design[[i]]$treat[j] == trts[k])
            dat.design[[i]]$id.t2[j] <- k
      ##
      ##
      ## Create y, the vector of treatment effects from each study. Only
      ## effects vs the first treatment are needed. y is coded as OR
      ## 1vsX
      ##
      for (j in seq_len(n.d[i] - 1))
        y[j + counter1] <- L.bar.xy[[i]][1, j + 1]
      ##
      counter1 <- counter1 + n.d[i] - 1
      ##
      ## Create V, the matrix of covariances of treatment effects from
      ## each study. Only covariances of effects vs the first treatment
      ## are needed.
      ##
      for (j in seq_len(n.d[i] - 1))
        for (k in seq_len(n.d[i] - 1))
          V1[j + counter2, k + counter2] <- (k == j) * Var.Lbar[[i]][1, j + 1]
      ##
      counter2 <- counter2 + n.d[i] - 1
      ##
      for (j in 2:(n.d[i]))
        for (k in 2:(n.d[i]))
          V2[j + counter3 - 1, k + counter3 - 1] <-
            (k != j) * CoVar.Lbar[[i]][j, k]
      ##
      counter3 <- counter3 + n.d[i] - 1
      ##
      for (j in seq_len(n.d[i] - 1)) {
        list1[N.j[i] + j, 1] <- dat.design[[i]]$id.t2[[1]]
        list1[N.j[i] + j, 2] <- dat.design[[i]]$id.t2[[j + 1]]
      }
      ##
      ## Drop unnecessary variables
      ##
      dat.design[[i]]$id.s <- NULL
      dat.design[[i]]$id.t <- NULL
      dat.design[[i]]$id.t2 <- NULL
      dat.design[[i]]$n.study <- NULL
    }
    ##
    V <- V1 + V2
    ##
    basic.contrasts <- 2:n.treat
    ##
    for (i in 1:length.y)
      for (k in 1:(n.treat - 1)) {
        if (list1[i, 1] == basic.contrasts[k])
          X[i, k] = -1
        if (list1[i, 2] == basic.contrasts[k])
          X[i, k] = 1
      }
    ##
    W <- solve(V)
    ##
    TE.basic <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
  }
  else if (method == "NCH") {
    ##
    ## NCH method
    ##
    dat.long$ttt <- 0
    ##
    for (k in 1:n.treat)
      for (i in seq_along(dat.long$studlab))
        if (dat.long$treat[i] == trts[k])
          dat.long$ttt[i] <- k
    ##
    dat.long$treat <- as.numeric(dat.long$ttt)
    dat.long$ttt <- NULL
    ##
    dat.long$treat <- as.numeric(dat.long$treat)
    ##
    dat.long <- dat.long[order(dat.long$studlab, dat.long$treat), ] # necessary ???
    ##
    for (i in unique(dat.long$studlab)) {
      counter <- 1
      for (j in seq_along(dat.long$studlab))
        if (dat.long$studlab[j] == i) {
          dat.long$count[j] <- counter
          counter <- counter + 1
        }
    }
    ##
    max.arms <- max(dat.long$count)
    ##
    d1 <- reshape(dat.long, idvar = "studlab",
                  timevar = "count", direction = "wide")
    ##
    d1 <- d1[order(d1$studlab), ]
    ##
    d1$N1 <- 0 # d1$n.all = 0
    d1$narms <- 0
    ##
    for (i in seq_len(max.arms)) {
      ev <- paste0("event.", i)
      tot <- paste0("n.", i)
      d1$N1 <- rowSums(cbind(d1$N1, d1[, colnames(d1) == ev]),
                       na.rm = TRUE)
      ## d1$n.all <- rowSums(cbind(d1$n.all, d1[, colnames(d1) == tot]), na.rm = TRUE)
    }
    ##
    for (i in seq_along(d1$studlab))
      for (k in 1:max.arms) {
        ev <- paste0("event.", k)
        if (!is.na(d1[, colnames(d1) == ev][i]))
          d1$narms[i] <- k
      }
    ##
    for (i in seq_along(colnames(d1)))
      for (k in seq_along(colnames(d1))) {
        if (colnames(d1)[k] == paste0("treat.", i))
          colnames(d1)[k] = paste0("t", i)
        ##
        if (colnames(d1)[k] == paste0("event.", i))
          colnames(d1)[k] = paste0("r", i)
        ##
        if (colnames(d1)[k] == paste0("n.", i))
          colnames(d1)[k] = paste0("n", i)
      }
    ##
    ## Likelihood function
    ##
    myLik1 <- function(mypar) {
      x <- mypar
      ##
      myLogLik1 <- myLogLik2 <- myLogLik3 <- rep(0, length(d1$studlab))
      ##
      for (i in seq_along(d1$studlab)) {
        for (k in 2:d1$narms[i]) {
          myLogLik1[i] <-
            myLogLik1[i] +
            d1[, colnames(d1) == paste0("r", k)][i] *
            (x[d1[, colnames(d1) == paste0("t", k)][i]] -
             x[d1$t1[i]] * (d1$t1[i] != 1))
          ##
          myLogLik2[i] <-
            myLogLik2[i] +
            d1[, colnames(d1) == paste0("n", k)][i] *
            exp(x[d1[, colnames(d1) == paste0("t", k)][i]] -
                x[d1$t1[i]] * (d1$t1[i] != 1))
          ##
          myLogLik3[i] <- -d1$N1[i] * log(d1$n1[i] + myLogLik2[i])
        }
      }
      ##
      myLogLik <- sum(myLogLik1 + myLogLik3)
      ##
      myLogLik
    }
    ##
    opt <- optim(rep(0, n.treat), myLik1, method = "L-BFGS-B",
                 lower = -Inf, upper = Inf,
                 control = list(fnscale = -1, maxit = 10000),
                 hessian = TRUE)
    ##
    W <- solve(-opt$hessian[2:n.treat, 2:n.treat])
    ##
    TE.basic <- -(opt$par)[2:n.treat]
  }
  else if (is.lrp) {
    #
    # Penalized logistic regression
    #
    is_installed_package("brglm2", "netmetabin", "method", " = \"LRP\"")
    #
    phi <- 1
    text.formula <- "cbind(event, non.event) ~ as.factor(treat)"
    #
    if (length(unique(dat.long$studlab)) > 1)
      text.formula <- paste(text.formula, "+ as.factor(studlab)")
    #
    res.lrp <- glm(as.formula(text.formula), data = dat.long,
                   family = binomial(link = "logit"),
                   method = brglm2::brglmFit, type = "MPL_Jeffreys",
                   ...)
    #
    if (length(unique(dat.long$studlab)) > 1)
      phi <- phi(res.lrp)
    #
    # Basic estimates (all treats vs reference)
    #
    ests.lrp <- summary(res.lrp)$coefficients
    ests.lrp <-
      ests.lrp[grepl("as.factor(treat)", rownames(ests.lrp), fixed = TRUE), ]
    rownames(ests.lrp) <- 
      gsub("as.factor(treat)", "", rownames(ests.lrp), fixed = TRUE)
    #
    trts.basic <- rownames(ests.lrp)
    #
    # Get the variance-covariance matrix for the basic estimates
    #
    var_cov_common <- vcov(res.lrp)
    var_cov_common <-
      var_cov_common[
        grepl("as.factor(treat)", rownames(var_cov_common), fixed = TRUE),
        grepl("as.factor(treat)", colnames(var_cov_common), fixed = TRUE)]
    #
    rownames(var_cov_common) <-
      gsub("as.factor(treat)", "", rownames(var_cov_common), fixed = TRUE)
    colnames(var_cov_common) <-
      gsub("as.factor(treat)", "", colnames(var_cov_common), fixed = TRUE)
    #
    # Inflate the elements of the whole variance covariance matrix from the
    # common effect model to get the respective matrix for the random effects
    # model
    #
    var_cov_random <- var_cov_common * phi
    #
    TE.basic <- ests.lrp[, 1]
    #
    # Use the basic2all() function to calculate the functional parameters using
    # the basic parameters
    #
    all_ests_common <- basic2all(trts, TE.basic, var_cov_common)
    #
    all_ests_random <- basic2all(trts, TE.basic, phi * var_cov_common)
    #
    TE.common <- all_ests_common$TEs
    seTE.common <- all_ests_common$seTEs
    #
    TE.random <- all_ests_random$TEs
    seTE.random <- all_ests_random$seTE
  }
  #
  if (!is.lrp) {
    #
    # Define H matrix
    #
    H <- matrix(0,
                nrow = n.treat * ((n.treat - 1) / 2),
                ncol = n.treat - 1)
    #
    diag(H) <- 1
    #
    if (n.treat > 2) {
      t1 <- c()
      t2 <- c()
      for (i in 2:(n.treat - 1))
        for (j in (i + 1):(n.treat)) {
          t1 <- rbind(t1, i)
          t2 <- rbind(t2, j)
        }
      #
      h1 <- matrix(c(t1, t2), nrow = nrow(t1))
      #
      for (i in 1:((n.treat - 1) * (n.treat - 2) / 2))
        for (j in 1:(n.treat - 1))
          H[i + n.treat - 1, j] <- -(h1[i, 1] == j + 1) + (h1[i, 2] == j + 1)
    }
    #
    colnames(H) <- trts[-1]
    #
    # Common effects matrix
    #
    d.hat <- H %*% TE.basic
    #
    TE.common[lower.tri(TE.common, diag = FALSE)] <-  d.hat
    TE.common <- t(TE.common)
    TE.common[lower.tri(TE.common, diag = FALSE)] <- -d.hat
    diag(TE.common) <- 0
    #
    # Matrix with standard errors
    #
    if (method == "MH")
      cov.d.hat <- H %*% solve(t(X) %*% W %*% X) %*% t(H)
    else
      cov.d.hat <- H %*% W %*% t(H)
    ##
    seTE.common[lower.tri(seTE.common, diag = FALSE)] <- sqrt(diag(cov.d.hat))
    seTE.common <- t(seTE.common)
    seTE.common[lower.tri(seTE.common, diag = FALSE)] <- sqrt(diag(cov.d.hat))
    diag(seTE.common) <- 0
  }
  else
    H <- matrix(NA,
                nrow = n.treat * ((n.treat - 1) / 2),
                ncol = n.treat - 1)
  #
  # Confidence intervals
  #
  ci.f <- ci(TE.common, seTE.common, level = level.ma)
  #
  if (is.lrp)
    ci.r <- ci(TE.random, seTE.random, level = level.ma)
  ##
  ## Inconsistency global
  ##
  if (method == "MH") {
    Q <- as.vector(t(y - X %*% TE.basic) %*% solve(V) %*%
                    (y - X %*% TE.basic))
    df.Q <- sum(n.d - 1) - (n.treat - 1)
    pval.Q <- pvalQ(Q, df.Q)
  }
  else {
    Q <- NA
    df.Q <- NA
    pval.Q <- NA
  }
  
  
  ##
  ##
  ## (9) Inconsistency evaluation: direct MH estimates
  ##
  ##
  B.matrix <- createB(treat1.pos, treat2.pos)
  B.matrix <- unique(B.matrix)
  colnames(B.matrix) <- trts
  ##
  for (i in seq_len(nrow(B.matrix))) {
    sel.treat1 <- colnames(B.matrix)[B.matrix[i, ] ==  1]
    sel.treat2 <- colnames(B.matrix)[B.matrix[i, ] == -1]
    selstud <- treat1 == sel.treat1 & treat2 == sel.treat2
    ##
    if (!is.lrp)
      m.i <- metabin(event1, n1, event2, n2, subset = selstud,
                     method = "MH", sm = "OR",
                     incr = incr,
                     method.incr = method.incr,
                     allstudies = allstudies, MH.exact = !cc.pooled,
                     method.tau = "DL", method.tau.ci = "",
                     Q.Cochrane = FALSE,
                     warn.deprecated = FALSE)
    else
      m.i <- metabin(event1, n1, event2, n2, subset = selstud,
                     method = "LRP",
                     warn.deprecated = FALSE, ...)
    #
    TE.direct.common[sel.treat1, sel.treat2]   <- m.i$TE.common
    seTE.direct.common[sel.treat1, sel.treat2] <- m.i$seTE.common
    #
    TE.direct.common[sel.treat2, sel.treat1]   <- -m.i$TE.common
    seTE.direct.common[sel.treat2, sel.treat1] <- m.i$seTE.common
    #
    Q.direct[sel.treat1, sel.treat2] <- m.i$Q
    tau2.direct[sel.treat1, sel.treat2] <- m.i$tau2
    I2.direct[sel.treat1, sel.treat2] <- m.i$I2
    #
    Q.direct[sel.treat2, sel.treat1] <- m.i$Q
    tau2.direct[sel.treat2, sel.treat1] <- m.i$tau2
    I2.direct[sel.treat2, sel.treat1] <- m.i$I2
    #
    if (is.lrp) {
      TE.direct.random[sel.treat1, sel.treat2]   <- m.i$TE.random
      seTE.direct.random[sel.treat1, sel.treat2] <- m.i$seTE.random
      #
      TE.direct.random[sel.treat2, sel.treat1]   <- -m.i$TE.random
      seTE.direct.random[sel.treat2, sel.treat1] <- m.i$seTE.random
      #
      phi.direct[sel.treat1, sel.treat2] <- m.i$phi
      phi.direct[sel.treat2, sel.treat1] <- m.i$phi
    }
  }
  #
  rm(sel.treat1, sel.treat2, selstud, m.i)
  #
  ci.direct.c <- ci(TE.direct.common, seTE.direct.common, level = level.ma)
  #
  if (is.lrp)
    ci.direct.r <- ci(TE.direct.random, seTE.direct.random, level = level.ma)
  
  
  labels <- sort(unique(c(treat1, treat2)))
  ##
  if (!is.null(seq))
    seq <- setseq(seq, labels)
  else {
    seq <- labels
    if (is.numeric(seq))
      seq <- as.character(seq)
  }
  
  secondfirst <- function(x, sep)
    paste0(x[2], sep, x[1])
  #
  trts.list <- lapply(compsplit(rownames(net.iv$Cov.common), sep.trts),
                      secondfirst, sep = sep.trts)
  #
  rn <- rep_len("", nrow(H))
  #
  for (i in seq_along(trts.list))
    rn[i] <- trts.list[[i]]
  #
  rownames(H) <- rn
  
  
  res <- list(studlab = studlab,
              treat1 = treat1,
              treat2 = treat2,
              ##
              TE = data$.TE[!data$.drop],
              seTE = data$.seTE[!data$.drop],
              seTE.adj = rep(NA, sum(!data$.drop)),
              seTE.adj.common = rep(NA, sum(!data$.drop)),
              seTE.adj.random = rep(NA, sum(!data$.drop)),
              ##
              design = designs(treat1, treat2, studlab)$design,
              ##
              event1 = event1,
              event2 = event2,
              n1 = n1,
              n2 = n2,
              ##
              k = n.studies,
              m = m,
              n = length(trts),
              d = d,
              ##
              trts = trts,
              k.trts = rowSums(A),
              n.trts = NA,
              events.trts = NA,
              ##
              studies = NA,
              narms = NA,
              ##
              designs = designs,
              ##
              TE.common = TE.common,
              seTE.common = seTE.common,
              lower.common = ci.f$lower,
              upper.common = ci.f$upper,
              statistic.common = ci.f$statistic,
              pval.common = ci.f$p,
              ##
              TE.random = if (is.lrp) TE.random else NAmatrix,
              seTE.random = if (is.lrp) seTE.random else NAmatrix,
              lower.random = if (is.lrp) ci.r$lower else NAmatrix,
              upper.random = if (is.lrp) ci.r$upper else NAmatrix,
              statistic.random = if (is.lrp) ci.r$statistic else NAmatrix,
              pval.random = if (is.lrp) ci.r$p else NAmatrix,
              ##
              seTE.predict = NAmatrix,
              lower.predict = NAmatrix,
              upper.predict = NAmatrix,
              method.predict = "V",
              ##
              prop.direct.common = NA,
              prop.direct.random = NA,
              ##
              TE.direct.common = TE.direct.common,
              seTE.direct.common = seTE.direct.common,
              lower.direct.common = ci.direct.c$lower,
              upper.direct.common = ci.direct.c$upper,
              statistic.direct.common = ci.direct.c$statistic,
              pval.direct.common = ci.direct.c$p,
              ##
              TE.direct.random = if (is.lrp) TE.direct.random else NAmatrix,
              seTE.direct.random = if (is.lrp) seTE.direct.random else NAmatrix,
              lower.direct.random = if (is.lrp) ci.direct.c$lower else NAmatrix,
              upper.direct.random = if (is.lrp) ci.direct.c$upper else NAmatrix,
              statistic.direct.random =
                if (is.lrp) ci.direct.c$statistic else NAmatrix,
              pval.direct.random = if (is.lrp) ci.direct.c$p else NAmatrix,
              #
              Q.direct = Q.direct,
              tau.direct = sqrt(tau2.direct),
              tau2.direct = tau2.direct,
              I2.direct = I2.direct,
              #
              phi = if (is.lrp) phi else NA,
              phi.direct = if (is.lrp) phi.direct else NAmatrix,
              #
              TE.indirect.common = NA,
              seTE.indirect.common = NA,
              lower.indirect.common = NA,
              upper.indirect.common = NA,
              statistic.indirect.common = NA,
              pval.indirect.common = NA,
              ##
              TE.indirect.random = NA,
              seTE.indirect.random = NA,
              lower.indirect.random = NA,
              upper.indirect.random = NA,
              statistic.indirect.random = NA,
              pval.indirect.random = NA,
              ##
              Q = NA,
              df.Q = NA,
              pval.Q = NA,
              I2 = NA, lower.I2 = NA, upper.I2 = NA,
              tau = NA,
              ##
              Q.heterogeneity = NA,
              df.Q.heterogeneity = NA,
              pval.Q.heterogeneity = NA,
              Q.inconsistency = Q,
              df.Q.inconsistency = df.Q,
              pval.Q.inconsistency = pval.Q,
              ##
              Q.decomp = NA,
              ##
              A.matrix = A,
              H.matrix.common = H,
              H.matrix.random = NA,
              ##
              n.matrix = NA,
              events.matrix = NA,
              ##
              P.common = NA,
              P.random = NA,
              ##
              Cov.common = net.iv$Cov.common,
              Cov.random = net.iv$Cov.random,
              ##
              sm = sm,
              method = method,
              ##
              incr = incr,
              method.incr = method.incr,
              allstudies = allstudies,
              cc.pooled = cc.pooled,
              ##
              level = level,
              level.ma = level.ma,
              common = common,
              random = random,
              ##
              prediction = prediction,
              level.predict = level.predict,
              ##
              reference.group = reference.group,
              baseline.reference = baseline.reference,
              small.values = small.values,
              all.treatments = all.treatments,
              seq = seq,
              ##
              method.tau = "DL",
              tau.preset = tau.preset,
              ##
              tol.multiarm = tol.multiarm,
              tol.multiarm.se = tol.multiarm.se,
              details.chkmultiarm = details.chkmultiarm,
              details.chkdata = details.chkdata,
              ##
              sep.trts = sep.trts,
              nchar.trts = nchar.trts,
              ##
              func.inverse = deparse(substitute(func.inverse)),
              #
              overall.hetstat = overall.hetstat,
              backtransf = backtransf,
              #
              title = title,
              ##
              data = if (keepdata) data else NULL,
              data.design = dat.design,
              ##
              warn = warn,
              call = match.call(),
              version = packageDescription("netmeta")$Version,
              #
              addincr = addincr,
              allincr = allincr
              )
  ##
  class(res) <- c("netmetabin", "netmeta")
  ##
  ## Add results for indirect treatment estimates
  ##
  n <- res$n
  ##
  res$prop.direct.common <-
    netmeasures(res, random = FALSE, warn = warn)$proportion
  if (is.logical(res$prop.direct.common))
    res$prop.direct.common <- as.numeric(res$prop.direct.common)
  #
  if (is.lrp) {
    res$prop.direct.random <-
      netmeasures(res, random = TRUE, warn = warn)$proportion
    if (is.logical(res$prop.direct.random))
      res$prop.direct.random <- as.numeric(res$prop.direct.random)
    #
    res$Cov.common <- matrix(NA, nrow = length(res$prop.direct.common),
                             ncol = length(res$prop.direct.common))
    res$Cov.random <- matrix(NA, nrow = length(res$prop.direct.random),
                             ncol = length(res$prop.direct.random))
    #
    rownames(res$Cov.common) <- colnames(res$Cov.common) <-
      names(res$prop.direct.common)
    rownames(res$Cov.random) <- colnames(res$Cov.random) <-
      names(res$prop.direct.random)
  }
  #
  P.common <- P.random <- matrix(NA, n, n)
  colnames(P.common) <- rownames(P.common) <-
    colnames(P.random) <- rownames(P.random) <- trts
  ##
  if (n == 2) {
    ##
    ## For two treatments only direct evidence is available
    ##
    res$prop.direct.common <- 1
    names(res$prop.direct.common) <- paste(labels, collapse = sep.trts)
    ##
    sel <- row(P.common) != col(P.common)
    P.common[sel] <- 1
    #
    if (is.lrp) {
      res$prop.direct.random <- 1
      names(res$prop.direct.random) <- paste(labels, collapse = sep.trts)
      ##
      sel <- row(P.random) != col(P.random)
      P.random[sel] <- 1
    }
  }
  else {
    k <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        k <- k + 1
        P.common[i, j] <- P.common[j, i] <- res$prop.direct.common[k]
        #
        if (is.lrp)
          P.random[i, j] <- P.random[j, i] <- res$prop.direct.random[k]
      }
    }
  }
  ##
  ## Set direct evidence estimates to 0 if only indirect evidence is available
  ## (otherwise indirect estimates would be NA as direct estimates are NA)
  ##
  TE.direct.common <- res$TE.direct.common
  ##
  TE.direct.common[is_zero(P.common)] <- 0
  ##
  ## Indirect estimate is NA if only direct evidence is available
  ##
  res$P.common <- P.common
  ##
  P.common[is_zero(P.common - 1)] <- NA
  P.common[P.common > 1] <- NA
  #
  if (is.lrp) {
    TE.direct.random <- res$TE.direct.random
    TE.direct.random[is_zero(P.random)] <- 0
    #
    # Indirect estimate is NA if only direct evidence is available
    #
    res$P.random <- P.random
    #
    P.random[is_zero(P.random - 1)] <- NA
    P.random[P.random > 1] <- NA
  }
  ##
  ## Common effects model
  ##
  ci.indirect.c <-
    ci((res$TE.common - P.common * TE.direct.common) / (1 - P.common),
       sqrt(res$seTE.common^2 / (1 - P.common)), level = level)
  ##
  res$TE.indirect.common   <- ci.indirect.c$TE
  res$seTE.indirect.common <- ci.indirect.c$seTE
  ##
  res$lower.indirect.common <- ci.indirect.c$lower
  res$upper.indirect.common <- ci.indirect.c$upper
  ##
  res$statistic.indirect.common <- ci.indirect.c$statistic
  res$pval.indirect.common <- ci.indirect.c$p
  ##
  ## Random effects model
  ##
  if (!is.lrp) {
    res$prop.direct.random <- res$prop.direct.common
    res$prop.direct.random[!is.na(res$prop.direct.random)] <- NA
    res$P.random <- P.random
    ##
    res$TE.indirect.random <- res$seTE.indirect.random <-
      res$lower.indirect.random <- res$upper.indirect.random <-
      res$statistic.indirect.random <- res$pval.indirect.random <-
      NAmatrix
  }
  else {
    ci.indirect.r <-
      ci((res$TE.random - P.random * TE.direct.random) / (1 - P.random),
         sqrt(res$seTE.random^2 / (1 - P.random)), level = level)
    ##
    res$TE.indirect.random   <- ci.indirect.r$TE
    res$seTE.indirect.random <- ci.indirect.r$seTE
    ##
    res$lower.indirect.random <- ci.indirect.r$lower
    res$upper.indirect.random <- ci.indirect.r$upper
    ##
    res$statistic.indirect.random <- ci.indirect.r$statistic
    res$pval.indirect.random <- ci.indirect.r$p
    #
    res$P.random <- P.random
  }
  
  
  ##
  ## Study overview
  ##
  p0 <- prepare(rep(1, nrow(dat.wide)),
                rep(1, nrow(dat.wide)),
                dat.wide$treat1,
                dat.wide$treat2,
                dat.wide$studlab,
                func.inverse = func.inverse)
  ##
  tdata <- data.frame(studies = p0$studlab, narms = p0$narms,
                      stringsAsFactors = FALSE)
  ##
  tdata <- tdata[!duplicated(tdata[, c("studies", "narms")]), , drop = FALSE]
  res$studies <- tdata$studies
  res$narms <- tdata$narms
  ##
  if (all(suppressWarnings(!is.na(as.numeric(res$studies))))) {
    o <- order(as.numeric(res$studies))
    res$studies <- res$studies[o]
    res$narms <- res$narms[o]
  }
  ##
  res$data <- merge(res$data,
                    data.frame(.studlab = res$studies,
                               .narms = res$narms),
                    by = ".studlab", all.x = TRUE)
  
  
  res$data <- res$data[order(res$data$.order), ]
  res$data$.order <- NULL
  rownames(res$data) <- seq_len(nrow(res$data))
  
  
  res$n.matrix <- netmatrix(res, n1 + n2, func = "sum")
  ##
  dat.n <- data.frame(studlab = c(studlab, studlab),
                      treat = c(treat1, treat2),
                      n = c(n1, n2))
  dat.n <- dat.n[!duplicated(dat.n[, c("studlab", "treat")]), ]
  dat.n <- by(dat.n$n, dat.n$treat, sum, na.rm = TRUE)
  res$n.trts <- as.vector(dat.n[trts])
  names(res$n.trts) <- trts
  ##
  res$events.matrix <- netmatrix(res, event1 + event2, func = "sum")
  ##
  dat.e <- data.frame(studlab = c(studlab, studlab),
                      treat = c(treat1, treat2),
                      n = c(event1, event2))
  dat.e <- dat.e[!duplicated(dat.e[, c("studlab", "treat")]), ]
  dat.e <- by(dat.e$n, dat.e$treat, sum, na.rm = TRUE)
  res$events.trts <- as.vector(dat.e[trts])
  names(res$events.trts) <- trts
  ##
  ## Backward compatibility
  ##
  res$fixed <- res$common
  res$comb.fixed <- res$common
  res$comb.random <- res$random
  ##
  res$seTE.adj.fixed <- res$seTE.adj.common
  ##
  res$TE.fixed <- res$TE.common
  res$seTE.fixed <- res$seTE.common
  res$lower.fixed <- res$lower.common
  res$upper.fixed <- res$upper.common
  res$statistic.fixed <- res$statistic.common
  res$pval.fixed <- res$pval.common
  ##
  res$prop.direct.fixed <- res$prop.direct.common
  ##
  res$TE.direct.fixed <- res$TE.direct.common
  res$seTE.direct.fixed <- res$seTE.direct.common
  res$lower.direct.fixed <- res$lower.direct.common
  res$upper.direct.fixed <- res$upper.direct.common
  res$statistic.direct.fixed <- res$statistic.direct.common
  res$pval.direct.fixed <- res$pval.direct.common
  ##
  res$TE.indirect.fixed <- res$TE.indirect.common
  res$seTE.indirect.fixed <- res$seTE.indirect.common
  res$lower.indirect.fixed <- res$lower.indirect.common
  res$upper.indirect.fixed <- res$upper.indirect.common
  res$statistic.indirect.fixed <- res$statistic.indirect.common
  res$pval.indirect.fixed <- res$pval.indirect.common
  ##
  res$P.fixed <- res$P.common
  res$Cov.fixed <- res$Cov.common               
  
  res
}
