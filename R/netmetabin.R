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
#'   \code{"MH"}, or \code{"NCH"}, can be abbreviated.
#' @param cc.pooled A logical indicating whether \code{incr} should be
#'   used as a continuity correction, when calculating the network
#'   meta-analysis estimates.
#' @param incr A numerical value which is added to each cell count,
#'   i.e., to the numbers of events and non-events, of all treatment
#'   arms in studies with zero events or non-events in any of the
#'   treatment arms ("continuity correction").
#' @param allincr A logical indicating whether \code{incr} should be
#'   added to each cell count of all studies if a continuity
#'   correction was used for at least one study (only considered if
#'   \code{method = "Inverse"}). If FALSE (default), \code{incr} is
#'   used as continuity correction only for studies with zero events
#'   or zero non-events in any of the treatment arms.
#' @param addincr A logical indicating whether \code{incr} should be
#'   added to each cell count of all studies, irrespective of zero
#'   cell counts (only considered if \code{method = "Inverse"}).
#' @param allstudies A logical indicating whether studies with zero
#'   events or non-events in all treatment arms should be included in
#'   an inverse variance meta-analysis (applies only if \code{method =
#'   "Inverse"} and \code{sm} is equal to either \code{"RR"} or
#'   \code{"OR"}).
#' @param level The level used to calculate confidence intervals for
#'   individual studies.
#' @param level.comb The level used to calculate confidence intervals
#'   for pooled estimates.
#' @param comb.fixed A logical indicating whether a fixed effects
#'   (common effects) network meta-analysis should be conducted.
#' @param comb.random A logical indicating whether a random effects
#'   network meta-analysis should be conducted.
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
#' @param all.treatments A logical or \code{"NULL"}. If \code{TRUE},
#'   matrices with all treatment effects, and confidence limits will
#'   be printed.
#' @param seq A character or numerical vector specifying the sequence
#'   of treatments in printouts.
#' @param tau.preset An optional value for manually setting the
#'   square-root of the between-study variance \eqn{\tau^2} (only
#'   considered if \code{method = "Inverse"}).
#' @param tol.multiarm A numeric for the tolerance for consistency of
#'   treatment estimates and corresponding variances in multi-arm
#'   studies which are consistent by design (only considered if
#'   \code{method = "Inverse"}).
#' @param details.chkmultiarm A logical indicating whether treatment
#'   estimates and / or variances of multi-arm studies with
#'   inconsistent results or negative multi-arm variances should be
#'   printed (only considered if \code{method = "Inverse"}).
#' @param sep.trts A character used in comparison names as separator
#'   between treatment labels.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param title Title of meta-analysis / systematic review.
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in netmeta object.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if studies are excluded from meta-analysis due to zero
#'   standard errors).
#' 
#' @details
#' This function implements three models for the network meta-analysis
#' of binary data:
#' \itemize{
#' \item The Mantel-Haenszel network meta-analysis model, as described
#'   in Efthimiou et al. (2018) (\code{method = "MH"});
#' \item a network meta-analysis model using the non-central
#'   hypergeometric distribution with the Breslow approximation, as
#'   described in Stijnen et al. (2010) (\code{method = "NCH"});
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
#' arm separately, function \code{\link{pairwise}} can be used to
#' transform the data to \emph{contrast-based} format (see help page
#' of function \code{\link{pairwise}}).
#' 
#' Note, all pairwise comparisons must be provided for a multi-arm
#' study. Consider a multi-arm study of \emph{p} treatments with known
#' variances. For this study, the number of events and observations
#' must be provided for each treatment, for each of \emph{p}(\emph{p}
#' - 1) / 2 possible comparisons in separate lines in the data. For
#' instance, a three-arm study contributes three pairwise comparisons,
#' a four-arm study even six pairwise comparisons. Function
#' \code{\link{pairwise}} automatically calculates all pairwise
#' comparisons for multi-arm studies.
#' 
#' For \code{method = "Inverse"}, both fixed effects and random
#' effects models are calculated regardless of values choosen for
#' arguments \code{comb.fixed} and \code{comb.random}. Accordingly,
#' the network estimates for the random effects model can be extracted
#' from component \code{TE.random} of an object of class
#' \code{"netmeta"} even if argument \code{comb.random = FALSE}.
#' However, all functions in R package \bold{netmeta} will adequately
#' consider the values for \code{comb.fixed} and
#' \code{comb.random}. E.g. function
#' \code{\link{print.summary.netmeta}} will not print results for the
#' random effects model if \code{comb.random = FALSE}.
#' 
#' For \code{method = "MH"} and \code{method = "NCH"}, only a fixed
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
#' used in the covariance matrices \code{Cov.fixed} and
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
#' \item{TE.fixed, seTE.fixed}{\emph{n}x\emph{n} matrix with estimated
#'   overall treatment effects and standard errors for fixed effects
#'   model.}
#' \item{lower.fixed, upper.fixed}{\emph{n}x\emph{n} matrices with
#'   lower and upper confidence interval limits for fixed effects
#'   model.}
#' \item{zval.fixed, pval.fixed}{\emph{n}x\emph{n} matrices with
#'   z-value and p-value for test of overall treatment effect under
#'   fixed effects model.}
#' \item{TE.random, seTE.random}{\emph{n}x\emph{n} matrix with
#'   estimated overall treatment effects and standard errors for
#'   random effects model (only available if \code{method =
#'   "Inverse"}).}
#' \item{lower.random, upper.random}{\emph{n}x\emph{n} matrices with
#'   lower and upper confidence interval limits for random effects
#'   model (only available if \code{method = "Inverse"}).}
#' \item{zval.random, pval.random}{\emph{n}x\emph{n} matrices with
#'   z-value and p-value for test of overall treatment effect under
#'   random effects model (only available if \code{method =
#'   "Inverse"}).}
#' \item{TE.direct.fixed, seTE.direct.fixed}{\emph{n}x\emph{n} matrix
#'   with estimated treatment effects and standard errors from direct
#'   evidence under fixed effects model.}
#' \item{lower.direct.fixed, upper.direct.fixed}{\emph{n}x\emph{n}
#'   matrices with lower and upper confidence interval limits from
#'   direct evidence under fixed effects model.}
#' \item{zval.direct.fixed, pval.direct.fixed}{\emph{n}x\emph{n}
#'   matrices with z-value and p-value for test of overall treatment
#'   effect from direct evidence under fixed effects model.}
#' \item{TE.direct.random, seTE.direct.random}{\emph{n}x\emph{n}
#'   matrix with estimated treatment effects and standard errors from
#'   direct evidence under random effects model (only available if
#'   \code{method = "Inverse"}).}
#' \item{lower.direct.random, upper.direct.random}{\emph{n}x\emph{n}
#'   matrices with lower and upper confidence interval limits from
#'   direct evidence under random effects model (only available if
#'   \code{method = "Inverse"}).}
#' \item{zval.direct.random, pval.direct.random}{\emph{n}x\emph{n}
#'   matrices with z-value and p-value for test of overall treatment
#'   effect from direct evidence under random effects model (only
#'   available if \code{method = "Inverse"}).}
#' \item{Q}{Overall heterogeneity / inconsistency statistic.}
#' \item{df.Q}{Degrees of freedom for test of heterogeneity /
#'   inconsistency.}
#' \item{pval.Q}{P-value for test of heterogeneity / inconsistency.}
#' \item{I2}{I-squared (only available if \code{method = "Inverse"}).}
#' \item{tau}{Square-root of between-study variance (only available if
#'   \code{method = "Inverse"}).}
#' \item{A.matrix}{Adjacency matrix (\emph{n}x\emph{n}).}
#' \item{H.matrix}{Hat matrix (\emph{m}x\emph{m})}
#' \item{n.matrix}{\emph{n}x\emph{n} matrix with number of
#'   observations in direct comparisons.}
#' \item{events.matrix}{\emph{n}x\emph{n} matrix with number of events
#'   in direct comparisons.}
#' \item{sm, method, level, level.comb}{As defined above.}
#' \item{incr, allincr, addincr, allstudies, cc.pooled}{As defined
#'   above.}
#' \item{comb.fixed, comb.random}{As defined above.}
#' \item{prediction, level.predict}{As defined above.}
#' \item{reference.group, baseline.reference, all.treatments}{As
#'   defined above.}
#' \item{seq, tau.preset, tol.multiarm, details.chkmultiarm}{As
#'   defined above.}
#' \item{sep.trts, nchar.trts}{As defined above.}
#' \item{backtransf, title, warn}{As defined above.}
#' \item{data}{Data set (in contrast-based format).}
#' \item{data.design}{List with data in arm-based format (each list
#'   element corresponds to a single design).}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Orestis Efthimiou \email{oremiou@@gmail.com}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{pairwise}}, \code{\link{netmeta}}
#' 
#' @references
#' Efthimiou O, Rücker G, Schwarzer G, Higgins J, Egger M, Salanti G
#' (2018):
#' A Mantel-Haenszel model for network meta-analysis of rare events.
#' \emph{Manuscript submitted for publication}
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
#' data(Dong2013)
#' 
#' # Only consider first ten studies (to reduce runtime of example)
#' #
#' first10 <- subset(Dong2013, id <= 10)
#' 
#' # Transform data from long arm-based format to contrast-based
#' # format. Argument 'sm' has to be used for odds ratio as summary
#' # measure; by default the risk ratio is used in the metabin
#' # function called internally.
#' #
#' p1 <- pairwise(treatment, death, randomized, studlab = id,
#'                data = first10, sm = "OR")
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
#'                data = Dong2013, sm = "OR")
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
#' data(Gurusamy2011)
#' 
#' p3 <- pairwise(treatment, death, n, studlab = study,
#'                data = Gurusamy2011, sm = "OR")
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
                       incr = gs("incr"),
                       allincr = gs("allincr"), addincr = gs("addincr"),
                       allstudies = gs("allstudies"),
                       level = gs("level"),
                       level.comb = gs("level.comb"),
                       comb.fixed = gs("comb.fixed"),
                       comb.random = method == "Inverse" &
                         (gs("comb.random") | !is.null(tau.preset)),
                       ##
                       prediction = FALSE,
                       level.predict = gs("level.predict"),
                       ##
                       reference.group = "",
                       baseline.reference = TRUE,
                       all.treatments = NULL,
                       seq = NULL,
                       ##
                       tau.preset = NULL,
                       ##
                       tol.multiarm = 0.0005,
                       details.chkmultiarm = FALSE,
                       ##
                       sep.trts = ":",
                       nchar.trts = 666,
                       ##
                       backtransf = gs("backtransf"),
                       ##
                       title = "",
                       keepdata = gs("keepdata"),
                       warn = TRUE) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkchar <- meta:::chkchar
  chklevel <- meta:::chklevel
  chklogical <- meta:::chklogical
  chknull <- meta:::chknull
  chknumeric <- meta:::chknumeric
  setchar <- meta:::setchar
  pvalQ <- meta:::pvalQ
  ##
  modtext <- paste0("must be equal to 'Inverse' (classic network meta-analysis), ",
                    "'MH' (Mantel-Haenszel, the default) or ",
                    "'NCH' (common-effects non-central hypergeometric).")
  method <- meta:::setchar(method, c("Inverse", "MH", "NCH"), modtext)
  ##
  chklogical(allincr)
  chklogical(addincr)
  chklogical(allstudies)
  chklogical(cc.pooled)
  ##
  chklevel(level)
  chklevel(level.comb)
  chklevel(level.predict)
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
  ##
  if (method != "Inverse" & !comb.fixed) {
    warning("Argument 'comb.fixed' set to TRUE for Mantel-Haenszel ",
            "method and non-central hypergeometric distribution.",
            call. = FALSE)
    comb.fixed <- TRUE
  }
  ##
  if (method != "Inverse" & comb.random) {
    warning("Argument 'comb.random' set to FALSE for Mantel-Haenszel ",
            "method and non-central hypergeometric distribution.",
            call. = FALSE)
    comb.random <- FALSE
  }
  ##
  if (method != "Inverse" & prediction) {
    warning("Argument 'prediction' set to FALSE for Mantel-Haenszel ",
            "method and non-central hypergeometric distribution.",
            call. = FALSE)
    prediction <- FALSE
  }
  ##
  chklogical(baseline.reference)
  ##
  if (!is.null(all.treatments))
    chklogical(all.treatments)
  ##
  if (!is.null(tau.preset))
    chknumeric(tau.preset, min = 0, single = TRUE)
  ##
  chknumeric(tol.multiarm, min = 0, single = TRUE)
  chklogical(details.chkmultiarm)
  ##
  missing.sep.trts <- missing(sep.trts)
  chkchar(sep.trts)
  chknumeric(nchar.trts, min = 1, single = TRUE)
  ##
  chklogical(backtransf)
  ##
  chkchar(title)
  chklogical(keepdata)
  chklogical(warn)
  ##
  ## Check value for reference group
  ##
  if (is.null(all.treatments))
    if (reference.group == "")
      all.treatments <- TRUE
    else
      all.treatments <- FALSE
  ##
  chklogical(baseline.reference)


  ##
  ##
  ## (2) Read data
  ##
  ##
  nulldata <- is.null(data)
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch 'event1', 'event2', 'n1', 'n2', 'treat1', 'treat2', and
  ## 'studlab' from data:
  ##
  event1 <- eval(mf[[match("event1", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  if (inherits(event1, "pairwise")) {
    is.pairwise <- TRUE
    ##
    if (missing(sm) & method == "Inverse")
      sm <- attr(event1, "sm")
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
  }
  else {
    is.pairwise <- FALSE
    ##
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
    n1 <- eval(mf[[match("n1", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
    ##
    event2 <- eval(mf[[match("event2", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
    ##
    n2 <- eval(mf[[match("n2", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
    ##
    treat1 <- eval(mf[[match("treat1", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
    ##
    treat2 <- eval(mf[[match("treat2", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
    ##
    studlab <- eval(mf[[match("studlab", names(mf))]],
                    data, enclos = sys.frame(sys.parent()))
  }
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
  ## Keep original order of studies
  ##
  .order <- seq_along(studlab)
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  missing.subset <- is.null(subset)
  
  
  ##
  ##
  ## (2b) Store complete dataset in list object data
  ##      (if argument keepdata is TRUE)
  ##
  ##
  if (keepdata) {
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
  }
  
  
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
    if (warn)
      warning("Note, treatments within a comparison have been ",
              "re-sorted in increasing order.",
              call. = FALSE)
    ##
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
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)
      abs(x - round(x)) < tol
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
  ##
  ## Catch 'incr' from data:
  ##
  if (!missing(incr))
    incr <- eval(mf[[match("incr", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknumeric(incr, min = 0)
  ##
  if (method != "Inverse" & length(incr) > 1) {
    warning("Argument 'incr' must be a single value for ",
            "Mantel-Haenszel and common-effects non-central ",
            "hypergeometric method. Set to zero (default).",
            call. = FALSE)
    incr <- 0
  }
  ##
  if (all(incr == 0) & method != "Inverse" & cc.pooled == TRUE)
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
  get.designs <- function(x) {
    ##
    res <- nma.krahn(netmeta(pairwise(studlab = x$studlab,
                                      treat = x$treat,
                                      event = x$event,
                                      n = x$n,
                                      sm = "RD"),
                             warn = FALSE)
                     )$studies[, c("studlab", "design")]
    ##
    if (is.null(res)) {
      net1 <- netmeta(pairwise(studlab = x$studlab,
                               treat = x$treat,
                               event = x$event,
                               n = x$n,
                               sm = "RD"),
                      warn = FALSE)
      ##
      res <- data.frame(studlab = net1$studlab, design = net1$designs)
    }
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
  dat.long <- merge(dat.long, tdat.design, by = "studlab")
  dat.wide <- merge(dat.wide, tdat.design, by = "studlab")
  ##
  dat.long <- dat.long[order(dat.long$.order), ]
  dat.wide <- dat.wide[order(dat.wide$.order), ]
  ##
  names(tdat.design) <- c("studlab", ".design")
  ##
  data <- merge(data, tdat.design, by = "studlab", all.x = TRUE)
  data <- data[order(data$.order), ]
  ##
  rm(tdat.design)
  
  
  ##
  ##
  ## (6) Stage I: setting up the data (Efthimiou et al., 2018)
  ##
  ##
  ## Step i. Remove all-zero studies (only MH and NCH methods)
  ##
  if (method != "Inverse") {
    n.events <- with(dat.long, tapply(event, studlab, sum))
    ##
    if (any(n.events == 0)) {
      allzero <- n.events == 0
      ##
      if (warn)
        if (sum(allzero) == 1)
          warning("Study '", names(n.events)[allzero],
                  "' without any events excluded from network meta-analysis.",
                  call. = FALSE)
        else
          warning("Studies without any events excluded ",
                  "from network meta-analysis: ",
                  paste(paste0("'", names(n.events)[allzero], "'"),
                        collapse = " - "),
                  call. = FALSE)
      ##
      dat.long <- dat.long[dat.long$studlab %in% names(n.events)[!allzero], ,
                           drop = FALSE]
      dat.wide <- dat.wide[dat.wide$studlab %in% names(n.events)[!allzero], ,
                           drop = FALSE]
      ##
      data$.drop <- data$.drop | data$studlab %in% names(n.events)[allzero]
      ##
      rm(allzero)
    }
    ##
    rm(n.events)
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
  ## Step iii. Drop treatment arms without events from individual
  ##           designs (argument 'cc.pooled' is FALSE) or add
  ##           increment if argument 'cc.pooled' is TRUE and argument
  ##           'incr' is larger than zero
  ##
  if (method != "Inverse") {
    ##
    d.events <- with(dat.long, tapply(event, list(design, treat), sum))
    ##
    if (any(d.events == 0, na.rm = TRUE)) {
      zero <- d.events == 0
      zerocells <- as.data.frame(which(zero, arr.ind = TRUE),
                                 stringsAsFactors = FALSE)
      ##
      zerocells$design <- rownames(zero)[zerocells$row]
      zerocells$treat <- colnames(zero)[zerocells$col]
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
      if (!cc.pooled) {
        if (warn)
          if (sum(zero, na.rm = TRUE) == 1)
            warning("Treatment arm '", zerocells$treat,
                    "' without events in design '",
                    zerocells$design, "' excluded from network meta-analysis.",
                    call. = FALSE)
          else
            warning("Treatment arms without events in a design excluded ",
                    "from network meta-analysis:\n    ",
                    paste0("'",
                           paste0(paste0(zerocells$treat, " in "),
                                  zerocells$design),
                           "'",
                           collapse = " - "),
                    call. = FALSE)
        ##
        dat.long <- dat.long[!zero.long, , drop = FALSE]
        dat.wide <- dat.wide[!zero.wide, , drop = FALSE]
        data$.drop <- data$.drop | zero.data
      }
      else {
        dat.long$incr[zero.long] <- incr
        ##
        dat.long$event[zero.long] <- dat.long$event[zero.long] + incr
        dat.long$non.event[zero.long] <- dat.long$non.event[zero.long] + incr
        dat.long$n[zero.long] <- dat.long$n[zero.long] + 2 * incr
        ##
        dat.wide$incr[zero.wide] <- incr
        ##
        data$.incr[zero.data] <- incr
      }
      ##
      rm(zero, zerocells, zero.long, zero.wide, zero.data)
    }
    ##
    rm(d.events)
  }
  ##
  ## Step iv. Remove designs with single treatment arm from dataset
  ##
  if (method != "Inverse") {
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
  Z <- matrix(NA, nrow = n.treat, ncol = n.treat)
  rownames(Z) <- colnames(Z) <- trts
  ##
  TE.fixed <- seTE.fixed <- Z
  TE.direct.fixed <- seTE.direct.fixed <- Z
  ##
  rm(Z)


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
                   allincr = allincr, addincr = addincr,
                   allstudies = allstudies.iv)
  ##
  net.iv <- netmeta(p.iv,
                    level = level, level.comb = level.comb,
                    comb.fixed = comb.fixed, comb.random = comb.random,
                    prediction = prediction, level.predict = level.predict,
                    reference.group = reference.group,
                    baseline.reference = baseline.reference,
                    all.treatments = all.treatments,
                    seq = seq,
                    tau.preset = tau.preset,
                    tol.multiarm = tol.multiarm,
                    details.chkmultiarm = details.chkmultiarm,
                    sep.trts = sep.trts,
                    nchar.trts = nchar.trts,
                    backtransf = backtransf,
                    title = title,
                    keepdata = keepdata,
                    warn = warn.iv)
  ##
  if (method == "Inverse")
    return(net.iv)
  else {
    net.iv$Cov.fixed[!is.na(net.iv$Cov.fixed)] <- NA
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
      ## Calculate total patients per studies
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
  ##
  ## Define H matrix
  ##
  H <- matrix(0,
              nrow = n.treat * ((n.treat - 1) / 2),
              ncol = n.treat - 1)
  ##
  diag(H) <- 1
  ##
  if (n.treat > 2) {
    t1 <- c()
    t2 <- c()
    for (i in 2:(n.treat - 1))
      for (j in (i + 1):(n.treat)) {
        t1 <- rbind(t1, i)
        t2 <- rbind(t2, j)
      }
    ##
    h1 <- matrix(c(t1, t2), nrow = nrow(t1))
    ##
    for (i in 1:((n.treat - 1) * (n.treat - 2) / 2))
      for (j in 1:(n.treat - 1))
        H[i + n.treat - 1, j] <- -(h1[i, 1] == j + 1) + (h1[i, 2] == j + 1)
  }
  ##
  ## Fixed effects matrix
  ##
  d.hat <- H %*% TE.basic
  ##
  TE.fixed[lower.tri(TE.fixed, diag = FALSE)] <-  d.hat
  TE.fixed <- t(TE.fixed)
  TE.fixed[lower.tri(TE.fixed, diag = FALSE)] <- -d.hat
  diag(TE.fixed) <- 0
  ##
  ## Matrix with standard errors
  ##
  if (method == "MH")
    cov.d.hat <- H %*% solve(t(X) %*% W %*% X) %*% t(H)
  else
    cov.d.hat <- H %*% W %*% t(H)
  ##
  seTE.fixed[lower.tri(seTE.fixed, diag = FALSE)] <- sqrt(diag(cov.d.hat))
  seTE.fixed <- t(seTE.fixed)
  seTE.fixed[lower.tri(seTE.fixed, diag = FALSE)] <- sqrt(diag(cov.d.hat))
  diag(seTE.fixed) <- 0
  ##
  ## Confidence intervals
  ##
  ci.f <- ci(TE.fixed, seTE.fixed, level = level.comb)
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
    m.i <- metabin(event1, n1, event2, n2, subset = selstud,
                   sm = "OR",
                   incr = incr, allincr = allincr, addincr = addincr,
                   allstudies = allstudies, MH.exact = !cc.pooled)
    ##
    TE.i   <- m.i$TE.fixed
    seTE.i <- m.i$seTE.fixed
    ##
    TE.direct.fixed[sel.treat1, sel.treat2]   <- TE.i
    seTE.direct.fixed[sel.treat1, sel.treat2] <- seTE.i
    ##
    TE.direct.fixed[sel.treat2, sel.treat1]   <- -TE.i
    seTE.direct.fixed[sel.treat2, sel.treat1] <- seTE.i
  }
  ##
  rm(sel.treat1, sel.treat2, selstud, m.i, TE.i, seTE.i)
  ##
  ci.d <- meta::ci(TE.direct.fixed, seTE.direct.fixed, level = level.comb)
  
  
  labels <- sort(unique(c(treat1, treat2)))
  ##
  if (!is.null(seq))
    seq <- setseq(seq, labels)
  else {
    seq <- labels
    if (is.numeric(seq))
      seq <- as.character(seq)
  }
  
  
  res <- list(studlab = studlab,
              treat1 = treat1,
              treat2 = treat2,
              ##
              TE = data$TE[!data$.drop],
              seTE = data$seTE[!data$.drop],
              seTE.adj = rep(NA, sum(!data$.drop)),
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
              TE.fixed = TE.fixed,
              seTE.fixed = seTE.fixed,
              lower.fixed = ci.f$lower,
              upper.fixed = ci.f$upper,
              zval.fixed = ci.f$z,
              pval.fixed = ci.f$p,
              ##
              TE.random = NA,
              seTE.random = NA,
              lower.random = NA,
              upper.random = NA,
              zval.random = NA,
              pval.random = NA,
              ##
              prop.direct.fixed = NA,
              prop.direct.random = NA,
              ##
              TE.direct.fixed = TE.direct.fixed,
              seTE.direct.fixed = seTE.direct.fixed,
              lower.direct.fixed = ci.d$lower,
              upper.direct.fixed = ci.d$upper,
              zval.direct.fixed = ci.d$z,
              pval.direct.fixed = ci.d$p,
              ##
              TE.indirect.fixed = NA,
              seTE.indirect.fixed = NA,
              lower.indirect.fixed = NA,
              upper.indirect.fixed = NA,
              zval.indirect.fixed = NA,
              pval.indirect.fixed = NA,
              ##
              Q = Q,
              df.Q = df.Q,
              pval.Q = pval.Q,
              I2 = NA,
              tau = NA,
              ##
              Q.heterogeneity = NA,
              df.Q.heterogeneity = NA,
              pval.Q.heterogeneity = NA,
              Q.inconsistency = NA,
              df.Q.inconsistency = NA,
              pval.Q.inconsistency = NA,
              ##
              Q.decomp = NA,
              ##
              A.matrix = A,
              H.matrix = H,
              ##
              n.matrix = NA,
              events.matrix = NA,
              ##
              P.fixed = NA,
              P.random = NA,
              ##
              Cov.fixed = net.iv$Cov.fixed,
              Cov.random = net.iv$Cov.random,
              ##
              sm = sm,
              method = method,
              ##
              incr = incr,
              allincr = allincr,
              addincr = addincr,
              allstudies = allstudies,
              cc.pooled = cc.pooled,
              ##
              level = level,
              level.comb = level.comb,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              ##
              prediction = prediction,
              level.predict = level.predict,
              ##
              reference.group = reference.group,
              baseline.reference = baseline.reference,
              all.treatments = all.treatments,
              seq = seq,
              ##
              tau.preset = tau.preset,
              ##
              tol.multiarm = tol.multiarm,
              details.chkmultiarm = details.chkmultiarm,
              ##
              sep.trts = sep.trts,
              nchar.trts = nchar.trts,
              ##
              backtransf = backtransf,
              ##
              title = title,
              ##
              data = data,
              data.design = dat.design,
              ##
              warn = warn,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- c("netmetabin", "netmeta")
  
  
  ##
  ## Study overview
  ##
  p0 <- prepare(rep(1, nrow(dat.wide)),
                rep(1, nrow(dat.wide)),
                dat.wide$treat1,
                dat.wide$treat2,
                dat.wide$studlab)
  tdata <- data.frame(studies = p0$studlab, narms = p0$narms,
                      stringsAsFactors = FALSE)
  tdata <- unique(tdata[order(tdata$studies, tdata$narms), ])
  res$studies <- tdata$studies
  res$narms <- tdata$narms
  ##
  res$data <- merge(res$data,
                    data.frame(.studlab = res$studies,
                               .narms = res$narms),
                    by = ".studlab", all.x = TRUE)
  
  
  res$data <- res$data[order(res$data$.order), ]
  res$data$.order <- NULL
  rownames(res$data) <- seq_len(nrow(res$data))
  
  
  res$events.matrix <- netmatrix(res, event1 + event2, func = "sum")
  ##
  dat.e <- bySummary(c(event1, event2), c(treat1, treat2), long = FALSE)
  rownames(dat.e) <- dat.e$indices
  res$events.trts <- dat.e[trts, "sum"]
  names(res$events.trts) <- trts
  ##
  res$n.matrix <- netmatrix(res, n1 + n2, func = "sum")
  ##
  dat.n <- bySummary(c(n1, n2), c(treat1, treat2), long = FALSE)
  rownames(dat.n) <- dat.n$indices
  res$n.trts <- dat.n[trts, "sum"]
  names(res$n.trts) <- trts
  
  
  res
}
