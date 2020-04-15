#' Network meta-analysis using graph-theoretical method
#' 
#' @description
#' Network meta-analysis is a generalisation of pairwise meta-analysis
#' that compares all pairs of treatments within a number of treatments
#' for the same condition. The graph-theoretical approach for network
#' meta-analysis uses methods that were originally developed in
#' electrical network theory. It has been found to be equivalent to
#' the frequentist approach to network meta-analysis which is based on
#' weighted least squares regression (Rücker, 2012).
#' 
#' @param TE Estimate of treatment effect, i.e. difference between
#'   first and second treatment (e.g. log odds ratio, mean difference,
#'   or log hazard ratio).
#' @param seTE Standard error of treatment estimate.
#' @param treat1 Label/Number for first treatment.
#' @param treat2 Label/Number for second treatment.
#' @param studlab An optional - but important! - vector with study
#'   labels (see Details).
#' @param data An optional data frame containing the study
#'   information.
#' @param subset An optional vector specifying a subset of studies to
#'   be used.
#' @param sm A character string indicating underlying summary measure,
#'   e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
#'   \code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param level The level used to calculate confidence intervals for
#'   individual comparisons.
#' @param level.comb The level used to calculate confidence intervals
#'   for pooled estimates.
#' @param comb.fixed A logical indicating whether a fixed effects
#'   (common effects) network meta-analysis should be conducted.
#' @param comb.random A logical indicating whether a random effects
#'   network meta-analysis should be conducted.
#' @param prediction A logical indicating whether prediction intervals
#'   should be printed.
#' @param level.predict The level used to calculate prediction
#'   intervals for a new study.
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
#'   square-root of the between-study variance \eqn{\tau^2}.
#' @param tol.multiarm A numeric for the tolerance for consistency of
#'   treatment estimates in multi-arm studies which are consistent by
#'   design.
#' @param tol.multiarm.se A numeric for the tolerance for consistency
#'   of standard errors in multi-arm studies which are consistent by
#'   design.
#' @param details.chkmultiarm A logical indicating whether treatment
#'   estimates and / or variances of multi-arm studies with
#'   inconsistent results or negative multi-arm variances should be
#'   printed.
#' @param sep.trts A character used in comparison names as separator
#'   between treatment labels.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
#' @param n1 Number of observations in first treatment group.
#' @param n2 Number of observations in second treatment group.
#' @param event1 Number of events in first treatment group.
#' @param event2 Number of events in second treatment group.
#' @param title Title of meta-analysis / systematic review.
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in netmeta object.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if studies are excluded from meta-analysis due to zero
#'   standard errors).
#' 
#' @details
#' Network meta-analysis using R package \bold{netmeta} is described
#' in detail in Schwarzer et al. (2015), Chapter 8.
#' 
#' Let \emph{n} be the number of different treatments (nodes,
#' vertices) in a network and let \emph{m} be the number of existing
#' comparisons (edges) between the treatments. If there are only
#' two-arm studies, \emph{m} is the number of studies. Let TE and seTE
#' be the vectors of observed effects and their standard errors. Let W
#' be the \emph{m}x\emph{m} diagonal matrix that contains the inverse
#' variance 1 / seTE^2.
#' 
#' The given comparisons define the network structure. Therefrom an
#' \emph{m}x\emph{n} design matrix X (edge-vertex incidence matrix) is
#' formed; for more precise information, see Rücker (2012). Moreover,
#' the \emph{n}x\emph{n} Laplacian matrix L and its Moore-Penrose
#' pseudoinverse L+ are calculated (both matrices play an important
#' role in graph theory and electrical network theory). Using these
#' matrices, the variances based on both direct and indirect
#' comparisons can be estimated. Moreover, the hat matrix H can be
#' estimated by \strong{H = XL+X^tW = X(X^t W X)^+X^tW} and finally
#' consistent treatment effects can be estimated by applying the hat
#' matrix to the observed (potentially inconsistent) effects. H is a
#' projection matrix which maps the observed effects onto the
#' consistent (n-1)-dimensional subspace. This is the Aitken estimator
#' (Senn et al., 2013). As in pairwise meta-analysis, the Q statistic
#' measures the deviation from consistency. Q can be separated into
#' parts for each pairwise meta-analysis and a part for remaining
#' inconsistency between comparisons.
#' 
#' Often multi-arm studies are included in a network meta-analysis.
#' In multi-arm studies, the treatment effects on different
#' comparisons are not independent, but correlated. This is accounted
#' for by reweighting all comparisons of each multi-arm study. The
#' method is described in Rücker (2012) and Rücker and Schwarzer
#' (2014).
#' 
#' Comparisons belonging to multi-arm studies are identified by
#' identical study labels (argument \code{studlab}). It is therefore
#' important to use identical study labels for all comparisons
#' belonging to the same multi-arm study, e.g., study label
#' "Willms1999" for the three-arm study in the data example (Senn et
#' al., 2013). The function netmeta then automatically accounts for
#' within-study correlation by reweighting all comparisons of each
#' multi-arm study.
#' 
#' Data entry for this function is in \emph{contrast-based} format,
#' that is, data are given as contrasts (differences) between two
#' treatments (argument \code{TE}) with standard error (argument
#' \code{seTE}). In principle, meta-analysis functions from R package
#' \bold{meta}, e.g. \code{\link{metabin}} for binary outcomes or
#' \code{\link{metacont}} for continuous outcomes, can be used to
#' calculate treatment effects separately for each treatment
#' comparison which is a rather tedious enterprise. If data are
#' provided in \emph{arm-based} format, that is, data are given for
#' each treatment arm separately (e.g. number of events and
#' participants for binary outcomes), a much more convenient way to
#' transform data into contrast-based form is available. Function
#' \code{\link{pairwise}} can automatically transform data with binary
#' outcomes (using the \code{\link{metabin}} function from R package
#' \bold{meta}), continuous outcomes (\code{\link{metacont}}
#' function), incidence rates (\code{\link{metainc}} function), and
#' generic outcomes (\code{\link{metagen}} function). Additional
#' arguments of these functions can be provided, e.g., to calculate
#' Hedges' \emph{g} or Cohen's \emph{d} for continuous outcomes (see
#' help page of function \code{\link{pairwise}}).
#' 
#' Note, all pairwise comparisons must be provided for a multi-arm
#' study. Consider a multi-arm study of \emph{p} treatments with known
#' variances. For this study, treatment effects and standard errors
#' must be provided for each of \emph{p}(\emph{p} - 1) / 2 possible
#' comparisons. For instance, a three-arm study contributes three
#' pairwise comparisons, a four-arm study even six pairwise
#' comparisons. Function \code{\link{pairwise}} automatically
#' calculates all pairwise comparisons for multi-arm studies.
#' 
#' A simple random effects model assuming that a constant
#' heterogeneity variance is added to each comparison of the network
#' can be defined via a generalised methods of moments estimate of the
#' between-studies variance \eqn{\tau^2} (Jackson et al., 2012). This
#' is added to the observed sampling variance \code{seTE^2} of each
#' comparison in the network (before appropriate adjustment for
#' multi-arm studies). Then, as in standard pairwise meta-analysis,
#' the procedure is repeated with the resulting enlarged standard
#' errors.
#' 
#' For the random-effects model, the direct treatment estimates are
#' based on the common between-study variance \eqn{\tau^2} from the
#' network meta-analysis.
#' 
#' Internally, both fixed effects and random effects models are
#' calculated regardless of values choosen for arguments
#' \code{comb.fixed} and \code{comb.random}. Accordingly, the network
#' estimates for the random effects model can be extracted from
#' component \code{TE.random} of an object of class \code{"netmeta"}
#' even if argument \code{comb.random = FALSE}. However, all functions
#' in R package \bold{netmeta} will adequately consider the values for
#' \code{comb.fixed} and \code{comb.random}. E.g. function
#' \code{\link{print.summary.netmeta}} will not print results for the
#' random effects model if \code{comb.random = FALSE}.
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
#' An object of class \code{netmeta} with corresponding \code{print},
#' \code{summary}, \code{forest}, and \code{netrank} functions. The
#' object is a list containing the following components:
#' \item{studlab, treat1, treat2, TE, seTE}{As defined above.}
#' \item{seTE.adj}{Standard error of treatment estimate, adjusted for
#'   multi-arm studies.}
#' \item{n1, n2, event1, event2}{As defined above.}
#' \item{k}{Total number of studies.}
#' \item{m}{Total number of pairwise comparisons.}
#' \item{n}{Total number of treatments.}
#' \item{d}{Total number of designs (corresponding to the unique set
#'   of treatments compared within studies).}
#' \item{trts}{Treatments included in network meta-analysis.}
#' \item{k.trts}{Number of studies evaluating a treatment.}
#' \item{n.trts}{Number of observations receiving a treatment (if
#'   arguments \code{n1} and \code{n2} are provided).}
#' \item{events.trts}{Number of events observed for a treatment (if
#'   arguments \code{event1} and \code{event2} are provided).}
#' \item{multiarm}{Logical vector to identify pairwise comparisons
#'   from multi-arm studies.}
#' \item{n.arms}{Number of treatment arms in study providing pairwise
#'   comparison.}
#' \item{studies}{Vector with unique study labels.}
#' \item{narms}{Number of arms for each study.}
#' \item{designs}{Vector with unique designs present in the network. A
#'   design corresponds to the set of treatments compared within a
#'   study.}
#' \item{TE.nma.fixed, TE.nma.random}{A vector of length \emph{m} of
#'   consistent treatment effects estimated by network meta-analysis
#'   (nma) (fixed effects / random effects model).}
#' \item{seTE.nma.fixed, seTE.nma.random}{A vector of length \emph{m}
#'   of effective standard errors estimated by network meta-analysis
#'   (fixed effects / random effects model).}
#' \item{lower.nma.fixed, lower.nma.random}{A vector of length
#'   \emph{m} of lower confidence interval limits for consistent
#'   treatment effects estimated by network meta-analysis (fixed
#'   effects / random effects model).}
#' \item{upper.nma.fixed, upper.nma.random}{A vector of length
#'   \emph{m} of upper confidence interval limits for the consistent
#'   treatment effects estimated by network meta-analysis (fixed
#'   effects / random effects model).}
#' \item{zval.nma.fixed, zval.nma.random}{A vector of length \emph{m}
#'   of z-values for test of treatment effect for individual
#'   comparisons (fixed effects / random effects model).}
#' \item{pval.nma.fixed, pval.nma.random}{A vector of length \emph{m}
#'   of p-values for test of treatment effect for individual
#'   comparisons (fixed effects / random effects model).}
#' \item{leverage.fixed}{A vector of length \emph{m} of leverages,
#'   interpretable as factors by which variances are reduced using
#'   information from the whole network.}
#' \item{w.fixed, w.random}{A vector of length \emph{m} of weights of
#'   individual studies (fixed effects / random effects model).}
#' \item{Q.fixed}{A vector of length \emph{m} of contributions to
#'   total heterogeneity / inconsistency statistic.}
#' \item{TE.fixed, TE.random}{\emph{n}x\emph{n} matrix with estimated
#'   overall treatment effects (fixed effects / random effects model).}
#' \item{seTE.fixed, seTE.random}{\emph{n}x\emph{n} matrix with
#'   standard errors (fixed effects / random effects model).}
#' \item{lower.fixed, upper.fixed, lower.random,
#'   upper.random}{\emph{n}x\emph{n} matrices with lower and upper
#'   confidence interval limits (fixed effects / random effects
#'   model).}
#' \item{zval.fixed, pval.fixed, zval.random,
#'   pval.random}{\emph{n}x\emph{n} matrices with z-value and p-value
#'   for test of overall treatment effect (fixed effects / random
#'   effects model).}
#' \item{seTE.predict}{\emph{n}x\emph{n} matrix with standard errors
#'   for prediction intervals.}
#' \item{lower.predict, upper.predict}{\emph{n}x\emph{n} matrices with
#'   lower and upper prediction interval limits.}
#' \item{prop.direct.fixed, prop.direct.random}{A named vector of the
#'   direct evidence proportion of each network estimate. (fixed
#'   effects / random effects model).}
#' \item{TE.direct.fixed, TE.direct.random}{\emph{n}x\emph{n} matrix
#'   with estimated treatment effects from direct evidence (fixed
#'   effects / random effects model).}
#' \item{seTE.direct.fixed, seTE.direct.random}{\emph{n}x\emph{n}
#'   matrix with estimated standard errors from direct evidence (fixed
#'   effects / random effects model).}
#' \item{lower.direct.fixed, upper.direct.fixed, lower.direct.random,
#'   }{\emph{n}x\emph{n} matrices with lower and upper confidence
#'   interval limits from direct evidence (fixed effects / random
#'   effects model).}
#' \item{ upper.direct.random}{\emph{n}x\emph{n} matrices with lower
#'   and upper confidence interval limits from direct evidence (fixed
#'   effects / random effects model).}
#' \item{zval.direct.fixed, pval.direct.fixed, zval.direct.random,
#'   }{\emph{n}x\emph{n} matrices with z-value and p-value for test of
#'   overall treatment effect from direct evidence (fixed effects /
#'   random effects model).}
#' \item{ pval.direct.random}{\emph{n}x\emph{n} matrices with z-value
#'   and p-value for test of overall treatment effect from direct
#'   evidence (fixed effects / random effects model).}
#' \item{TE.indirect.fixed, TE.indirect.random}{\emph{n}x\emph{n}
#'   matrix with estimated treatment effects from indirect evidence
#'   (fixed effects / random effects model).}
#' \item{seTE.indirect.fixed, seTE.indirect.random}{\emph{n}x\emph{n}
#'   matrix with estimated standard errors from indirect evidence
#'   (fixed effects / random effects model).}
#' \item{lower.indirect.fixed, upper.indirect.fixed,
#'   lower.indirect.random, }{\emph{n}x\emph{n} matrices with lower
#'   and upper confidence interval limits from indirect evidence
#'   (fixed effects / random effects model).}
#' \item{ upper.indirect.random}{\emph{n}x\emph{n} matrices with lower
#'   and upper confidence interval limits from indirect evidence
#'   (fixed effects / random effects model).}
#' \item{zval.indirect.fixed, pval.indirect.fixed,
#'   zval.indirect.random, }{\emph{n}x\emph{n} matrices with z-value
#'   and p-value for test of overall treatment effect from indirect
#'   evidence (fixed effects / random effects model).}
#' \item{pval.indirect.random}{\emph{n}x\emph{n} matrices with z-value
#'   and p-value for test of overall treatment effect from indirect
#'   evidence (fixed effects / random effects model).}
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
#' \item{Q.decomp}{Data frame with columns 'treat1', 'treat2', 'Q',
#'   'df' and 'pval.Q', providing heterogeneity statistics for each
#'   pairwise meta-analysis of direct comparisons.}
#' \item{A.matrix}{Adjacency matrix (\emph{n}x\emph{n}).}
#' \item{X.matrix}{Design matrix (\emph{m}x\emph{n}).}
#' \item{B.matrix}{Edge-vertex incidence matrix (\emph{m}x\emph{n}).}
#' \item{L.matrix}{Laplacian matrix (\emph{n}x\emph{n}).}
#' \item{Lplus.matrix}{Moore-Penrose pseudoinverse of the Laplacian
#'   matrix (\emph{n}x\emph{n}).}
#' \item{Q.matrix}{Matrix of heterogeneity statistics for pairwise
#'   meta-analyses, where direct comparisons exist
#'   (\emph{n}x\emph{n}).}
#' \item{G.matrix}{Matrix with variances and covariances of
#'   comparisons (\emph{m}x\emph{m}). G is defined as
#'   \strong{BL+B^t}.}
#' \item{H.matrix}{Hat matrix (\emph{m}x\emph{m}), defined as
#'   \strong{H = GW = BL+B^tW}.}
#' \item{n.matrix}{\emph{n}x\emph{n} matrix with number of
#'   observations in direct comparisons (if arguments \code{n1} and
#'   \code{n2} are provided).}
#' \item{events.matrix}{\emph{n}x\emph{n} matrix with number of events
#'   in direct comparisons (if arguments \code{event1} and
#'   \code{event2} are provided).}
#' \item{P.fixed, P.random}{\emph{n}x\emph{n} matrix with direct
#'   evidence proportions (fixed effects / random effects model).}
#' \item{Cov.fixed}{Variance-covariance matrix (fixed effects model)}
#' \item{Cov.random}{Variance-covariance matrix (random effects
#'   model)}
#' \item{sm, level, level.comb}{As defined above.}
#' \item{comb.fixed, comb.random}{As defined above.}
#' \item{prediction, level.predict}{As defined above.}
#' \item{reference.group, baseline.reference, all.treatments}{As
#'   defined above.}
#' \item{seq, tau.preset, tol.multiarm, tol.multiarm.se}{As defined
#'   above.}
#' \item{details.chkmultiarm, sep.trts, nchar.trts}{As defined above.}
#' \item{backtransf, title, warn}{As defined above.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{pairwise}}, \code{\link{forest.netmeta}},
#'   \code{\link{netrank}}, \code{\link{metagen}}
#' 
#' @references
#' Jackson D, White IR, Riley RD (2012):
#' Quantifying the impact of between-study heterogeneity in
#' multivariate meta-analyses.
#' \emph{Statistics in Medicine},
#' \bold{31}, 3805--20
#' 
#' Rücker G (2012):
#' Network meta-analysis, electrical networks and graph theory.
#' \emph{Research Synthesis Methods},
#' \bold{3}, 312--24
#' 
#' Rücker G, Schwarzer G (2014):
#' Reduce dimension or reduce weights? Comparing two approaches to
#' multi-arm studies in network meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{33}, 4353--69
#' 
#' Schwarzer G, Carpenter JR, Rücker G (2015):
#' \emph{Meta-Analysis with R (Use-R!)}.
#' Springer International Publishing, Switzerland
#' 
#' Senn S, Gavini F, Magrez D, Scheen A (2013):
#' Issues in performing a network meta-analysis.
#' \emph{Statistical Methods in Medical Research},
#' \bold{22}, 169--89
#' 
#' @examples
#' data(Senn2013)
#' 
#' # Conduct fixed effects network meta-analysis
#' #
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'                 data = Senn2013, sm = "MD",
#'                 comb.random = FALSE)
#' net1
#' net1$Q.decomp
#' 
#' # Comparison with reference group
#' #
#' print(net1, reference = "plac")
#'
#' \dontrun{
#' # Conduct random effects network meta-analysis
#' #
#' net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'                 data = Senn2013, sm = "MD",
#'                 comb.fixed = FALSE)
#' net2
#' 
#' # Change printing order of treatments with placebo last and use
#' # long treatment names
#' #
#' trts <- c("acar", "benf", "metf", "migl", "piog",
#'           "rosi", "sita", "sulf", "vild", "plac")
#' net3 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'                 data = Senn2013, sm = "MD", comb.fixed = FALSE,
#'                 seq = trts, reference = "Placebo")
#' print(summary(net3), digits = 2)
#' }
#' 
#' @export netmeta


netmeta <- function(TE, seTE,
                    treat1, treat2, studlab,
                    data = NULL, subset = NULL,
                    sm,
                    level = gs("level"),
                    level.comb = gs("level.comb"),
                    comb.fixed = gs("comb.fixed"),
                    comb.random = gs("comb.random") | !is.null(tau.preset),
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
                    tol.multiarm = 0.001,
                    tol.multiarm.se = tol.multiarm,
                    details.chkmultiarm = FALSE,
                    ##
                    sep.trts = ":",
                    nchar.trts = 666,
                    ##
                    n1 = NULL,
                    n2 = NULL,
                    event1 = NULL,
                    event2 = NULL,
                    ##
                    backtransf = gs("backtransf"),
                    ##
                    title = "",
                    keepdata = gs("keepdata"),
                    warn = TRUE
                    ) {


  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkchar <- meta:::chkchar
  chklevel <- meta:::chklevel
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  ##
  chklevel(level)
  chklevel(level.comb)
  chklevel(level.predict)
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
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
  chknumeric(tol.multiarm.se, min = 0, single = TRUE)
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
  ## Catch TE, treat1, treat2, seTE, studlab from data:
  ##
  TE <- eval(mf[[match("TE", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  ##
  if (inherits(TE, "pairwise")) {
    is.pairwise <- TRUE
    ##
    sm <- attr(TE, "sm")
    ##
    seTE <- TE$seTE
    treat1 <- TE$treat1
    treat2 <- TE$treat2
    studlab <- TE$studlab
    ##
    if (!is.null(TE$n1))
      n1 <- TE$n1
    if (!is.null(TE$n2))
      n2 <- TE$n2
    if (!is.null(TE$event1))
      event1 <- TE$event1
    if (!is.null(TE$event2))
      event2 <- TE$event2
    ##
    pairdata <- TE
    data <- TE
    ##
    TE <- TE$TE
  }
  else {
    is.pairwise <- FALSE
    if (missing(sm))
      if (!is.null(data) && !is.null(attr(data, "sm")))
        sm <- attr(data, "sm")
      else
        sm <- ""
    ##
    seTE <- eval(mf[[match("seTE", names(mf))]],
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
    ##
    n1 <- eval(mf[[match("n1", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
    ##
    n2 <- eval(mf[[match("n2", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
    ##
    event1 <- eval(mf[[match("event1", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
    ##
    event2 <- eval(mf[[match("event2", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
  }
  ##
  chknumeric(TE)
  chknumeric(seTE)
  ##
  if (!any(!is.na(TE) & !is.na(seTE)))
    stop("Missing data for estimates (argument 'TE') and ",
         "standard errors (argument 'seTE') in all studies.\n  ",
         "No network meta-analysis possible.",
         call. = FALSE)
  ##
  k.Comp <- length(TE)
  ##
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  ##
  if (length(studlab) == 0) {
    if (warn)
      warning("No information given for argument 'studlab'. ",
              "Assuming that comparisons are from independent studies.",
              call. = FALSE)
    studlab <- seq(along = TE)
  }
  studlab <- as.character(studlab)
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  missing.subset <- is.null(subset)
  ##
  if (!is.null(event1) & !is.null(event2))
    available.events <- TRUE
  else
    available.events <- FALSE
  ##
  if (!is.null(n1) & !is.null(n2))
    available.n <- TRUE
  else
    available.n <- FALSE


  ##
  ##
  ## (2b) Store complete dataset in list object data
  ##      (if argument keepdata is TRUE)
  ##
  ##
  if (keepdata) {
    if (nulldata & !is.pairwise)
      data <- data.frame(.studlab = studlab, stringsAsFactors = FALSE)
    else if (nulldata & is.pairwise) {
      data <- pairdata
      data$.studlab <- studlab
    }
    else
      data$.studlab <- studlab
    ##
    data$.order <- seq_along(studlab)
    ##
    data$.treat1 <- treat1
    data$.treat2 <- treat2
    ##
    data$.TE <- TE
    data$.seTE <- seTE
    ##
    data$.event1 <- event1
    data$.n1 <- n1
    data$.event2 <- event2
    data$.n2 <- n2
    ##
    ## Check for correct treatment order within comparison
    ##
    wo <- data$.treat1 > data$.treat2
    ##
    if (any(wo)) {
      data$.TE[wo] <- -data$.TE[wo]
      ttreat1 <- data$.treat1
      data$.treat1[wo] <- data$.treat2[wo]
      data$.treat2[wo] <- ttreat1[wo]
      ##
      if (meta:::isCol(data, ".n1") & meta:::isCol(data, ".n2")) {
        tn1 <- data$.n1
        data$.n1[wo] <- data$.n2[wo]
        data$.n2[wo] <- tn1[wo]
      }
      ##
      if (meta:::isCol(data, ".event1") & meta:::isCol(data, ".event2")) {
        tevent1 <- data$.event1
        data$.event1[wo] <- data$.event2[wo]
        data$.event2[wo] <- tevent1[wo]
      }
    }
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
    TE <- TE[subset]
    seTE <- seTE[subset]
    treat1 <- treat1[subset]
    treat2 <- treat2[subset]
    studlab <- studlab[subset]
    ##
    if (!is.null(n1))
      n1 <- n1[subset]
    if (!is.null(n2))
      n2 <- n2[subset]
    if (!is.null(event1))
      event1 <- event1[subset]
    if (!is.null(event2))
      event2 <- event2[subset]
  }
  ##
  labels <- sort(unique(c(treat1, treat2)))
  ##
  if (compmatch(labels, sep.trts)) {
    if (!missing.sep.trts)
      warning("Separator '", sep.trts,
              "' used in at least one treatment label. ",
              "Try to use predefined separators: ",
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
      stop("All predefined separators (':', '-', '_', '/', '+', '.', '|', '*') are used in at least one treatment label.",
           "\n   Please specify a different character that should be used as separator (argument 'sep.trts').",
           call. = FALSE)
  }
  ##
  if (!is.null(seq))
    seq <- setseq(seq, labels)
  else {
    seq <- labels
    if (is.numeric(seq))
      seq <- as.character(seq)
  }
  ##
  if (reference.group != "")
    reference.group <- setref(reference.group, labels)


  ##
  ##
  ## (4) Additional checks
  ##
  ##
  if (any(treat1 == treat2))
    stop("Treatments must be different (arguments 'treat1' and 'treat2').",
         call. = FALSE)
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
    stop(paste("Study '", names(tabnarms)[sel.narms],
               "' has a wrong number of comparisons.",
               "\n  Please provide data for all treatment comparisons",
               " (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
               sep = ""),
         call. = FALSE)
  if (sum(sel.narms) > 1)
    stop(paste("The following studies have a wrong number of comparisons: ",
               paste(paste("'", names(tabnarms)[sel.narms], "'", sep = ""),
                     collapse = ", "),
               "\n  Please provide data for all treatment comparisons",
               " (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
               sep = ""),
         call. = FALSE)
  ##
  ## Check number of subgraphs
  ##
  n.subnets <- netconnection(treat1, treat2, studlab)$n.subnets
  ##
  if (n.subnets > 1)
    stop(paste("Network consists of ", n.subnets, " separate sub-networks.\n  ",
               "Use R function 'netconnection' to identify sub-networks.",
               sep = ""),
         call. = FALSE)
  ##
  ## Check NAs and zero standard errors
  ##
  excl <- is.na(TE) | is.na(seTE) | seTE <= 0
  ##
  if (any(excl)) {
    if (keepdata)
      data$.excl <- excl
    ##
    dat.NAs <- data.frame(studlab = studlab[excl],
                          treat1 = treat1[excl],
                          treat2 = treat2[excl],
                          TE = format(round(TE[excl], 4)),
                          seTE = format(round(seTE[excl], 4)),
                          stringsAsFactors = FALSE
                          )
    if (warn)
      warning("Comparison",
              if (sum(excl) > 1) "s",
              " with missing TE / seTE or zero seTE not considered ",
              "in network meta-analysis.",
              call. = FALSE)
    if (warn) {
      cat(paste("Comparison",
                if (sum(excl) > 1) "s",
                " not considered in network meta-analysis:\n", sep = ""))
      prmatrix(dat.NAs, quote = FALSE, right = TRUE,
               rowlab = rep("", sum(excl)))
      cat("\n")
    }
    ##
    studlab <- studlab[!(excl)]
    treat1  <- treat1[!(excl)]
    treat2  <- treat2[!(excl)]
    TE      <- TE[!(excl)]
    seTE    <- seTE[!(excl)]
    ##
    if (!is.null(n1))
      n1 <- n1[!excl]
    if (!is.null(n2))
      n2 <- n2[!excl]
    if (!is.null(event1))
      event1 <- event1[!excl]
    if (!is.null(event2))
      event2 <- event2[!excl]
    ##
    seq <- seq[seq %in% unique(c(treat1, treat2))]
    labels <- labels[labels %in% unique(c(treat1, treat2))]
  }
  ##
  ## Check for correct number of comparisons (after removing
  ## comparisons with missing data)
  ##
  tabnarms <- table(studlab)
  sel.narms <- !is.wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  ##
  if (sum(sel.narms) == 1)
    stop(paste("After removing comparisons with missing treatment effects",
               " or standard errors,\n  study '",
               names(tabnarms)[sel.narms],
               "' has a wrong number of comparisons.",
               " Please check data and\n  consider to remove study",
               " from network meta-analysis.",
               sep = ""),
         call. = FALSE)
  if (sum(sel.narms) > 1)
    stop(paste("After removing comparisons with missing treatment effects",
               " or standard errors,\n  the following studies have",
               " a wrong number of comparisons: ",
               paste(paste("'", names(tabnarms)[sel.narms], "'", sep = ""),
                     collapse = ", "),
               "\n  Please check data and consider to remove studies",
               " from network meta-analysis.",
               sep = ""),
         call. = FALSE)
  ##
  ## Check number of subgraphs
  ##
  n.subnets <- netconnection(treat1, treat2, studlab)$n.subnets
  ##
  if (n.subnets > 1)
    stop(paste("After removing comparisons with missing treatment effects",
               " or standard errors,\n  network consists of ",
               n.subnets, " separate sub-networks.\n  ",
               "Please check data and consider to remove studies",
               " from network meta-analysis.",
               sep = ""),
         call. = FALSE)
  ##
  ## Check for correct treatment order within comparison
  ##
  wo <- treat1 > treat2
  ##
  if (any(wo)) {
    TE[wo] <- -TE[wo]
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
    ##
    if (available.n) {
      tn1 <- n1
      n1[wo] <- n2[wo]
      n2[wo] <- tn1[wo]
    }
    ##
    if (available.events) {
      tevent1 <- event1
      event1[wo] <- event2[wo]
      event2[wo] <- tevent1[wo]
    }
  }


  ##
  ##
  ## (5) Generate analysis dataset
  ##
  ##
  ##
  ## Generate ordered data set, with added numbers of arms per study
  ##
  p0 <- prepare(TE, seTE, treat1, treat2, studlab)
  ##
  ## Check consistency of treatment effects and standard errors in
  ## multi-arm studies
  ##
  chkmultiarm(p0$TE, p0$seTE, p0$treat1, p0$treat2, p0$studlab,
              tol.multiarm = tol.multiarm, tol.multiarm.se = tol.multiarm.se,
              details = details.chkmultiarm)
  ##
  ## Study overview
  ##
  tdata <- data.frame(studies = p0$studlab, narms = p0$narms,
                      order = p0$order,
                      stringsAsFactors = FALSE)
  ##
  tdata <- tdata[!duplicated(tdata[, c("studies", "narms")]), , drop = FALSE]
  studies <- tdata$studies[order(tdata$order)]
  narms <- tdata$narms[order(tdata$order)]
  
  
  ##
  ##
  ## (6) Conduct network meta-analysis
  ##
  ##
  ## Fixed effects model
  ##
  res.f <- nma.ruecker(p0$TE, sqrt(1 / p0$weights),
                       p0$treat1, p0$treat2,
                       p0$treat1.pos, p0$treat2.pos,
                       p0$narms, p0$studlab,
                       sm,
                       level, level.comb,
                       p0$seTE, sep.trts = sep.trts)
  ##
  ## Random effects model
  ##
  if (is.null(tau.preset))
    tau <- res.f$tau
  else
    tau <- tau.preset
  ##
  p1 <- prepare(TE, seTE, treat1, treat2, studlab, tau)
  ##
  res.r <- nma.ruecker(p1$TE, sqrt(1 / p1$weights),
                       p1$treat1, p1$treat2,
                       p1$treat1.pos, p1$treat2.pos,
                       p1$narms, p1$studlab,
                       sm,
                       level, level.comb,
                       p1$seTE, tau, sep.trts = sep.trts)
  ##
  TE.random <- res.r$TE.pooled
  seTE.random <- res.r$seTE.pooled
  df.Q <- res.f$df
  ##
  ## Prediction intervals
  ##
  if (df.Q == 0)
    prediction <- FALSE
  ##
  if (df.Q >= 2) {
    seTE.predict <- sqrt(seTE.random^2 + tau^2)
    ci.p <- ci(TE.random, seTE.predict, level.predict, df.Q - 1)
    p.lower <- ci.p$lower
    p.upper <- ci.p$upper
    diag(p.lower) <- 0
    diag(p.upper) <- 0
  }
  else {
    seTE.predict <- p.lower <- p.upper <- seTE.random
    seTE.predict[!is.na(seTE.predict)] <- NA
    p.lower[!is.na(p.lower)] <- NA
    p.upper[!is.na(p.upper)] <- NA
  }


  ##
  ##
  ## (7) Generate R object
  ##
  ##
  trts <- rownames(res.f$A.matrix)
  ##
  o <- order(p0$order)
  ##
  res <- list(studlab = res.f$studlab[o],
              treat1 = res.f$treat1[o],
              treat2 = res.f$treat2[o],
              ##
              TE = res.f$TE[o],
              seTE = res.f$seTE.orig[o],
              seTE.adj = res.f$seTE[o],
              ##
              event1 = event1,
              event2 = event2,
              n1 = n1,
              n2 = n2,
              ##
              k = res.f$k,
              m = res.f$m,
              n = res.f$n,
              d = NA,
              ##
              trts = trts,
              k.trts = rowSums(res.f$A.matrix),
              n.trts = if (available.n) NA else NULL,
              events.trts = if (available.events) NA else NULL,
              ##
              n.arms = NA,
              multiarm = NA,
              ##
              studies = studies,
              narms = narms,
              ##
              designs = NA,
              ##
              TE.nma.fixed = res.f$TE.nma[o],
              seTE.nma.fixed = res.f$seTE.nma[o],
              lower.nma.fixed = res.f$lower.nma[o],
              upper.nma.fixed = res.f$upper.nma[o],
              zval.nma.fixed = res.f$zval.nma[o],
              pval.nma.fixed = res.f$pval.nma[o],
              ##
              leverage.fixed = res.f$leverage[o],
              w.fixed = res.f$w.pooled[o],
              Q.fixed = res.f$Q.pooled[o],
              ##
              TE.fixed = res.f$TE.pooled,
              seTE.fixed = res.f$seTE.pooled,
              lower.fixed = res.f$lower.pooled,
              upper.fixed = res.f$upper.pooled,
              zval.fixed = res.f$zval.pooled,
              pval.fixed = res.f$pval.pooled,
              ##
              TE.nma.random = res.r$TE.nma[o],
              seTE.nma.random = res.r$seTE.nma[o],
              lower.nma.random = res.r$lower.nma[o],
              upper.nma.random = res.r$upper.nma[o],
              zval.nma.random = res.r$zval.nma[o],
              pval.nma.random = res.r$pval.nma[o],
              ##
              w.random = res.r$w.pooled[o],
              ##
              TE.random = TE.random,
              seTE.random = seTE.random,
              lower.random = res.r$lower.pooled,
              upper.random = res.r$upper.pooled,
              zval.random = res.r$zval.pooled,
              pval.random = res.r$pval.pooled,
              ##
              seTE.predict = seTE.predict,
              lower.predict = p.lower,
              upper.predict = p.upper,
              ##
              prop.direct.fixed = NA,
              prop.direct.random = NA,
              ##
              TE.direct.fixed = res.f$TE.direct,
              seTE.direct.fixed = res.f$seTE.direct,
              lower.direct.fixed = res.f$lower.direct,
              upper.direct.fixed = res.f$upper.direct,
              zval.direct.fixed = res.f$zval.direct,
              pval.direct.fixed = res.f$pval.direct,
              ##
              TE.direct.random = res.r$TE.direct,
              seTE.direct.random = res.r$seTE.direct,
              lower.direct.random = res.r$lower.direct,
              upper.direct.random = res.r$upper.direct,
              zval.direct.random = res.r$zval.direct,
              pval.direct.random = res.r$pval.direct,
              ##
              TE.indirect.fixed = NA,
              seTE.indirect.fixed = NA,
              lower.indirect.fixed = NA,
              upper.indirect.fixed = NA,
              zval.indirect.fixed = NA,
              pval.indirect.fixed = NA,
              ##
              TE.indirect.random = NA,
              seTE.indirect.random = NA,
              lower.indirect.random = NA,
              upper.indirect.random = NA,
              zval.indirect.random = NA,
              pval.indirect.random = NA,
              ##
              Q = res.f$Q,
              df.Q = df.Q,
              pval.Q = res.f$pval.Q,
              I2 = res.f$I2,
              lower.I2 = res.f$lower.I2,
              upper.I2 = res.f$upper.I2,
              tau = tau,
              ##
              Q.heterogeneity = NA,
              df.Q.heterogeneity = NA,
              pval.Q.heterogeneity = NA,
              ##
              Q.inconsistency = NA,
              df.Q.inconsistency = NA,
              pval.Q.inconsistency = NA,
              ##
              Q.decomp = res.f$Q.decomp,
              ##
              A.matrix = res.f$A.matrix,
              X.matrix = res.f$B.matrix[o, ],
              B.matrix = res.f$B.matrix[o, ],
              L.matrix = res.f$L.matrix,
              Lplus.matrix = res.f$Lplus.matrix,
              Q.matrix = res.f$Q.matrix,
              ##
              G.matrix = res.f$G.matrix[o, o],
              H.matrix = res.f$H.matrix[o, o],
              ##
              n.matrix = if (available.n) NA else NULL,
              events.matrix = if (available.events) NA else NULL,
              ##
              P.fixed = NA,
              P.random = NA,
              ##
              Cov.fixed = res.f$Cov,
              Cov.random = res.r$Cov,
              ##
              treat1.pos = res.f$treat1.pos[o],
              treat2.pos = res.f$treat2.pos[o],
              ##
              sm = sm,
              method = "Inverse",
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
              tol.multiarm.se = tol.multiarm.se,
              details.chkmultiarm = details.chkmultiarm,
              ##
              sep.trts = sep.trts,
              nchar.trts = nchar.trts,
              ##
              backtransf = backtransf,
              ##
              title = title,
              ##
              warn = warn,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- "netmeta"
  ##
  ## Add results for indirect treatment estimates
  ##
  n <- res$n
  ##
  res$prop.direct.fixed  <- netmeasures(res, random = FALSE,
                                        warn = warn)$proportion
  ## Print warning(s) in call of netmeasures() once
  res$prop.direct.random <-
    suppressWarnings(netmeasures(res, random = TRUE,
                                 tau.preset = res$tau,
                                 warn = FALSE)$proportion)
  if (is.logical(res$prop.direct.fixed))
    res$prop.direct.fixed <- as.numeric(res$prop.direct.fixed)
  if (is.logical(res$prop.direct.random))
    res$prop.direct.random <- as.numeric(res$prop.direct.random)
  ##
  P.fixed <- P.random <- matrix(NA, n, n)
  colnames(P.fixed) <- rownames(P.fixed) <-
    colnames(P.random) <- rownames(P.random) <- trts
  ##
  if (n == 2) {
    ##
    ## For two treatments only direct evidence is available
    ##
    res$prop.direct.fixed <- 1
    res$prop.direct.random <- 1
    names(res$prop.direct.fixed) <-
      names(res$prop.direct.random) <- paste(labels, collapse = sep.trts)
    ##
    sel <- row(P.fixed) != col(P.fixed)
    P.fixed[sel] <- 1
    P.random[sel] <- 1
  }
  else {
    k <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        k <- k + 1
        P.fixed[i, j] <- P.fixed[j, i] <- res$prop.direct.fixed[k]
        P.random[i, j] <- P.random[j, i] <- res$prop.direct.random[k]
      }
    }
  }
  ##
  ## Set direct evidence estimates to 0 if only indirect evidence is available
  ## (otherwise indirect estimates would be NA as direct estimates are NA)
  ##
  TE.direct.fixed <- res$TE.direct.fixed
  TE.direct.random <- res$TE.direct.random
  ##
  TE.direct.fixed[abs(P.fixed) < .Machine$double.eps^0.5] <- 0
  TE.direct.random[abs(P.random) < .Machine$double.eps^0.5] <- 0
  ##
  ## Indirect estimate is NA if only direct evidence is available
  ##
  res$P.fixed <- P.fixed
  res$P.random <- P.random
  ##
  P.fixed[abs(P.fixed - 1) < .Machine$double.eps^0.5] <- NA
  P.random[abs(P.random - 1) < .Machine$double.eps^0.5] <- NA
  ##
  ## Fixed effects model
  ##
  ci.if <- ci((res$TE.fixed - P.fixed * TE.direct.fixed) / (1 - P.fixed),
              sqrt(res$seTE.fixed^2 / (1 - P.fixed)),
              level = level)
  ##
  res$TE.indirect.fixed   <- ci.if$TE
  res$seTE.indirect.fixed <- ci.if$seTE
  ##
  res$lower.indirect.fixed <- ci.if$lower
  res$upper.indirect.fixed <- ci.if$upper
  ##
  res$zval.indirect.fixed <- ci.if$z
  res$pval.indirect.fixed <- ci.if$p
  ##
  ## Random effects model
  ##
  ci.ir <- ci((res$TE.random - P.random * TE.direct.random) / (1 - P.random),
              sqrt(res$seTE.random^2 / (1 - P.random)),
              level = level)
  ##
  res$TE.indirect.random   <- ci.ir$TE
  res$seTE.indirect.random <- ci.ir$seTE
  ##
  res$lower.indirect.random <- ci.ir$lower
  res$upper.indirect.random <- ci.ir$upper
  ##
  res$zval.indirect.random <- ci.ir$z
  res$pval.indirect.random <- ci.ir$p
  ##
  ## Number of designs
  ##
  krahn <- nma.krahn(res)
  res$d <- krahn$d
  if (is.null(res$d))
    res$d <- 1
  ##
  if (is.null(krahn$design$design))
    res$designs <- rownames(res$Cov.fixed)
  else
    res$designs <- as.character(krahn$design$design)
  ##
  res$designs <- unique(res$designs)
  ##
  ##
  ##
  if (any(res$narms > 2)) {
    tdata1 <- data.frame(studlab = res$studlab,
                         .order = seq(along = res$studlab))
    tdata2 <- data.frame(studlab = as.character(res$studies),
                         narms = res$narms)
    ##
    tdata12 <- merge(tdata1, tdata2,
                     by = "studlab", all.x = TRUE, all.y = FALSE,
                     sort = FALSE)
    tdata12 <- tdata12[order(tdata12$.order), ]
    res$n.arms <- tdata12$narms
    res$multiarm <- tdata12$narms > 2
  }
  else {
    res$n.arms <- rep(2, length(res$studlab))
    res$multiarm <- rep(FALSE, length(res$studlab))
  }
  
  
  ##
  ## Calculate heterogeneity and inconsistency statistics
  ##
  if (res$d > 1) {
    dd <- decomp.design(res, warn = FALSE)
    res$Q.heterogeneity <- dd$Q.decomp$Q[2]
    res$Q.inconsistency <- dd$Q.decomp$Q[3]
    ##
    res$df.Q.heterogeneity <- dd$Q.decomp$df[2]
    res$df.Q.inconsistency <- dd$Q.decomp$df[3]
    ##
    res$pval.Q.heterogeneity <- dd$Q.decomp$pval[2]
    res$pval.Q.inconsistency <- dd$Q.decomp$pval[3]
  }


  if (keepdata) {
    if (is.null(krahn))
      ddat <- data.frame(.studlab = data$.studlab,
                         .design = paste(data$.treat1, data$.treat2,
                                         sep = sep.trts),
                         stringsAsFactors = FALSE)
    else {
      ddat <- unique(krahn$studies[, c("studlab", "design")])
      names(ddat) <- paste0(".", names(ddat))
    }

    data <- merge(data,
                  data.frame(.studlab = res$studies,
                             .narms = res$narms),
                  by = ".studlab",
                  stringsAsFactors = FALSE)
    ##
    res$data <- merge(data, ddat,
                      by = ".studlab",
                      suffixes = c(".orig", ""),
                      stringsAsFactors = FALSE)
    res$data$.design <- as.character(res$data$.design)
    res$data <- res$data[order(res$data$.order), ]
    res$data$.order <- NULL
  }
  ##
  if (available.n) {
    res$n.matrix <- netmatrix(res, n1 + n2, func = "sum")
    ##
    dat.n <- data.frame(studlab = c(studlab, studlab),
                        treat = c(treat1, treat2),
                        n = c(n1, n2))
    dat.n <- dat.n[!duplicated(dat.n[, c("studlab", "treat")]), ]
    dat.n <- by(dat.n$n, dat.n$treat, sum, na.rm = TRUE)
    res$n.trts <- as.vector(dat.n[trts])
    names(res$n.trts) <- trts
  }
  ##
  if (available.events) {
    res$events.matrix <- netmatrix(res, event1 + event2, func = "sum")
    ##
    dat.e <- data.frame(studlab = c(studlab, studlab),
                        treat = c(treat1, treat2),
                        n = c(event1, event2))
    dat.e <- dat.e[!duplicated(dat.e[, c("studlab", "treat")]), ]
    dat.e <- by(dat.e$n, dat.e$treat, sum, na.rm = TRUE)
    res$events.trts <- as.vector(dat.e[trts])
    names(res$events.trts) <- trts
  }


  res
}
