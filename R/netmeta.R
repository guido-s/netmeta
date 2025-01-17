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
#'   or log hazard ratio). Or an R object created with
#'   \code{\link[meta]{pairwise}}.
#' @param seTE Standard error of treatment estimate.
#' @param treat1 Label/Number for first treatment.
#' @param treat2 Label/Number for second treatment.
#' @param studlab An optional - but important! - vector with study
#'   labels (see Details).
#' @param data An optional data frame containing the study
#'   information.
#' @param subset An optional vector specifying a subset of studies to
#'   be used.
#' @param correlated An optional logical vector specifying whether
#'   treatment arms of a multi-arm study are correlated.
#' @param sm A character string indicating underlying summary measure,
#'   e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
#'   \code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param level The level used to calculate confidence intervals for
#'   individual comparisons.
#' @param level.ma The level used to calculate confidence intervals
#'   for network estimates.
#' @param common A logical indicating whether a common effects network
#'   meta-analysis should be conducted.
#' @param random A logical indicating whether a random effects network
#'   meta-analysis should be conducted.
#' @param prediction A logical indicating whether prediction intervals
#'   should be printed.
#' @param level.predict The level used to calculate prediction
#'   intervals for a new study.
#' @param reference.group Reference treatment (first treatment is used
#'   if argument is missing).
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
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau}. Either \code{"DL"}, \code{"REML"}, or
#'   \code{"ML"}, can be abbreviated.
#' @param tau.preset An optional value for manually setting the
#'   square-root of the between-study variance \eqn{\tau^2}.
#' @param tol.multiarm A numeric for the tolerance for consistency of
#'   treatment estimates in multi-arm studies which are consistent by
#'   design.
#' @param tol.multiarm.se A numeric for the tolerance for consistency
#'   of standard errors in multi-arm studies which are consistent by
#'   design. This check is not conducted if the argument is
#'   \code{NULL}.
#' @param details.chkmultiarm A logical indicating whether treatment
#'   estimates and / or variances of multi-arm studies with
#'   inconsistent results or negative multi-arm variances should be
#'   printed.
#' @param sep.trts A character used in comparison names as separator
#'   between treatment labels.
#' @param overall.hetstat A logical indicating whether to print heterogeneity
#'   measures.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
#' @param nchar.studlab A numeric defining the minimum number of
#'   characters used to create unique study labels.
#' @param func.inverse R function used to calculate the pseudoinverse
#'   of the Laplacian matrix L (see Details).
#' @param n1 Number of observations in first treatment group.
#' @param n2 Number of observations in second treatment group.
#' @param event1 Number of events in first treatment group.
#' @param event2 Number of events in second treatment group.
#' @param incr Numerical value added to cell frequencies (for details,
#'   see \code{\link[meta]{pairwise}}).
#' @param mean1 Mean in first treatment group.
#' @param mean2 Mean in second treatment group.
#' @param sd1 Standard deviation in first treatment group.
#' @param sd2 Standard deviation in second treatment group.
#' @param time1 Person time at risk in first treatment group.
#' @param time2 Person time at risk in second treatment group.
#' @param title Title of meta-analysis / systematic review.
#' @param keepdata A logical indicating whether original data(set)
#'   should be kept in netmeta object.
#' @param keeprma A logical indicating whether \code{\link[metafor]{rma.mv}}
#'   object should be stored.
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance \eqn{\tau^2}. This argument
#'   is passed on to \code{\link[metafor]{rma.mv}}.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if studies are excluded from meta-analysis due to zero
#'   standard errors).
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param nchar Deprecated argument (replaced by \code{nchar.trts}).
#' @param \dots Additional arguments (to catch deprecated arguments).
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
#' \bold{meta}, e.g. \code{\link[meta]{metabin}} for binary outcomes or
#' \code{\link[meta]{metacont}} for continuous outcomes, can be used to
#' calculate treatment effects separately for each treatment
#' comparison which is a rather tedious enterprise. If data are
#' provided in \emph{arm-based} format, that is, data are given for
#' each treatment arm separately (e.g. number of events and
#' participants for binary outcomes), a much more convenient way to
#' transform data into contrast-based form is available. Function
#' \code{\link[meta]{pairwise}} can automatically transform data with binary
#' outcomes (using the \code{\link[meta]{metabin}} function from R package
#' \bold{meta}), continuous outcomes (\code{\link[meta]{metacont}}
#' function), incidence rates (\code{\link[meta]{metainc}} function), and
#' generic outcomes (\code{\link[meta]{metagen}} function). Additional
#' arguments of these functions can be provided (see help page of
#' function \code{\link[meta]{pairwise}}).
#' 
#' Note, all pairwise comparisons must be provided for a multi-arm
#' study. Consider a multi-arm study of \emph{p} treatments with known
#' variances. For this study, treatment effects and standard errors
#' must be provided for each of \emph{p}(\emph{p} - 1) / 2 possible
#' comparisons. For instance, a three-arm study contributes three
#' pairwise comparisons, a four-arm study even six pairwise
#' comparisons. Function \code{\link[meta]{pairwise}} automatically
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
#' Internally, both common and random effects models are calculated
#' regardless of values choosen for arguments \code{common} and
#' \code{random}. Accordingly, the network estimates for the random
#' effects model can be extracted from component \code{TE.random} of
#' an object of class \code{"netmeta"} even if argument \code{random =
#' FALSE}. However, all functions in R package \bold{netmeta} will
#' adequately consider the values for \code{common} and
#' \code{random}. E.g. function \code{\link{print.summary.netmeta}}
#' will not print results for the random effects model if \code{random
#' = FALSE}.
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
#' @note
#' R function \code{\link[metafor]{rma.mv}} from R package
#' \pkg{metafor} (Viechtbauer 2010) is called internally to estimate
#' the between-study variance \eqn{\tau^2} for the (restricted)
#' maximum likelihood method. For binary outcomes, incidence rates,
#' and the mean difference, the variance-covariance matrix is
#' calculated if arguments \code{event1}, \code{event2}, \code{n1},
#' and \code{n2} (binary outcomes); \code{event1}, \code{event2},
#' \code{time1}, and \code{time2} (incidence rates); \code{n1},
#' \code{n2}, \code{sd1}, and \code{sd2} (mean difference) are
#' provided. For datasets preprocessed with \code{\link[meta]{pairwise}}
#' the respective variables are selected automatically.
#'
#' @return
#' An object of class \code{netmeta} with corresponding \code{print},
#' \code{summary}, \code{forest}, and \code{netrank} functions. The
#' object is a list containing the following components:
#' \item{studlab, treat1, treat2, TE, seTE}{As defined above.}
#' \item{seTE.adj.common, seTE.adj.random}{Standard error of treatment
#'   estimate, adjusted for multi-arm studies.}
#' \item{design}{Design of study providing pairwise comparison.}
#' \item{n1, n2, event1, event2, incr}{As defined above.}
#' \item{mean1, mean2, sd1, sd2, time1, time2}{As defined above.}
#' \item{sd1, sd2, time1, time2}{As defined above.}
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
#' \item{comparisons}{Vector with unique direct comparisons present in the
#'   network.}
#' \item{TE.nma.common, TE.nma.random}{A vector of length \emph{m} of
#'   consistent treatment effects estimated by network meta-analysis
#'   (nma) (common / random effects model).}
#' \item{seTE.nma.common, seTE.nma.random}{A vector of length \emph{m}
#'   of effective standard errors estimated by network meta-analysis
#'   (common / random effects model).}
#' \item{lower.nma.common, lower.nma.random}{A vector of length
#'   \emph{m} of lower confidence interval limits for consistent
#'   treatment effects estimated by network meta-analysis (common
#'   effects / random effects model).}
#' \item{upper.nma.common, upper.nma.random}{A vector of length
#'   \emph{m} of upper confidence interval limits for the consistent
#'   treatment effects estimated by network meta-analysis (common
#'   effects / random effects model).}
#' \item{statistic.nma.common, statistic.nma.random}{A vector of length
#'   \emph{m} of z-values for test of treatment effect for individual
#'   comparisons (common / random effects model).}
#' \item{pval.nma.common, pval.nma.random}{A vector of length \emph{m}
#'   of p-values for test of treatment effect for individual
#'   comparisons (common / random effects model).}
#' \item{leverage.common}{A vector of length \emph{m} of leverages,
#'   interpretable as factors by which variances are reduced using
#'   information from the whole network.}
#' \item{w.common, w.random}{A vector of length \emph{m} of weights of
#'   individual studies (common / random effects model).}
#' \item{Q.common}{A vector of length \emph{m} of contributions to
#'   total heterogeneity / inconsistency statistic.}
#' \item{TE.common, TE.random}{\emph{n}x\emph{n} matrix with estimated
#'   overall treatment effects (common / random effects model).}
#' \item{seTE.common, seTE.random}{\emph{n}x\emph{n} matrix with
#'   standard errors (common / random effects model).}
#' \item{lower.common, upper.common, lower.random,
#'   upper.random}{\emph{n}x\emph{n} matrices with lower and upper
#'   confidence interval limits (common / random effects
#'   model).}
#' \item{statistic.common, pval.common, statistic.random,
#'   pval.random}{\emph{n}x\emph{n} matrices with z-value and p-value
#'   for test of overall treatment effect (common / random
#'   effects model).}
#' \item{seTE.predict}{\emph{n}x\emph{n} matrix with standard errors
#'   for prediction intervals.}
#' \item{lower.predict, upper.predict}{\emph{n}x\emph{n} matrices with
#'   lower and upper prediction interval limits.}
#' \item{prop.direct.common, prop.direct.random}{A named vector of the
#'   direct evidence proportion of each network estimate. (common
#'   effects / random effects model).}
#' \item{TE.direct.common, TE.direct.random}{\emph{n}x\emph{n} matrix
#'   with estimated treatment effects from direct evidence (common
#'   effects / random effects model).}
#' \item{seTE.direct.common, seTE.direct.random}{\emph{n}x\emph{n}
#'   matrix with estimated standard errors from direct evidence (common
#'   effects / random effects model).}
#' \item{lower.direct.common, upper.direct.common, lower.direct.random,
#'   }{\emph{n}x\emph{n} matrices with lower and upper confidence
#'   interval limits from direct evidence (common / random
#'   effects model).}
#' \item{upper.direct.random}{\emph{n}x\emph{n} matrices with lower
#'   and upper confidence interval limits from direct evidence (common
#'   effects / random effects model).}
#' \item{statistic.direct.common, pval.direct.common,
#'   statistic.direct.random, }{\emph{n}x\emph{n} matrices with
#'   z-value and p-value for test of overall treatment effect from
#'   direct evidence (common / random effects model).}
#' \item{pval.direct.random}{\emph{n}x\emph{n} matrices with z-value
#'   and p-value for test of overall treatment effect from direct
#'   evidence (common / random effects model).}
#' \item{TE.indirect.common, TE.indirect.random}{\emph{n}x\emph{n}
#'   matrix with estimated treatment effects from indirect evidence
#'   (common / random effects model).}
#' \item{seTE.indirect.common, seTE.indirect.random}{\emph{n}x\emph{n}
#'   matrix with estimated standard errors from indirect evidence
#'   (common / random effects model).}
#' \item{lower.indirect.common, upper.indirect.common,
#'   lower.indirect.random, }{\emph{n}x\emph{n} matrices with lower
#'   and upper confidence interval limits from indirect evidence
#'   (common / random effects model).}
#' \item{upper.indirect.random}{\emph{n}x\emph{n} matrices with lower
#'   and upper confidence interval limits from indirect evidence
#'   (common / random effects model).}
#' \item{statistic.indirect.common, pval.indirect.common,
#'   statistic.indirect.random, }{\emph{n}x\emph{n} matrices with
#'   z-value and p-value for test of overall treatment effect from
#'   indirect evidence (common / random effects model).}
#' \item{pval.indirect.random}{\emph{n}x\emph{n} matrices with z-value
#'   and p-value for test of overall treatment effect from indirect
#'   evidence (common / random effects model).}
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
#' \item{L.matrix.common, L.matrix.random}{Laplacian matrix
#'   (\emph{n}x\emph{n}).}
#' \item{Lplus.matrix.common, Lplus.matrix.random}{Moore-Penrose
#'   pseudoinverse of the Laplacian matrix (\emph{n}x\emph{n}).}
#' \item{Q.matrix}{Matrix of heterogeneity statistics for pairwise
#'   meta-analyses, where direct comparisons exist
#'   (\emph{n}x\emph{n}).}
#' \item{G.matrix}{Matrix with variances and covariances of
#'   comparisons (\emph{m}x\emph{m}). G is defined as
#'   \strong{BL+B^t}.}
#' \item{H.matrix.common, H.matrix.random}{Hat matrix
#'   (\emph{m}x\emph{m}), defined as \strong{H = GW = BL+B^tW}.}
#' \item{n.matrix}{\emph{n}x\emph{n} matrix with number of
#'   observations in direct comparisons (if arguments \code{n1} and
#'   \code{n2} are provided).}
#' \item{events.matrix}{\emph{n}x\emph{n} matrix with number of events
#'   in direct comparisons (if arguments \code{event1} and
#'   \code{event2} are provided).}
#' \item{P.common, P.random}{\emph{n}x\emph{n} matrix with direct
#'   evidence proportions (common / random effects model).}
#' \item{Cov.common}{Variance-covariance matrix (common effects model)}
#' \item{Cov.random}{Variance-covariance matrix (random effects model)}
#' \item{sm, level, level.ma}{As defined above.}
#' \item{common, random}{As defined above.}
#' \item{prediction, level.predict}{As defined above.}
#' \item{reference.group, baseline.reference, small.values,
#'   all.treatments}{As defined above.}
#' \item{seq, tau.preset, tol.multiarm, tol.multiarm.se}{As defined
#'   above.}
#' \item{details.chkmultiarm, sep.trts, nchar.trts}{As defined above.}
#' \item{backtransf, title, warn, warn.deprecated}{As defined above.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create object.}
#'
#' In addition, the following component is stored if \bold{metafor} is
#' used to calculate the between-study variance and argument
#' \code{keeprma = TRUE}:
#' \item{rma.tau}{R object created with \code{\link[metafor]{rma.mv}}.}
#' 
#' @author Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}, Guido
#'   Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link[meta]{pairwise}}, \code{\link{forest.netmeta}},
#'   \code{\link{netrank}}, \code{\link[meta]{metagen}}
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
#' \emph{Meta-Analysis with R (Use R!)}.
#' Springer International Publishing, Switzerland
#' 
#' Senn S, Gavini F, Magrez D, Scheen A (2013):
#' Issues in performing a network meta-analysis.
#' \emph{Statistical Methods in Medical Research},
#' \bold{22}, 169--89
#' 
#' Viechtbauer W (2010):
#' Conducting Meta-Analyses in R with the metafor Package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#' 
#' @examples
#' data(smokingcessation)
#' 
#' # Transform data from arm-based format to contrast-based format
#' #
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' 
#' # Conduct random effects network meta-analysis
#' #
#' net1 <- netmeta(p1, common = FALSE)
#' net1
#' 
#' \dontrun{
#' data(Senn2013)
#' 
#' # Conduct common effects network meta-analysis
#' #
#' net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD", random = FALSE)
#' net2
#' net2$Q.decomp
#' 
#' # Comparison with reference group
#' #
#' print(net2, reference = "plac")
#'
#' # Conduct random effects network meta-analysis
#' #
#' net3 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD", common = FALSE)
#' net3
#' 
#' # Change printing order of treatments with placebo last and use
#' # long treatment names
#' #
#' trts <- c("acar", "benf", "metf", "migl", "piog",
#'   "rosi", "sita", "sulf", "vild", "plac")
#' net4 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'   data = Senn2013, sm = "MD", common = FALSE,
#'   seq = trts, reference = "Placebo")
#' print(net4, digits = 2)
#' }
#' 
#' @export netmeta

netmeta <- function(TE, seTE,
                    treat1, treat2, studlab,
                    data = NULL, subset = NULL,
                    correlated,
                    #
                    sm,
                    level = gs("level"),
                    level.ma = gs("level.ma"),
                    common = gs("common"),
                    random = gs("random") | !is.null(tau.preset),
                    #
                    prediction = gs("prediction"),
                    level.predict = gs("level.predict"),
                    #
                    reference.group,
                    baseline.reference = gs("baseline.reference"),
                    small.values = gs("small.values"),
                    all.treatments = gs("all.treatments"),
                    seq = gs("seq"),
                    #
                    method.tau = gs("method.tau.netmeta"),
                    tau.preset = NULL,
                    #
                    tol.multiarm = gs("tol.multiarm"),
                    tol.multiarm.se = gs("tol.multiarm.se"),
                    details.chkmultiarm = gs("details.chkmultiarm"),
                    #
                    sep.trts = gs("sep.trts"),
                    nchar.trts = gs("nchar.trts"),
                    nchar.studlab = gs("nchar.studlab"),
                    #
                    func.inverse = invmat,
                    #
                    n1 = NULL,
                    n2 = NULL,
                    event1 = NULL,
                    event2 = NULL,
                    incr = NULL,
                    mean1 = NULL,
                    mean2 = NULL,
                    sd1 = NULL,
                    sd2 = NULL,
                    time1 = NULL,
                    time2 = NULL,
                    #
                    overall.hetstat = gs("overall.hetstat"),
                    backtransf = gs("backtransf"),
                    #
                    title = gs("title"),
                    keepdata = gs("keepdata"),
                    keeprma = gs("keeprma"),
                    control = NULL,
                    #
                    warn = gs("warn"), warn.deprecated = gs("warn.deprecated"),
                    #
                    nchar = nchar.trts,
                    ...) {


  #
  #
  # (1) Check arguments
  #
  #
  chklevel(level)
  chklevel(level.predict)
  #
  chklogical(prediction)
  #
  missing.reference.group <- missing(reference.group)
  #
  baseline.reference <- replaceNULL(baseline.reference, TRUE)
  chklogical(baseline.reference)
  #
  small.values <- setsv(replaceNULL(small.values, "desirable"))
  #
  if (!is.null(all.treatments))
    chklogical(all.treatments)
  #
  method.tau <- setchar(method.tau, c("DL", "ML", "REML"))
  #
  if (!is.null(tau.preset))
    chknumeric(tau.preset, min = 0, length = 1)
  #
  tol.multiarm <- replaceNULL(tol.multiarm, 0.001)
  chknumeric(tol.multiarm, min = 0, length = 1)
  if (!is.null(tol.multiarm.se))
    chknumeric(tol.multiarm.se, min = 0, length = 1)
  #
  details.chkmultiarm <- replaceNULL(details.chkmultiarm, FALSE)
  chklogical(details.chkmultiarm)
  #
  missing.sep.trts <- missing(sep.trts)
  sep.trts <- replaceNULL(sep.trts, ":")
  chkchar(sep.trts, length = 1)
  #
  nchar.studlab <- replaceNULL(nchar.studlab, 666)
  chknumeric(nchar.studlab, length = 1)
  #
  overall.hetstat <- replaceNULL(overall.hetstat, TRUE)
  chklogical(overall.hetstat)
  chklogical(backtransf)
  #
  chkchar(title)
  chklogical(keepdata)
  chklogical(keeprma)
  chklogical(warn)
  #
  chklogical(baseline.reference)
  #
  # Check for deprecated arguments in '...'
  #
  args <- list(...)
  chklogical(warn.deprecated)
  #
  level.ma <- deprecated(level.ma, missing(level.ma), args, "level.comb",
                         warn.deprecated)
  chklevel(level.ma)
  #
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  #
  random <-
    deprecated(random, missing(random), args, "comb.random", warn.deprecated)
  chklogical(random)
  #
  missing.nchar.trts <- missing(nchar.trts)
  nchar.trts <- replaceNULL(nchar.trts, 666)
  nchar.trts <-
    deprecated2(nchar.trts, missing.nchar.trts, nchar, missing(nchar),
                warn.deprecated)
  chknumeric(nchar.trts, min = 1, length = 1)
  
  
  #
  #
  # (2) Read data
  #
  #
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  if (nulldata)
    data <- sfsp
  #
  # Catch TE, treat1, treat2, seTE, studlab from data:
  #
  TE <- catch("TE", mc, data, sfsp)
  #
  avail.reference.group.pairwise <- FALSE
  #
  if (is.data.frame(TE) &&
      (!is.null(attr(TE, "pairwise")) ||
       inherits(TE, "pairwise"))) {
    is.pairwise <- TRUE
    #
    sm <- attr(TE, "sm")
    #
    if (missing.reference.group) {
      reference.group <- attr(TE, "reference.group")
      #
      if (is.null(reference.group))
        reference.group <- ""
      else
        avail.reference.group.pairwise <- TRUE
    }
    #
    keep.all.comparisons <- attr(TE, "keep.all.comparisons")
    if (!is.null(keep.all.comparisons) && !keep.all.comparisons)
      stop("First argument is a pairwise object created with ",
           "'keep.all.comparisons = FALSE'.",
           call. = TRUE)
    #
    if (is.null(attr(TE, "varnames")))
      seTE <- TE$seTE
    else
      seTE <- TE[[attr(TE, "varnames")[2]]]
    #
    treat1 <- TE$treat1
    treat2 <- TE$treat2
    studlab <- TE$studlab
    #
    if (!is.null(TE$n1))
      n1 <- TE$n1
    if (!is.null(TE$n2))
      n2 <- TE$n2
    if (!is.null(TE$event1))
      event1 <- TE$event1
    if (!is.null(TE$event2))
      event2 <- TE$event2
    if (!is.null(TE$incr))
      incr <- TE$incr
    if (!is.null(TE$mean1))
      mean1 <- TE$mean1
    if (!is.null(TE$mean2))
      mean2 <- TE$mean2
    if (!is.null(TE$sd1))
      sd1 <- TE$sd1
    if (!is.null(TE$sd2))
      sd2 <- TE$sd2
    if (!is.null(TE$time1))
      time1 <- TE$time1
    if (!is.null(TE$time2))
      time2 <- TE$time2
    #
    pairdata <- TE
    data <- TE
    #
    if (is.null(attr(TE, "varnames")))
      TE <- TE$TE
    else
      TE <- TE[[attr(TE, "varnames")[1]]]
  }
  else {
    is.pairwise <- FALSE
    if (missing(sm))
      if (!is.null(data) && !is.null(attr(data, "sm")))
        sm <- attr(data, "sm")
      else
        sm <- ""
    #
    seTE <- catch("seTE", mc, data, sfsp)
    #
    treat1 <- catch("treat1", mc, data, sfsp)
    treat2 <- catch("treat2", mc, data, sfsp)
    #
    studlab <- catch("studlab", mc, data, sfsp)
    #
    n1 <- catch("n1", mc, data, sfsp)
    n2 <- catch("n2", mc, data, sfsp)
    #
    event1 <- catch("event1", mc, data, sfsp)
    event2 <- catch("event2", mc, data, sfsp)
    #
    incr <- catch("incr", mc, data, sfsp)
    #
    mean1 <- catch("mean1", mc, data, sfsp)
    mean2 <- catch("mean2", mc, data, sfsp)
    #
    sd1 <- catch("sd1", mc, data, sfsp)
    sd2 <- catch("sd2", mc, data, sfsp)
    #
    time1 <- catch("time1", mc, data, sfsp)
    time2 <- catch("time2", mc, data, sfsp)
  }
  #
  chknumeric(TE)
  chknumeric(seTE)
  #
  if (!any(!is.na(TE) & !is.na(seTE)))
    stop("Missing data for estimates (argument 'TE') and ",
         "standard errors (argument 'seTE') in all studies.\n  ",
         "No network meta-analysis possible.",
         call. = FALSE)
  #
  k.Comp <- length(TE)
  #
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  #
  # Remove leading and trailing whitespace
  #
  treat1 <- rmSpace(rmSpace(treat1, end = TRUE))
  treat2 <- rmSpace(rmSpace(treat2, end = TRUE))
  #
  if (length(studlab) == 0) {
    if (warn)
      warning("No information given for argument 'studlab'. ",
              "Assuming that comparisons are from independent studies.",
              call. = FALSE)
    studlab <- seq(along = TE)
  }
  studlab <- as.character(studlab)
  #
  subset <- catch("subset", mc, data, sfsp)
  missing.subset <- is.null(subset)
  #
  correlated <- catch("correlated", mc, data, sfsp)
  if (is.null(correlated))
    correlated <- FALSE
  if (!is.logical(correlated))
    stop("Argument 'correlated' must be a logical vector.", call. = FALSE)
  if (length(correlated) == 1)
    correlated <- rep(correlated, length(TE))
  else if (length(correlated) != length(TE))
    stop("Different length for arguments 'TE' and 'correlated'.", call. = FALSE)
  #
  if (!is.null(event1) & !is.null(event2))
    available.events <- TRUE
  else
    available.events <- FALSE
  #
  if (!is.null(n1) & !is.null(n2))
    available.n <- TRUE
  else
    available.n <- FALSE
  #
  if (available.events & is.null(incr))
    incr <- rep(0, length(event2))
  #
  if (!is.null(mean1) & !is.null(mean2))
    available.means <- TRUE
  else
    available.means <- FALSE
  #
  if (!is.null(sd1) & !is.null(sd2))
    available.sds <- TRUE
  else
    available.sds <- FALSE
  #
  if (!is.null(time1) & !is.null(time2))
    available.times <- TRUE
  else
    available.times <- FALSE
  
  
  #
  #
  # (2b) Store complete dataset in list object data
  #      (if argument keepdata is TRUE)
  #
  #
  if (keepdata) {
    if (nulldata & !is.pairwise)
      data <- data.frame(.studlab = studlab, stringsAsFactors = FALSE)
    else if (nulldata & is.pairwise) {
      data <- pairdata
      data$.studlab <- studlab
    }
    else
      data$.studlab <- studlab
    #
    data$.order <- seq_along(studlab)
    #
    data$.treat1 <- treat1
    data$.treat2 <- treat2
    #
    data$.TE <- TE
    data$.seTE <- seTE
    #
    data$.correlated <- correlated
    #
    data$.event1 <- event1
    data$.n1 <- n1
    data$.event2 <- event2
    data$.n2 <- n2
    data$.incr <- incr
    #
    data$.mean1 <- mean1
    data$.sd1 <- sd1
    data$.mean2 <- mean2
    data$.sd2 <- sd2
    #
    data$.time1 <- time1
    data$.time2 <- time2
    #
    # Check for correct treatment order within comparison
    #
    wo <- data$.treat1 > data$.treat2
    #
    if (any(wo)) {
      data$.TE[wo] <- -data$.TE[wo]
      ttreat1 <- data$.treat1
      data$.treat1[wo] <- data$.treat2[wo]
      data$.treat2[wo] <- ttreat1[wo]
      #
      if (isCol(data, ".n1") & isCol(data, ".n2")) {
        tn1 <- data$.n1
        data$.n1[wo] <- data$.n2[wo]
        data$.n2[wo] <- tn1[wo]
      }
      #
      if (isCol(data, ".event1") & isCol(data, ".event2")) {
        tevent1 <- data$.event1
        data$.event1[wo] <- data$.event2[wo]
        data$.event2[wo] <- tevent1[wo]
      }
      #
      if (isCol(data, ".mean1") & isCol(data, ".mean2")) {
        tmean1 <- data$.mean1
        data$.mean1[wo] <- data$.mean2[wo]
        data$.mean2[wo] <- tmean1[wo]
      }
      #
      if (isCol(data, ".sd1") & isCol(data, ".sd2")) {
        tsd1 <- data$.sd1
        data$.sd1[wo] <- data$.sd2[wo]
        data$.sd2[wo] <- tsd1[wo]
      }
      #
      if (isCol(data, ".time1") & isCol(data, ".time2")) {
        ttime1 <- data$.time1
        data$.time1[wo] <- data$.time2[wo]
        data$.time2[wo] <- ttime1[wo]
      }
    }
    #
    if (!missing.subset) {
      if (length(subset) == dim(data)[1])
        data$.subset <- subset
      else {
        data$.subset <- FALSE
        data$.subset[subset] <- TRUE
      }
    }
  }
  
  
  #
  #
  # (3) Use subset for analysis
  #
  #
  if (!missing.subset) {
    if ((is.logical(subset) & (sum(subset) > k.Comp)) ||
        (length(subset) > k.Comp))
      stop("Length of subset is larger than number of studies.",
           call. = FALSE)
    #
    TE <- TE[subset]
    seTE <- seTE[subset]
    treat1 <- treat1[subset]
    treat2 <- treat2[subset]
    studlab <- studlab[subset]
    #
    correlated <- correlated[subset]
    #
    if (!is.null(n1))
      n1 <- n1[subset]
    if (!is.null(n2))
      n2 <- n2[subset]
    if (!is.null(event1))
      event1 <- event1[subset]
    if (!is.null(event2))
      event2 <- event2[subset]
    if (!is.null(incr))
      incr <- incr[subset]
    if (!is.null(mean1))
      mean1 <- mean1[subset]
    if (!is.null(mean2))
      mean2 <- mean2[subset]
    if (!is.null(sd1))
      sd1 <- sd1[subset]
    if (!is.null(sd2))
      sd2 <- sd2[subset]
    if (!is.null(time1))
      time1 <- time1[subset]
    if (!is.null(time2))
      time2 <- time2[subset]
  }
  #
  labels <- sort(unique(c(treat1, treat2)))
  #
  if (!is.null(seq))
    seq <- setseq(seq, labels)
  else {
    seq <- labels
    if (is.numeric(seq))
      seq <- as.character(seq)
  }
  #
  sep.trts <- setsep(labels, sep.trts, missing = missing.sep.trts)
  
  
  #
  #
  # (4) Additional checks
  #
  #
  if (any(treat1 == treat2))
    stop("Treatments must be different (arguments 'treat1' and 'treat2').",
         call. = FALSE)
  #
  # Check for correct number of comparisons
  #
  tabnarms <- table(studlab)
  sel.narms <- !is_wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  #
  if (sum(sel.narms) == 1)
    stop("Study '", names(tabnarms)[sel.narms],
         "' has a wrong number of comparisons.",
         "\n  Please provide data for all treatment comparisons",
         " (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
         call. = FALSE)
  if (sum(sel.narms) > 1)
    stop("The following studies have a wrong number of comparisons: ",
         paste(paste0("'", names(tabnarms)[sel.narms], "'"),
               collapse = ", "),
         "\n  Please provide data for all treatment comparisons",
         " (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
         call. = FALSE)
  #
  # Check number of subgraphs
  #
  n.subnets <- netconnection(treat1, treat2, studlab)$n.subnets
  #
  if (n.subnets > 1)
    stop("Network consists of ", n.subnets, " separate sub-networks.\n  ",
         "Use R function 'netconnection' to identify sub-networks.",
         call. = FALSE)
  #
  # Check NAs and zero standard errors
  #
  excl <- is.na(TE) | is.na(seTE) | seTE <= 0
  #
  if (any(excl)) {
    if (keepdata)
      data$.excl <- excl
    #
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
      cat("Comparison",
          if (sum(excl) > 1) "s",
          " not considered in network meta-analysis:\n",
          sep = "")
      prmatrix(dat.NAs, quote = FALSE, right = TRUE,
               rowlab = rep("", sum(excl)))
      cat("\n")
    }
    #
    studlab <- studlab[!excl]
    treat1  <- treat1[!excl]
    treat2  <- treat2[!excl]
    TE      <- TE[!excl]
    seTE    <- seTE[!excl]
    #
    correlated <- correlated[!excl]
    #
    if (!is.null(n1))
      n1 <- n1[!excl]
    if (!is.null(n2))
      n2 <- n2[!excl]
    if (!is.null(event1))
      event1 <- event1[!excl]
    if (!is.null(event2))
      event2 <- event2[!excl]
    if (!is.null(incr))
      incr <- incr[!excl]
    if (!is.null(mean1))
      mean1 <- mean1[!excl]
    if (!is.null(mean2))
      mean2 <- mean2[!excl]
    if (!is.null(sd1))
      sd1 <- sd1[!excl]
    if (!is.null(sd2))
      sd2 <- sd2[!excl]
    if (!is.null(time1))
      time1 <- time1[!excl]
    if (!is.null(time2))
      time2 <- time2[!excl]
    #
    seq <- seq[seq %in% unique(c(treat1, treat2))]
    labels <- labels[labels %in% unique(c(treat1, treat2))]
  }
  #
  # Check for correct number of comparisons (after removing
  # comparisons with missing data)
  #
  tabnarms <- table(studlab)
  sel.narms <- !is_wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  #
  if (sum(sel.narms) == 1)
    stop("After removing comparisons with missing treatment effects",
         " or standard errors,\n  study '",
         names(tabnarms)[sel.narms],
         "' has a wrong number of comparisons.",
         " Please check data and\n  consider to remove study",
         " from network meta-analysis.",
         call. = FALSE)
  if (sum(sel.narms) > 1)
    stop("After removing comparisons with missing treatment effects",
         " or standard errors,\n  the following studies have",
         " a wrong number of comparisons: ",
         paste(paste0("'", names(tabnarms)[sel.narms], "'"), collapse = ", "),
         "\n  Please check data and consider to remove studies",
         " from network meta-analysis.",
         call. = FALSE)
  #
  # Check number of subgraphs
  #
  n.subnets <- netconnection(treat1, treat2, studlab)$n.subnets
  #
  if (n.subnets > 1)
    stop("After removing comparisons with missing treatment effects",
         " or standard errors,\n  network consists of ",
         n.subnets, " separate sub-networks.\n  ",
         "Please check data and consider to remove studies",
         " from network meta-analysis.",
         call. = FALSE)
  #
  # Check for correct treatment order within comparison
  #
  wo <- treat1 > treat2
  #
  if (any(wo)) {
    TE[wo] <- -TE[wo]
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
    #
    if (available.n) {
      tn1 <- n1
      n1[wo] <- n2[wo]
      n2[wo] <- tn1[wo]
    }
    #
    if (available.events) {
      tevent1 <- event1
      event1[wo] <- event2[wo]
      event2[wo] <- tevent1[wo]
    }
    #
    if (available.means) {
      tmean1 <- mean1
      mean1[wo] <- mean2[wo]
      mean2[wo] <- tmean1[wo]
    }
    #
    if (available.sds) {
      tsd1 <- sd1
      sd1[wo] <- sd2[wo]
      sd2[wo] <- tsd1[wo]
    }
    #
    if (available.times) {
      ttime1 <- time1
      time1[wo] <- time2[wo]
      time2[wo] <- ttime1[wo]
    }
  }
  #
  # Set reference group
  #
  if (missing.reference.group & !avail.reference.group.pairwise) {
    go.on <- TRUE
    i <- 0
    while (go.on) {
      i <- i + 1
      sel.i <-
        !is.na(TE) & !is.na(seTE) &
        (treat1 == labels[i] | treat2 == labels[i])
      if (sum(sel.i) > 0) {
        go.on <- FALSE
        reference.group <- labels[i]
      }
      else if (i == length(labels)) {
        go.on <- FALSE
        reference.group <- ""
      }
    }
  }
  #
  if (is.null(all.treatments))
    if (reference.group == "")
      all.treatments <- TRUE
    else
      all.treatments <- FALSE
  #
  # Check reference group
  #
  if (reference.group != "")
    reference.group <- setref(reference.group, labels)
  
  
  #
  #
  # (5) Generate analysis dataset
  #
  #
  #
  # Calculate weight matrix and generate ordered dataset, with added numbers
  # of arms per study
  #
  p0 <- prepare2(TE, seTE, treat1, treat2, studlab, correlated = correlated,
                 func.inverse = func.inverse)
  #
  W.matrix.common <- p0$W
  dat.c <- p0$data
  #
  # Check consistency of treatment effects and standard errors in
  # multi-arm studies
  #
  chkmultiarm(dat.c$TE, dat.c$seTE, dat.c$treat1, dat.c$treat2, dat.c$studlab,
              dat.c$correlated,
              tol.multiarm = tol.multiarm, tol.multiarm.se = tol.multiarm.se,
              details = details.chkmultiarm)
  #
  # Study overview
  #
  tdata <- data.frame(studies = dat.c$studlab, narms = dat.c$narms,
                      order = dat.c$order,
                      stringsAsFactors = FALSE)
  #
  tdata <- tdata[!duplicated(tdata[, c("studies", "narms")]), , drop = FALSE]
  studies <- tdata$studies[order(tdata$order)]
  narms <- tdata$narms[order(tdata$order)]
  
  
  #
  #
  # (6) Conduct network meta-analysis
  #
  
  # Common effects model
  #
  res.c <- nma_ruecker(dat.c$TE,
                       as.matrix(W.matrix.common),
                       sqrt(1 / dat.c$weights),
                       dat.c$treat1, dat.c$treat2,
                       dat.c$treat1.pos, dat.c$treat2.pos,
                       dat.c$narms, dat.c$studlab,
                       sm,
                       level, level.ma,
                       dat.c$seTE, 0, sep.trts,
                       method.tau,
                       func.inverse, Cov0 = p0$Cov)
  trts <- rownames(res.c$A.matrix)
  #
  #
  # Random effects model
  #
  if (is.null(tau.preset)) {
    if (method.tau %in% c("ML", "REML")) {
      #
      dat.tau <-
        data.frame(studlab = studlab,
                   treat1 = treat1, treat2 = treat2,
                   TE = TE, seTE = seTE)
      if (available.n) {
        dat.tau$n1 <- n1
        dat.tau$n2 <- n2
      }
      if (available.events) {
        dat.tau$event1 <- event1
        dat.tau$event2 <- event2
        dat.tau$incr <- incr
      }
      if (available.means) {
        dat.tau$mean1 <- mean1
        dat.tau$mean2 <- mean2
      }
      if (available.sds) {
        dat.tau$sd1 <- sd1
        dat.tau$sd2 <- sd2
      }
      if (available.times) {
        dat.tau$time1 <- time1
        dat.tau$time2 <- time2
      }
      #
      #dat.tau <- dat.tau[order(dat.tau$studlab,
      #                         dat.tau$treat1, dat.tau$treat2), , drop = FALSE]
      #
      keep <- logical(0)
      wo <- logical(0)
      #
      for (i in unique(dat.tau$studlab)) {
        d.i <- dat.tau[dat.tau$studlab == i, , drop = FALSE]
        trts.i <- unique(sort(c(d.i$treat1, d.i$treat2)))
        if (reference.group %in% trts.i)
          ref.i <- reference.group
        else
          ref.i <- rev(trts.i)[1]
        #
        keep.i <- !(d.i$treat1 != ref.i & d.i$treat2 != ref.i)
        wo.i <- d.i$treat1 == ref.i
        #
        keep <- c(keep, keep.i)
        wo <- c(wo, wo.i)
      }
      #
      dat.tau <- dat.tau[keep, , drop = FALSE]
      #
      wo <- wo[keep]
      if (sum(wo) > 0) {
        dat.tau$TE[wo] <- -dat.tau$TE[wo]
        #
        t2.i <- dat.tau$treat2
        e2.i <- dat.tau$event2
        n2.i <- dat.tau$n2
        mean2.i <- dat.tau$mean2
        sd2.i <- dat.tau$sd2
        time2.i <- dat.tau$time2
        #
        dat.tau$treat2[wo] <- dat.tau$treat1[wo]
        dat.tau$event2[wo] <- dat.tau$event1[wo]
        dat.tau$n2[wo] <- dat.tau$n1[wo]
        dat.tau$mean2[wo] <- dat.tau$mean1[wo]
        dat.tau$sd2[wo] <- dat.tau$sd1[wo]
        dat.tau$time2[wo] <- dat.tau$time2[wo]
        #
        dat.tau$treat1[wo] <- t2.i[wo]
        dat.tau$event1[wo] <- e2.i[wo]
        dat.tau$n1[wo] <- n2.i[wo]
        dat.tau$mean1[wo] <- mean2.i[wo]
        dat.tau$sd1[wo] <- sd2.i[wo]
        dat.tau$time1[wo] <- time2.i[wo]
      }
      #
      ncols1 <- ncol(dat.tau)
      dat.tau <- contrmat(dat.tau, grp1 = "treat1", grp2 = "treat2")
      ncols2 <- ncol(dat.tau)
      newnames <- paste0("V", seq_len(ncols2 - ncols1))
      names(dat.tau)[(ncols1 + 1):ncols2] <- newnames
      #
      dat.tau <- dat.tau[order(dat.tau$studlab), ]
      #
      trts.tau <- newnames[-length(newnames)]
      #
      formula.trts <-
        as.formula(paste("~ ", paste(trts.tau, collapse = " + "), " - 1"))
      #
      # Calculate Variance-Covariance matrix
      #
      if (available.n &
          (available.events | available.times | (available.sds))) {
        V <- bldiag(lapply(split(dat.tau, dat.tau$studlab), calcV, sm = sm))
      }
      else
        V <- dat.tau$seTE^2
      #
      dat.tau.TE <- dat.tau$TE
      dat.tau$comparison <- paste(dat.tau$treat1, dat.tau$treat2, sep = " vs ")
      #
      if (length(dat.tau.TE) == 1) {
        rma1 <- runNN(rma.uni,
                      list(yi = dat.tau.TE, vi = V, data = dat.tau,
                           method = method.tau, control = control))
        #
        tau <- NA
      }
      else {
        rma1 <- runNN(rma.mv,
                      list(yi = dat.tau.TE, V = V,
                           data = dat.tau,
                           mods = formula.trts,
                           random = as.call(~ factor(comparison) | studlab),
                           rho = 0.5,
                           method = method.tau, control = control))
        #
        tau <- sqrt(rma1$tau2)
      }
    }
    else
      tau <- res.c$tau
  }
  else
    tau <- tau.preset
  #
  p1 <- prepare2(TE, seTE, treat1, treat2, studlab, tau, correlated,
                 func.inverse)
  #
  W.matrix.random <- p1$W
  dat.r <- p1$data
  #
  res.r <- nma_ruecker(dat.r$TE,
                       as.matrix(W.matrix.random),
                       sqrt(1 / dat.r$weights),
                       dat.r$treat1, dat.r$treat2,
                       dat.r$treat1.pos, dat.r$treat2.pos,
                       dat.r$narms, dat.r$studlab,
                       sm,
                       level, level.ma,
                       dat.r$seTE, tau, sep.trts,
                       method.tau,
                       func.inverse, Cov0 = p1$Cov)
  #
  TE.random <- res.r$TE.pooled
  seTE.random <- res.r$seTE.pooled
  df.Q <- res.c$df
  #
  # Prediction intervals
  #
  if (df.Q == 0)
    prediction <- FALSE
  #
  if (df.Q >= 2) {
    seTE.predict <- sqrt(seTE.random^2 + tau^2)
    ci.p <- ci(TE.random, seTE.predict, level.predict, df.Q)
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


  #
  #
  # (7) Generate R object
  #
  #
  o <- order(dat.c$order)
  #
  designs <- designs(res.c$treat1, res.c$treat2, res.c$studlab,
                     sep.trts = sep.trts)
  #
  W.matrix.common <- W.matrix.common[o, o, drop = FALSE]
  rownames(W.matrix.common) <- colnames(W.matrix.common) <- res.c$studlab[o]
  #
  W.matrix.random <- W.matrix.random[o, o, drop = FALSE]
  rownames(W.matrix.random) <- colnames(W.matrix.random) <- res.c$studlab[o]
  #
  res <- list(studlab = res.c$studlab[o],
              treat1 = res.c$treat1[o],
              treat2 = res.c$treat2[o],
              #
              TE = res.c$TE[o],
              seTE = res.c$seTE.orig[o],
              seTE.adj = res.c$seTE[o],
              seTE.adj.common = res.c$seTE[o],
              seTE.adj.random = res.r$seTE[o],
              correlated = correlated,
              #
              design = designs$design[o],
              #
              event1 = event1,
              event2 = event2,
              n1 = n1,
              n2 = n2,
              incr = incr,
              #
              mean1 = mean1,
              mean2 = mean2,
              sd1 = sd1,
              sd2 = sd2,
              #
              time1 = time1,
              time2 = time2,
              #
              k = res.c$k,
              m = res.c$m,
              n = res.c$n,
              d = length(unique(designs$design)),
              #
              trts = trts,
              k.trts = NA,
              n.trts = if (available.n) NA else NULL,
              events.trts = if (available.events) NA else NULL,
              #
              n.arms = NA,
              multiarm = NA,
              #
              studies = studies,
              narms = narms,
              #
              designs = unique(sort(designs$design)),
              comparisons = "",
              #
              TE.nma.common = res.c$TE.nma[o],
              seTE.nma.common = res.c$seTE.nma[o],
              lower.nma.common = res.c$lower.nma[o],
              upper.nma.common = res.c$upper.nma[o],
              statistic.nma.common = res.c$statistic.nma[o],
              pval.nma.common = res.c$pval.nma[o],
              #
              leverage.common = res.c$leverage[o],
              w.common = res.c$w.pooled[o],
              Q.common = res.c$Q.pooled[o],
              #
              TE.common = res.c$TE.pooled,
              seTE.common = res.c$seTE.pooled,
              lower.common = res.c$lower.pooled,
              upper.common = res.c$upper.pooled,
              statistic.common = res.c$statistic.pooled,
              pval.common = res.c$pval.pooled,
              #
              TE.nma.random = res.r$TE.nma[o],
              seTE.nma.random = res.r$seTE.nma[o],
              lower.nma.random = res.r$lower.nma[o],
              upper.nma.random = res.r$upper.nma[o],
              statistic.nma.random = res.r$statistic.nma[o],
              pval.nma.random = res.r$pval.nma[o],
              #
              w.random = res.r$w.pooled[o],
              #
              TE.random = TE.random,
              seTE.random = seTE.random,
              lower.random = res.r$lower.pooled,
              upper.random = res.r$upper.pooled,
              statistic.random = res.r$statistic.pooled,
              pval.random = res.r$pval.pooled,
              #
              seTE.predict = seTE.predict,
              lower.predict = p.lower,
              upper.predict = p.upper,
              method.predict = "V",
              #
              prop.direct.common = NA,
              prop.direct.random = NA,
              #
              TE.direct.common = res.c$TE.direct,
              seTE.direct.common = res.c$seTE.direct,
              lower.direct.common = res.c$lower.direct,
              upper.direct.common = res.c$upper.direct,
              statistic.direct.common = res.c$statistic.direct,
              pval.direct.common = res.c$pval.direct,
              #
              TE.direct.random = res.r$TE.direct,
              seTE.direct.random = res.r$seTE.direct,
              lower.direct.random = res.r$lower.direct,
              upper.direct.random = res.r$upper.direct,
              statistic.direct.random = res.r$statistic.direct,
              pval.direct.random = res.r$pval.direct,
              #
              Q.direct = res.r$Q.direct,
              tau.direct = sqrt(res.r$tau2.direct),
              tau2.direct = res.r$tau2.direct,
              I2.direct = res.r$I2.direct,
              #
              TE.indirect.common = NA,
              seTE.indirect.common = NA,
              lower.indirect.common = NA,
              upper.indirect.common = NA,
              statistic.indirect.common = NA,
              pval.indirect.common = NA,
              #
              TE.indirect.random = NA,
              seTE.indirect.random = NA,
              lower.indirect.random = NA,
              upper.indirect.random = NA,
              statistic.indirect.random = NA,
              pval.indirect.random = NA,
              #
              Q = res.c$Q,
              df.Q = df.Q,
              pval.Q = res.c$pval.Q,
              I2 = res.c$I2,
              lower.I2 = res.c$lower.I2,
              upper.I2 = res.c$upper.I2,
              tau = tau,
              tau2 = tau^2,
              #
              Q.heterogeneity = NA,
              df.Q.heterogeneity = NA,
              pval.Q.heterogeneity = NA,
              #
              Q.inconsistency = NA,
              df.Q.inconsistency = NA,
              pval.Q.inconsistency = NA,
              #
              Q.decomp = res.c$Q.decomp,
              #
              W.matrix.common = W.matrix.common,
              W.matrix.random = W.matrix.random,
              #
              A.matrix = res.c$A.matrix,
              X.matrix = res.c$B.matrix[o, ],
              B.matrix = res.c$B.matrix[o, ],
              #
              L.matrix.common = res.c$L.matrix,
              Lplus.matrix.common = res.c$Lplus.matrix,
              L.matrix.random = res.r$L.matrix,
              Lplus.matrix.random = res.r$Lplus.matrix,
              #
              Q.matrix = res.c$Q.matrix,
              #
              G.matrix = res.c$G.matrix[o, o, drop = FALSE],
              #
              H.matrix.common = res.c$H.matrix[o, o, drop = FALSE],
              H.matrix.random = res.r$H.matrix[o, o, drop = FALSE],
              #
              n.matrix = if (available.n) NA else NULL,
              events.matrix = if (available.events) NA else NULL,
              #
              P.common = NA,
              P.random = NA,
              #
              Cov.common = res.c$Cov,
              Cov.random = res.r$Cov,
              #
              treat1.pos = res.c$treat1.pos[o],
              treat2.pos = res.c$treat2.pos[o],
              #
              sm = sm,
              method = "Inverse",
              level = level,
              level.ma = level.ma,
              common = common,
              random = random,
              #
              prediction = prediction,
              level.predict = level.predict,
              #
              reference.group = reference.group,
              baseline.reference = baseline.reference,
              all.treatments = all.treatments,
              seq = seq,
              #
              method.tau = method.tau,
              tau.preset = tau.preset,
              #
              tol.multiarm = tol.multiarm,
              tol.multiarm.se = tol.multiarm.se,
              details.chkmultiarm = details.chkmultiarm,
              #
              func.inverse = deparse(substitute(func.inverse)),
              #
              sep.trts = sep.trts,
              nchar.trts = nchar.trts,
              nchar.studlab = nchar.studlab,
              #
              overall.hetstat = overall.hetstat,
              backtransf = backtransf,
              #
              title = title,
              #
              warn = warn,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  #
  class(res) <- "netmeta"
  #
  # Add results for indirect treatment estimates
  #
  n <- res$n
  #
  res$prop.direct.common <-
    netmeasures(res, random = FALSE, warn = warn)$proportion
  # Print warning(s) in call of netmeasures() once
  res$prop.direct.random <-
    suppressWarnings(netmeasures(res, random = TRUE,
                                 tau.preset = res$tau,
                                 warn = FALSE)$proportion)
  if (is.logical(res$prop.direct.common))
    res$prop.direct.common <- as.numeric(res$prop.direct.common)
  if (is.logical(res$prop.direct.random))
    res$prop.direct.random <- as.numeric(res$prop.direct.random)
  #
  res$comparisons <-
    names(res$prop.direct.random)[!is_zero(res$prop.direct.random)]
  #
  # Add P.common and P.random
  #
  P.common <- P.random <- matrix(NA, n, n)
  colnames(P.common) <- rownames(P.common) <-
    colnames(P.random) <- rownames(P.random) <- trts
  #
  if (n == 2) {
    #
    # For two treatments only direct evidence is available
    #
    res$prop.direct.common <- 1
    res$prop.direct.random <- 1
    names(res$prop.direct.common) <-
      names(res$prop.direct.random) <- paste(labels, collapse = sep.trts)
    #
    sel <- row(P.common) != col(P.common)
    P.common[sel] <- 1
    P.random[sel] <- 1
  }
  else {
    k <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        k <- k + 1
        P.common[i, j] <- P.common[j, i] <- res$prop.direct.common[k]
        P.random[i, j] <- P.random[j, i] <- res$prop.direct.random[k]
      }
    }
  }
  #
  # Set direct evidence estimates to 0 if only indirect evidence is available
  # (otherwise indirect estimates would be NA as direct estimates are NA)
  #
  TE.direct.common <- res$TE.direct.common
  TE.direct.random <- res$TE.direct.random
  #
  TE.direct.common[abs(P.common) < .Machine$double.eps^0.5] <- 0
  TE.direct.random[abs(P.random) < .Machine$double.eps^0.5] <- 0
  #
  # Indirect estimate is NA if only direct evidence is available
  #
  res$P.common <- P.common
  res$P.random <- P.random
  #
  P.common[abs(P.common - 1) < .Machine$double.eps^0.5] <- NA
  P.common[P.common > 1] <- NA
  P.random[abs(P.random - 1) < .Machine$double.eps^0.5] <- NA
  P.random[P.random > 1] <- NA
  #
  # Common effects model
  #
  ci.if <- ci((res$TE.common - P.common * TE.direct.common) / (1 - P.common),
              sqrt(res$seTE.common^2 / (1 - P.common)),
              level = level)
  #
  res$TE.indirect.common   <- ci.if$TE
  res$seTE.indirect.common <- ci.if$seTE
  #
  res$lower.indirect.common <- ci.if$lower
  res$upper.indirect.common <- ci.if$upper
  #
  res$statistic.indirect.common <- ci.if$statistic
  res$pval.indirect.common <- ci.if$p
  #
  # Random effects model
  #
  ci.ir <- ci((res$TE.random - P.random * TE.direct.random) / (1 - P.random),
              sqrt(res$seTE.random^2 / (1 - P.random)),
              level = level)
  #
  res$TE.indirect.random   <- ci.ir$TE
  res$seTE.indirect.random <- ci.ir$seTE
  #
  res$lower.indirect.random <- ci.ir$lower
  res$upper.indirect.random <- ci.ir$upper
  #
  res$statistic.indirect.random <- ci.ir$statistic
  res$pval.indirect.random <- ci.ir$p
  #
  # Additional assignments
  #
  res$small.values <- small.values
  #
  l1 <- length(res$treat1)
  tab.trts <-
    table(longarm(res$treat1, res$treat2,
                  rep(1, l1), rep(2, l1), rep(1, l1), rep(2, l1),
                  studlab = res$studlab)$treat)
  res$k.trts <- as.numeric(tab.trts)
  names(res$k.trts) <- names(tab.trts)
  #
  if (any(res$narms > 2)) {
    tdata1 <- data.frame(studlab = res$studlab,
                         .order = seq(along = res$studlab))
    tdata2 <- data.frame(studlab = as.character(res$studies),
                         narms = res$narms)
    #
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
  #
  # Set leverage of multi-arm studies to NA
  #
  if (any(res$multiarm))
    res$leverage.common[res$multiarm] <- NA
  
  
  #
  # Calculate heterogeneity and inconsistency statistics
  #
  if (res$d > 1) {
    dd <- decomp.design(res, warn = FALSE)
    res$Q.heterogeneity <- dd$Q.decomp$Q[2]
    res$Q.inconsistency <- dd$Q.decomp$Q[3]
    #
    res$df.Q.heterogeneity <- dd$Q.decomp$df[2]
    res$df.Q.inconsistency <- dd$Q.decomp$df[3]
    #
    res$pval.Q.heterogeneity <- dd$Q.decomp$pval[2]
    res$pval.Q.inconsistency <- dd$Q.decomp$pval[3]
  }
  
  
  if (keepdata) {
    data$.design <- designs(data$.treat1, data$.treat2, data$.studlab,
                            sep = sep.trts)$design
    #
    res$data <- merge(data,
                      data.frame(.studlab = res$studies,
                                 .narms = res$narms),
                      by = ".studlab",
                      stringsAsFactors = FALSE)
    #
    # Store adjusted standard errors in dataset
    #
    res$data <- merge(res$data,
                      data.frame(.studlab = res.c$studlab,
                                 .treat1 = res.c$treat1,
                                 .treat2 = res.c$treat2,
                                 .seTE.adj.common = res.c$seTE),
                      by = c(".studlab", ".treat1", ".treat2"),
                      stringsAsFactors = FALSE)
    #
    res$data <- merge(res$data,
                      data.frame(.studlab = res.r$studlab,
                                 .treat1 = res.r$treat1,
                                 .treat2 = res.r$treat2,
                                 .seTE.adj.random = res.r$seTE),
                      by = c(".studlab", ".treat1", ".treat2"),
                      stringsAsFactors = FALSE)
    #
    res$data <- res$data[order(res$data$.order), ]
    res$data$.order <- NULL
  }
  #
  if (available.n) {
    res$n.matrix <- netmatrix(res, n1 + n2, func = "sum")
    #
    dat.n <- data.frame(studlab = c(studlab, studlab),
                        treat = c(treat1, treat2),
                        n = c(n1, n2))
    dat.n <- dat.n[!duplicated(dat.n[, c("studlab", "treat")]), ]
    dat.n <- by(dat.n$n, dat.n$treat, sum, na.rm = TRUE)
    res$n.trts <- as.vector(dat.n[trts])
    names(res$n.trts) <- trts
  }
  #
  if (available.events) {
    res$events.matrix <- netmatrix(res, event1 + event2, func = "sum")
    #
    dat.e <- data.frame(studlab = c(studlab, studlab),
                        treat = c(treat1, treat2),
                        n = c(event1, event2))
    dat.e <- dat.e[!duplicated(dat.e[, c("studlab", "treat")]), ]
    dat.e <- by(dat.e$n, dat.e$treat, sum, na.rm = TRUE)
    res$events.trts <- as.vector(dat.e[trts])
    names(res$events.trts) <- trts
  }
  
  #
  # Add results from rma.mv()
  #
  if (keeprma & method.tau %in% c("ML", "REML"))
    res$rma.tau <- rma1
  
  #
  # Drop list element 'correlated' if no correlated outcomes are considered 
  #
  if (all(!correlated)) {
    res$correlated <- NULL
    res$data <- res$data[, names(res$data) != ".correlated"]
  }
  
  #
  # Backward compatibility
  #
  res$fixed <- res$common
  res$comb.fixed <- res$common
  res$comb.random <- res$random
  #
  res$seTE.adj.fixed <- res$seTE.adj.common
  res$TE.nma.fixed <- res$TE.nma.common
  res$seTE.nma.fixed <- res$seTE.nma.common
  res$lower.nma.fixed <- res$lower.nma.common
  res$upper.nma.fixed <- res$upper.nma.common
  res$statistic.nma.fixed <- res$statistic.nma.common
  res$pval.nma.fixed <- res$pval.nma.common
  res$leverage.fixed <- res$leverage.common
  res$w.fixed <- res$w.common
  res$Q.fixed <- res$Q.common
  #
  res$TE.fixed <- res$TE.common
  res$seTE.fixed <- res$seTE.common
  res$lower.fixed <- res$lower.common
  res$upper.fixed <- res$upper.common
  res$statistic.fixed <- res$statistic.common
  res$pval.fixed <- res$pval.common
  #
  res$prop.direct.fixed <- res$prop.direct.common
  #
  res$TE.direct.fixed <- res$TE.direct.common
  res$seTE.direct.fixed <- res$seTE.direct.common
  res$lower.direct.fixed <- res$lower.direct.common
  res$upper.direct.fixed <- res$upper.direct.common
  res$statistic.direct.fixed <- res$statistic.direct.common
  res$pval.direct.fixed <- res$pval.direct.common
  #
  res$TE.indirect.fixed <- res$TE.indirect.common
  res$seTE.indirect.fixed <- res$seTE.indirect.common
  res$lower.indirect.fixed <- res$lower.indirect.common
  res$upper.indirect.fixed <- res$upper.indirect.common
  res$statistic.indirect.fixed <- res$statistic.indirect.common
  res$pval.indirect.fixed <- res$pval.indirect.common
  #
  res$L.matrix.fixed <- res$L.matrix.common
  res$Lplus.matrix.fixed <- res$Lplus.matrix.common
  res$H.matrix.fixed <- res$H.matrix.common
  res$P.fixed <- res$P.common
  res$Cov.fixed <- res$Cov.common
  
  res
}
