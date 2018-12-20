\name{netmeta}
\alias{netmeta}

\title{Network meta-analysis using graph-theoretical method}

\description{
  Network meta-analysis is a generalisation of pairwise meta-analysis
  that compares all pairs of treatments within a number of treatments
  for the same condition. The graph-theoretical method for analysis of
  network meta-analyses uses graph-theoretical methods that were
  originally developed in electrical network theory. It has been found
  to be equivalent to the frequentist approach to network
  meta-analysis (Rücker, 2012).
}

\usage{
netmeta(TE, seTE, treat1, treat2, studlab, data=NULL, subset=NULL,
        sm, level=gs("level"), level.comb=gs("level.comb"),
        comb.fixed=gs("comb.fixed"),
        comb.random=gs("comb.random") | !is.null(tau.preset),
        prediction=FALSE, level.predict=gs("level.predict"),
        reference.group="", baseline.reference=TRUE,
        all.treatments=NULL, seq=NULL, tau.preset=NULL,
        tol.multiarm = 0.0005, details.chkmultiarm = FALSE,
        sep.trts=":", nchar.trts=666,
        n1=NULL, n2=NULL, event1=NULL, event2=NULL,
        backtransf=gs("backtransf"), title="",
        keepdata=gs("keepdata"), warn=TRUE)
}

\arguments{
  \item{TE}{Estimate of treatment effect, i.e. difference between
   first and second treatment (e.g. log odds ratio, mean difference,
   or log hazard ratio).}
  \item{seTE}{Standard error of treatment estimate.}
  \item{treat1}{Label/Number for first treatment.}
  \item{treat2}{Label/Number for second treatment.}
  \item{studlab}{An optional - but important! - vector with study labels (see Details).}
  \item{data}{An optional data frame containing the study information.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{sm}{A character string indicating underlying summary measure,
    e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
    \code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.}
  \item{level}{The level used to calculate confidence intervals for
    individual comparisons.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
  \item{prediction}{A logical indicating whether prediction intervals
    should be printed.}
  \item{level.predict}{The level used to calculate prediction
    intervals for a new study.}
  \item{reference.group}{Reference group.}
  \item{baseline.reference}{A logical indicating whether results
    should be expressed as comparisons of other treatments versus the
    reference treatment (default) or vice versa. This argument is only
    considered if \code{reference.group} has been specified.}
  \item{all.treatments}{A logical or value \code{"NULL"}. If
    \code{TRUE}, matrices with all treatment effects, and confidence
    limits will be printed.}
  \item{seq}{A character or numerical vector specifying the sequence
    of treatments in printouts.}
  \item{tau.preset}{An optional value for the square-root of the
    between-study variance \eqn{\tau^2}.}
  \item{tol.multiarm}{A numeric for the tolerance for consistency of
    treatment estimates and corresponding variances in multi-arm studies
    which are consistent by design.}
  \item{details.chkmultiarm}{A logical indicating whether treatment
    estimates and / or variances of multi-arm studies with inconsistent
    results or negative multi-arm variances should be printed.}
  \item{sep.trts}{A character used in comparison names as separator
    between treatment labels.}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in printouts and forest plots. If
    \code{backtransf=TRUE}, results for \code{sm="OR"} are presented
    as odds ratios rather than log odds ratios, for example.}
  \item{nchar.trts}{A numeric defining the minium number of characters
    used to create unique treatment names (see Details).}
  \item{n1}{Number of observations in first treatment group.}
  \item{n2}{Number of observations in second treatment group.}
  \item{event1}{Number of events in first treatment group.}
  \item{event2}{Number of events in second treatment group.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{keepdata}{A logical indicating whether original data (set)
    should be kept in netmeta object.}
  \item{warn}{A logical indicating whether warnings should be printed
    (e.g., if studies are excluded from meta-analysis due to zero
    standard errors).}
}

\details{
  Network meta-analysis using R package \bold{netmeta} is described in
  detail in Schwarzer et al. (2015), Chapter 8.
  
  Let \emph{n} be the number of different treatments (nodes, vertices)
  in a network and let \emph{m} be the number of existing comparisons
  (edges) between the treatments. If there are only two-arm studies,
  \emph{m} is the number of studies. Let TE and seTE be the vectors of
  observed effects and their standard errors. Let W be the
  \emph{m}x\emph{m} diagonal matrix that contains the inverse variance
  1/seTE^2.
  
  The given comparisons define the network structure. Therefrom an
  \emph{m}x\emph{n} design matrix B (edge-vertex incidence matrix) is
  formed; for more precise information, see Rücker (2012). Moreover,
  the \emph{n}x\emph{n} Laplacian matrix L and its Moore-Penrose
  pseudoinverse L+ are calculated (both matrices play an important
  role in graph theory and electrical network theory).  Using these
  matrices, the variances based on both direct and indirect
  comparisons can be estimated. Moreover, the hat matrix H can be
  estimated by \strong{H = BL+B^{t}W = B(B^{t} W B)^{+}B^{t}W} and
  finally consistent treatment effects can be estimated by applying
  the hat matrix to the observed (potentially inconsistent) effects. H
  is a projection matrix which maps the observed effects onto the
  consistent (n-1)-dimensional subspace. This is the Aitken estimator
  (Senn et al., 2013). As in pairwise meta-analysis, the Q statistic
  measures the deviation from consistency. Q can be separated into
  parts for each pairwise meta-analysis and a part for remaining
  inconsistency between comparisons.

  Often multi-arm studies are included in a network meta-analysis.  In
  multi-arm studies, the treatment effects on different comparisons
  are not independent, but correlated. This is accounted for by
  reweighting all comparisons of each multi-arm study. The method is
  described in Rücker (2012) and Rücker and Schwarzer (2014).

  Comparisons belonging to multi-arm studies are identified by
  identical study labels (argument \code{studlab}). It is therefore
  important to use identical study labels for all comparisons
  belonging to the same multi-arm study, e.g., study label
  "Willms1999" for the three-arm study in the data example (Senn et
  al., 2013). The function netmeta then automatically accounts for
  within-study correlation by reweighting all comparisons of each
  multi-arm study.

  Data entry for this function is in \emph{contrast-based} format,
  that is, data are given as contrasts (differences) between two
  treatments (argument \code{TE}) with standard error (argument
  \code{seTE}). In principle, meta-analysis functions from R package
  \bold{meta}, e.g. \code{\link{metabin}} for binary outcomes or
  \code{\link{metacont}} for continuous outcomes, can be used to
  calculate treatment effects separately for each treatment comparison
  which is a rather tedious enterprise. If data are provided in
  \emph{arm-based} format, that is, data are given for each treatment
  arm separately (e.g. number of events and participants for binary
  outcomes), a much more convenient way to transform data into
  contrast-based form is available. Function \code{\link{pairwise}}
  can automatically transform data with binary outcomes (using the
  \code{\link{metabin}} function from R package \bold{meta}),
  continuous outcomes (\code{\link{metacont}} function), incidence
  rates (\code{\link{metainc}} function), and generic outcomes
  (\code{\link{metagen}} function). Additional arguments of these
  functions can be provided, e.g., to calculate Hedges' \emph{g} or
  Cohen's \emph{d} for continuous outcomes (see help page of function
  \code{\link{pairwise}}).
  
  Note, all pairwise comparisons must be provided for a multi-arm
  study. Consider a multi-arm study of \emph{p} treatments with known
  variances.  For this study, treatment effects and standard errors
  must be provided for each of \emph{p}(\emph{p} - 1)/2 possible
  comparisons. For instance, a three-arm study contributes three
  pairwise comparisons, a four-arm study even six pairwise
  comparisons. Function \code{\link{pairwise}} automatically
  calculates all pairwise comparisons for multi-arm studies.
  
  A simple random effects model assuming that a constant heterogeneity
  variance is added to each comparison of the network can be defined
  via a generalised methods of moments estimate of the between-studies
  variance tau^2 (Jackson et al., 2012). This is added to the observed
  sampling variance seTE^2 of each comparison in the network (before
  appropriate adjustment for multi-arm studies). Then, as in standard
  pairwise meta-analysis, the procedure is repeated with the resulting
  enlarged standard errors.
  
  Internally, both fixed effects and random effects models are
  calculated regardless of values choosen for arguments
  \code{comb.fixed} and \code{comb.random}. Accordingly, the network
  estimates for the random effects model can be extracted from
  component \code{TE.random} of an object of class \code{"netmeta"}
  even if argument \code{comb.random=FALSE}. However, all functions in
  R package \bold{netmeta} will adequately consider the values for
  \code{comb.fixed} and \code{comb.random}. E.g. function
  \code{\link{print.summary.netmeta}} will not print results for the
  random effects model if \code{comb.random=FALSE}.

  By default, treatment names are not abbreviated in
  printouts. However, in order to get more concise printouts, argument
  \code{nchar.trts} can be used to define the minimum number of
  characters for abbreviated treatment names. R function
  \code{\link{treats}} is utilised internally to create abbreviated
  treatment names.

  Names of treatment comparisons are created by concatenating
  treatment labels of pairwise comparisons using \code{sep.trts} as
  separator (see \code{\link{paste}}). These comparison names are used
  in the covariance matrices \code{Cov.fixed} and \code{Cov.random}
  and in some R functions, e.g, \code{\link{decomp.design}}. By
  default, a colon is used as the separator. If any treatment label
  contains a colon the following characters are used as separator (in
  consecutive order): \code{"-"}, \code{"_"}, \code{"/"}, \code{"+"},
  \code{"."}, \code{"|"}, and \code{"*"}.  If all of these characters
  are used in treatment labels, a corresponding error message is
  printed asking the user to specify a different separator.
}

\value{
  An object of class \code{netmeta} with corresponding \code{print},
  \code{summary}, \code{forest}, and \code{netrank} functions. The
  object is a list containing the following components:
  \item{studlab, treat1, treat2, TE, seTE}{As defined above.}
  \item{seTE.adj}{Standard error of treatment estimate, adjusted for
    multi-arm studies.}
  \item{n1, n2, event1, event2}{As defined above.}
  \item{k}{Total number of studies.}
  \item{m}{Total number of pairwise comparisons.}
  \item{n}{Total number of treatments.}
  \item{d}{Total number of designs.}
  \item{trts}{Treatments included in network meta-analysis.}
  \item{k.trts}{Number of studies evaluating a treatment.}
  \item{n.trts}{Number of observations receiving a treatment (if
    arguments \code{n1} and \code{n2} are provided).}
  \item{events.trts}{Number of events observed for a treatment (if
    arguments \code{event1} and \code{event2} are provided).}
  \item{studies}{Study labels coerced into a factor with its levels
    sorted alphabetically.}
  \item{narms}{Number of arms for each study.}
  \item{TE.nma.fixed, TE.nma.random}{A vector of length \emph{m} of
    consistent treatment effects estimated by network meta-analysis
    (nma) (fixed effect / random effects model).}
  \item{seTE.nma.fixed, seTE.nma.random}{A vector of length \emph{m} of
    effective standard errors estimated by network meta-analysis
    (fixed effect / random effects model).}
  \item{lower.nma.fixed, lower.nma.random}{A vector of length \emph{m} of
    lower confidence interval limits for consistent treatment effects
    estimated by network meta-analysis (fixed effect / random effects
    model).}
  \item{upper.nma.fixed, upper.nma.random}{A vector of length \emph{m} of
    upper confidence interval limits for the consistent treatment
    effects estimated by network meta-analysis (fixed effect /random
    effects model).}
  \item{leverage.fixed}{A vector of length \emph{m} of leverages,
    interpretable as factors by which variances are reduced using
    information from the whole network.}
  \item{w.fixed, w.random}{A vector of length \emph{m} of weights of
    individual studies (fixed effect / random effects model).}
  \item{Q.fixed}{A vector of length \emph{m} of contributions to total
    heterogeneity / inconsistency statistic.}
  \item{TE.fixed, TE.random}{\emph{n}x\emph{n} matrix with estimated overall
    treatment effects (fixed effect / random effects model).}
  \item{seTE.fixed, seTE.random}{\emph{n}x\emph{n} matrix with standard errors
    (fixed effect / random effects model).}
  \item{lower.fixed, upper.fixed, lower.random, upper.random}{\emph{n}x\emph{n}
    matrices with lower and upper confidence interval limits (fixed
    effect / random effects model).}
  \item{zval.fixed, pval.fixed, zval.random, pval.random}{\emph{n}x\emph{n}
    matrices with z-value and p-value for test of overall treatment
    effect (fixed effect / random effects model).}
  \item{seTE.predict}{\emph{n}x\emph{n} matrix with standard errors
    for prediction intervals.}
  \item{lower.predict, upper.predict}{\emph{n}x\emph{n} matrices with
    lower and upper prediction interval limits.}
  \item{prop.direct.fixed, prop.direct.random}{A named vector of the
    direct evidence proportion of each network estimate. (fixed effect
    / random effects model).}
  \item{TE.direct.fixed, TE.direct.random}{\emph{n}x\emph{n} matrix
    with estimated treatment effects from direct evidence (fixed
    effect / random effects model).}
  \item{seTE.direct.fixed, seTE.direct.random}{\emph{n}x\emph{n}
    matrix with estimated standard errors from direct evidence (fixed
    effect / random effects model).}
  \item{lower.direct.fixed, upper.direct.fixed, lower.direct.random,
    upper.direct.random}{\emph{n}x\emph{n} matrices with lower and
    upper confidence interval limits from direct evidence (fixed
    effect / random effects model).}
  \item{zval.direct.fixed, pval.direct.fixed, zval.direct.random,
    pval.direct.random}{\emph{n}x\emph{n} matrices with z-value and
    p-value for test of overall treatment effect from direct evidence
    (fixed effect / random effects model).}
  \item{TE.indirect.fixed, TE.indirect.random}{\emph{n}x\emph{n} matrix
    with estimated treatment effects from indirect evidence (fixed
    effect / random effects model).}
  \item{seTE.indirect.fixed, seTE.indirect.random}{\emph{n}x\emph{n}
    matrix with estimated standard errors from indirect evidence (fixed
    effect / random effects model).}
  \item{lower.indirect.fixed, upper.indirect.fixed, lower.indirect.random,
    upper.indirect.random}{\emph{n}x\emph{n} matrices with lower and
    upper confidence interval limits from indirect evidence (fixed
    effect / random effects model).}
  \item{zval.indirect.fixed, pval.indirect.fixed, zval.indirect.random,
    pval.indirect.random}{\emph{n}x\emph{n} matrices with z-value and
    p-value for test of overall treatment effect from indirect evidence
    (fixed effect / random effects model).}
  \item{Q}{Overall heterogeneity / inconsistency statistic.}
  \item{df.Q}{Degrees of freedom for test of heterogeneity / inconsistency.}
  \item{pval.Q}{P-value for test of heterogeneity / inconsistency.}
  \item{I2}{I-squared.}
  \item{tau}{Square-root of between-study variance.}
  \item{Q.heterogeneity}{Overall heterogeneity statistic.}
  \item{df.Q.heterogeneity}{Degrees of freedom for test of overall
    heterogeneity.}
  \item{pval.Q.heterogeneity}{P-value for test of overall heterogeneity.}
  \item{Q.inconsistency}{Overall inconsistency statistic.}
  \item{df.Q.inconsistency}{Degrees of freedom for test of overall
    inconsistency.}
  \item{pval.Q.inconsistency}{P-value for test of overall inconsistency.}
  \item{Q.decomp}{Data frame with columns 'treat1', 'treat2', 'Q',
    'df' and 'pval.Q', providing heterogeneity statistics for each
    pairwise meta-analysis of direct comparisons.}
  \item{A.matrix}{Adjacency matrix (\emph{n}x\emph{n}).}
  \item{B.matrix}{Edge-vertex incidence matrix (\emph{m}x\emph{n}).}
  \item{L.matrix}{Laplacian matrix (\emph{n}x\emph{n}).}
  \item{Lplus.matrix}{Moore-Penrose pseudoinverse of the Laplacian
    matrix (\emph{n}x\emph{n}).}
  \item{Q.matrix}{Matrix of heterogeneity statistics for pairwise
    meta-analyses, where direct comparisons exist (\emph{n}x\emph{n}).}
  \item{G.matrix}{Matrix with variances and covariances of comparisons
    (\emph{m}x\emph{m}). G is defined as \strong{BL+B^{t}}.}
  \item{H.matrix}{Hat matrix (\emph{m}x\emph{m}), defined as
    \strong{H=GW=BL+B^{t}W}.}
  \item{n.matrix}{\emph{n}x\emph{n} matrix with number of observations
    in direct comparisons (if arguments \code{n1} and \code{n2} are
    provided).}
  \item{events.matrix}{\emph{n}x\emph{n} matrix with number of events
    in direct comparisons (if arguments \code{event1} and
    \code{event2} are provided).}
  \item{P.fixed, P.random}{\emph{n}x\emph{n} matrix with direct
    evidence proportions (fixed effect / random effects model).}
  \item{Cov.fixed}{Variance-covariance matrix (fixed effect model)}
  \item{Cov.random}{Variance-covariance matrix (random effects model)}
  \item{sm, level, level.comb}{As defined above.}
  \item{comb.fixed, comb.random}{As defined above.}
  \item{prediction, level.predict}{As defined above.}
  \item{reference.group, baseline.reference, all.treatments}{As
    defined above.}
  \item{seq, tau.preset, tol.multiarm, details.chkmultiarm}{As defined
    above.}
  \item{sep.trts, nchar.trts}{As defined above.}
  \item{backtransf, title, warn}{As defined above.}
  \item{call}{Function call.}
  \item{version}{Version of R package netmeta used to create object.}
}

\references{
  Jackson D, White IR, Riley RD (2012),
  Quantifying the impact of between-study heterogeneity in multivariate
  meta-analyses.
  \emph{Statistics in Medicine},
  \bold{31}(29), 3805--20

  Rücker G (2012),
  Network meta-analysis, electrical networks and graph theory.
  \emph{Research Synthesis Methods},
  \bold{3}, 312--24

  Rücker G, Schwarzer G (2014),
  Reduce dimension or reduce weights? Comparing two approaches to
  multi-arm studies in network meta-analysis.
  \emph{Statistics in Medicine},
  \bold{33}, 4353--69

  Schwarzer G, Carpenter JR, Rücker G (2015),
  \emph{Meta-Analysis with R (Use-R!)}.
  Springer International Publishing, Switzerland

  Senn S, Gavini F, Magrez D, Scheen A (2013),
  Issues in performing a network meta-analysis.
  \emph{Statistical Methods in Medical Research},
  \bold{22}(2), 169--89
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido
Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{pairwise}}, \code{\link{forest.netmeta}}, \code{\link{netrank}}, \code{\link{metagen}}}

\examples{
data(Senn2013)

#
# Fixed effect model (default)
#
net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD")
net1
net1$Q.decomp

#
# Comparison with reference group
#
print(net1, reference="plac")

#
# Random effects model
#
net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD", comb.random=TRUE)
net2

#
# Change printing order of treatments with placebo last and
# use long treatment names
#
trts <- c("acar", "benf", "metf", "migl", "piog",
          "rosi", "sita", "sulf", "vild", "plac")
net3 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
                data=Senn2013, sm="MD", comb.random=TRUE,
                seq=trts, reference="Placebo")
print(summary(net3), digits=2)
}
