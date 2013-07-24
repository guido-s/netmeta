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

\usage{netmeta(TE, seTE, treat1, treat2, studlab, data=NULL, subset=NULL,
        sm="", level=0.95, level.comb=0.95,
        comb.fixed=TRUE, comb.random=FALSE, reference.group="",
        all.treatments=NULL, tau.preset=NULL, title="", warn=TRUE)
}

\arguments{
  \item{TE}{Estimate of treatment effect.}
  \item{seTE}{Standard error of treatment estimate.}
  \item{treat1}{Label/Number for first treatment.}
  \item{treat2}{Label/Number for second treatment.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{sm}{A character string indicating underlying summary measure,
    e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"AS"},
    \code{"MD"}, \code{"SMD"}.}
  \item{level}{The level used to calculate confidence intervals for
    individual comparisons.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
  \item{reference.group}{Reference group.}
  \item{all.treatments}{A logical or value \code{"NULL"}. If
        \code{TRUE}, matrices with all treatment effects, and confidence
        limits will be printed.}
  \item{tau.preset}{Prespecified value for between-study variance tau-squared.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{warn}{A logical indicating whether warnings should be printed
    (e.g., if studies are excluded from meta-analysis due to zero
    standard errors).}
}

\details{
  Let \emph{n} be the number of different treatments in a network and
  let \emph{m} be the number of existing comparisons (edges) between
  the treatments. If there are only two-arm studies, \emph{m} is the
  number of studies. Let TE and seTE be the vectors of observed
  effects and their standard errors. Let W be the \emph{m}x\emph{m}
  diagonal matrix that contains the inverse variance 1/seTE^2.
  
  The given comparisons define the network structure. Therefrom an
  \emph{m}x\emph{n} design matrix B is formed; for more precise information, see
  Rücker (2012). Moreover, the \emph{n}x\emph{n} Laplacian matrix L and its
  Moore-Penrose pseudoinverse L+ are calculated (both matrices play an
  important role in graph theory and electrical network theory).
  Using these matrices, the variances based on both direct and
  indirect comparisons can be estimated. Moreover, the hat matrix H
  can be estimated by \strong{H = BL+B^{t}W = B(B^{t} W B)^{+}B^{t}W} and
  finally consistent treatment effects can be estimated by applying
  the hat matrix to the observed (potentially inconsistent) effects. H
  is a projection matrix which maps the observed effects onto the
  consistent (n-1)-dimensional subspace. This is the Aitken
  estimator (Senn et al., 2013). As in pairwise meta-analysis, the Q
  statistic measures the deviation from consistency. Q can be
  separated into parts for each pairwise meta-analysis and a part for
  remaining inconsistency between comparisons.

  Often multi-armed studies are included in a network meta-analysis.
  In multi-armed studies, the treatment effects on different
  comparisons are not independent, but correlated. This is accounted
  for by reweighting all edges of each multi-armed study. The method
  is described in Rücker (2012).

  A simple random effects model assuming that a constant heterogeneity
  variance is added to each edge of the network can be defined via a
  generalised methods of moments estimate of the between-studies
  variance tau^2 (Jackson et al., 2012). This is added to the observed
  sampling variance seTE^2 of each edge in the network (after
  appropriate adjustment for multi-armed studies). Then, as in
  standard pairwise meta-analysis, the procedure is repeated with the
  resulting enlarged standard errors.
}

\value{
  An object of class \code{c("netmeta")} with corresponding
  \code{print}, \code{summary}, \code{forest} function. The object is
  a list containing the following components:
  \item{TE, seTE, studlab, treat1, treat2, sm, level, level.comb}{As
    defined above.}
  \item{comb.fixed, comb.random, tau.preset, title, warn}{As defined
    above.}
  \item{seTE.adj}{Standard error of treatment estimate, adjusted for
    multi-armed studies.}
  \item{reference.group}{The name of the reference group, if
    specified, otherwise \code{c("")}.}
  \item{all.treatments}{A logical or value \code{"NULL"}. If
        \code{TRUE}, matrices with all treatment effects, and
        confidence limits will be printed.}
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
  \item{Q.fixed}{A vector of length \emph{m} of contributions to total
    heterogeneity / inconsistency statistic.}
  \item{k}{Total number of studies.}
  \item{m}{Total number of pairwise comparisons.}
  \item{Q}{Overall heterogeneity / inconsistency statistic.}
  \item{df}{Degrees of freedom for test of heterogeneity / inconsistency.}
  \item{pval.Q}{P-value for test of heterogeneity / inconsistency.}
  \item{I2}{I-squared.}
  \item{tau}{Square-root of between-study variance.}
  \item{Q.heterogeneity}{Overall heterogeneity statistic.}
  \item{Q.inconsistency}{Overall inconsistency statistic.}
  \item{A.matrix}{Adjacency matrix (\emph{n}x\emph{n}).}
  \item{L.matrix}{Laplacian matrix (\emph{n}x\emph{n}).}
  \item{Lplus.matrix}{Moore-Penrose pseudoinverse of the Laplacian
    matrix (\emph{n}x\emph{n}).}
  \item{Q.matrix}{Matrix of heterogeneity statistics for pairwise
    meta-analyses, where direct comparisons exist (\emph{n}x\emph{n}).}
  \item{G.matrix}{Matrix with variances and covariances of edges
    (\emph{m}x\emph{m}). G is defined as \strong{BL+B^{t}}.}
  \item{H.matrix}{Hat matrix (\emph{m}x\emph{m}), defined as \strong{H=GW=BL+B^{t}W}.}
  \item{Q.decomp}{Data frame with columns 'treat1', 'treat2', 'Q',
    'df' and 'pval.Q', providing heterogeneity statistics for each
    pairwise meta-analysis of direct comparisons.}
  \item{call}{Function call.}
  \item{version}{Version of R package netmeta used to create object.}
}

\references{
  Jackson D, White IR and Riley RD (2012),
  Quantifying the impact of between-study heterogeneity in multivariate
  meta-analyses.
  \emph{Statistics in Medicine},
  \bold{31}(29). 3805--3820.

  Rücker G (2012),
  Network meta-analysis, electrical networks and graph theory.
  \emph{Research Synthesis Methods},
  \bold{3}, 312--324.

  Senn S, Gavini F, Magrez D, and Scheen A (2013).
  Issues in performing a network meta-analysis.
  \emph{Statistical Methods in Medical Research},
  \bold{22}(2), 169--189.
  First published online 2012 Jan 3.
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido
Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{forest.netmeta}}, \code{\link{metagen}}}

\examples{

data(Senn2013)

#### Fixed effect model (default)

net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD")
net1
net1$Q.decomp

#### Comparison with reference group

netmeta(TE, seTE, treat1, treat2, studlab,
        data=Senn2013, sm="MD", reference="plac")

#### Random effects model

net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD", comb.random=TRUE)
net2

}
