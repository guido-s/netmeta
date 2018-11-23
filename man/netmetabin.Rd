\name{netmetabin}
\alias{netmetabin}

\title{Network meta-analysis of binary outcome data}

\description{
  Network meta-analysis of binary outcome data using the
  Mantel-Haenszel method or the non-central hypergeometric
  distribution distribution.
}

\usage{
netmetabin(event1, n1, event2, n2,
           treat1, treat2, studlab,
           data = NULL, subset = NULL,
           sm,
           method = "MH",
           cc.pooled = FALSE,
           incr = gs("incr"),
           allincr = gs("allincr"), addincr = gs("addincr"),
           allstudies = gs("allstudies"),
           level=gs("level"), level.comb=gs("level.comb"),
           comb.fixed = gs("comb.fixed"),
           comb.random = method == "Inverse" &
             (gs("comb.random") | !is.null(tau.preset)),
           ##
           prediction=FALSE, level.predict=gs("level.predict"),
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
           warn = TRUE)
}

\arguments{
  \item{event1}{Number of events (first treatment).}
  \item{n1}{Number of observations (first treatment).}
  \item{event2}{Number of events (second treatment).}
  \item{n2}{Number of observations (second treatment)}
  \item{treat1}{Label/Number for first treatment.}
  \item{treat2}{Label/Number for second treatment.}
  \item{studlab}{An optional - but important! - vector with study labels (see Details).}
  \item{data}{An optional data frame containing the study information.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{sm}{A character string indicating underlying summary measure,
    i.e., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"}.}
  \item{method}{A character string indicating which method is to be
    used for pooling of studies. One of \code{"Inverse"}, \code{"MH"},
    or \code{"NCH"}, can be abbreviated.}
  \item{cc.pooled}{A logical indicating if \code{incr} is considered as
    continuity correction to calculate the network meta-analysis
    estimates.}
  \item{incr}{A numerical value which is added to each cell frequency
    for studies with a zero cell count.}
  \item{allincr}{A logical indicating if \code{incr} is added to each
    cell frequency of all studies if at least one study has a zero cell
    count (only considered if \code{method="Inverse"}). If FALSE
    (default), \code{incr} is added only to each cell frequency of
    studies with a zero cell count.}
  \item{addincr}{A logical indicating if \code{incr} is added to each
    cell frequency of all studies irrespective of zero cell counts (only
    considered if \code{method="Inverse"}).}
  \item{allstudies}{A logical indicating if studies with zero or all
    events in both groups are to be included in an inverse variance
    meta-analysis (applies only if \code{sm} is equal to \code{"RR"} or
    \code{"OR"}).}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
  \item{prediction}{A logical indicating whether a prediction interval
    should be printed.}
  \item{level.predict}{The level used to calculate prediction interval
    for a new study.}
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
  \item{nchar.trts}{A numeric defining the minium number of characters
    used to create unique treatment names (see Details).}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in printouts and forest plots. If
    \code{backtransf=TRUE}, results for \code{sm="OR"} are presented
    as odds ratios rather than log odds ratios, for example.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{keepdata}{A logical indicating whether original data (set)
    should be kept in netmeta object.}
  \item{warn}{A logical indicating whether warnings should be printed
    (e.g., if studies are excluded from meta-analysis due to zero
    standard errors).}
}

\details{
  \bold{... more details ...}

  Comparisons belonging to multi-arm studies are identified by
  identical study labels (argument \code{studlab}). It is therefore
  important to use identical study labels for all comparisons
  belonging to the same multi-arm study, e.g., study label
  "Willms1999" for the three-arm study in Senn et al. (2013).

  Data entry for this function is in \emph{contrast-based} format,
  that is, each line corresponds to a single pairwise comparison
  between two treatments (arguments \code{treat1}, \code{treat2},
  \code{event1}, \code{n1}, \code{event2}, and \code{n2}). If data are
  provided in \emph{arm-based} format, that is, number of events and
  participants are given for each treatment arm separately, function
  \code{\link{pairwise}} can be used to transform the data to
  \emph{contrast-based} format (see help page of function
  \code{\link{pairwise}}).
  
  Note, all pairwise comparisons must be provided for a multi-arm
  study. Consider a multi-arm study of \emph{p} treatments with known
  variances. For this study, treatment effects and standard errors
  must be provided for each of \emph{p}(\emph{p} - 1)/2 possible
  comparisons. For instance, a three-arm study contributes three
  pairwise comparisons, a four-arm study even six pairwise
  comparisons. Function \code{\link{pairwise}} automatically
  calculates all pairwise comparisons for multi-arm studies.
  
  For \code{method = "Inverse"}, both fixed effects and random effects
  models are calculated regardless of values choosen for arguments
  \code{comb.fixed} and \code{comb.random}. Accordingly, the network
  estimates for the random effects model can be extracted from
  component \code{TE.random} of an object of class \code{"netmeta"}
  even if argument \code{comb.random=FALSE}. However, all functions in
  R package \bold{netmeta} will adequately consider the values for
  \code{comb.fixed} and \code{comb.random}. E.g. function
  \code{\link{print.summary.netmeta}} will not print results for the
  random effects model if \code{comb.random=FALSE}. By default
  (\code{comb.fixed=TRUE}), a fixed effects network meta-analysis is
  reported. Use argument \code{comb.random=TRUE} to report results of
  a random effects model.

  For \code{method = "MH"} and \code{method = "NCH"}, only a fixed
  effects model is available.

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
  An object of class \code{netmetabin} and \code{netmeta} with
  corresponding \code{print}, \code{summary}, \code{forest}, and
  \code{netrank} function. The object is a list containing the
  following components:
  \item{studlab, treat1, treat2}{As defined above.}
  \item{n1, n2, event1, event2}{As defined above.}
  \item{TE}{Estimate of treatment effect, i.e. difference between
   first and second treatment (e.g. log odds ratio).}
  \item{seTE}{Standard error of treatment estimate.}   
  \item{k}{Total number of studies.}
  \item{m}{Total number of pairwise comparisons.}
  \item{n}{Total number of treatments.}
  \item{d}{Total number of designs.}
  \item{trts}{Treatments included in network meta-analysis.}
  \item{k.trts}{Number of studies evaluating a treatment.}
  \item{n.trts}{Number of observations receiving a treatment.}
  \item{events.trts}{Number of events observed for a treatment.}
  \item{studies}{Study labels coerced into a factor with its levels
    sorted alphabetically.}
  \item{narms}{Number of arms for each study.}
  \item{designs}{Designs.}
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
  \item{Q}{Overall heterogeneity / inconsistency statistic.}
  \item{df.Q}{Degrees of freedom for test of heterogeneity / inconsistency.}
  \item{pval.Q}{P-value for test of heterogeneity / inconsistency.}
  \item{I2}{I-squared.}
  \item{tau}{Square-root of between-study variance.}
  \item{y}{\bold{??? REMOVE or DESCRIBE ???}}
  \item{V}{\bold{??? REMOVE or DESCRIBE ???}}
  \item{H}{\bold{??? REMOVE or DESCRIBE ???}}
  \item{X}{\bold{??? REMOVE or DESCRIBE ???}}
  \item{A.matrix}{Adjacency matrix (\emph{n}x\emph{n}).}
  \item{H.matrix}{Hat matrix (\emph{m}x\emph{m})}
  \item{n.matrix}{\emph{n}x\emph{n} matrix with number of observations
    in direct comparisons.}
  \item{events.matrix}{\emph{n}x\emph{n} matrix with number of events
    in direct comparisons.}
  \item{data}{Data (contrast-based format).}
  \item{data.wide}{Data in wide (contrast-based) format. \bold{KEEP OR DROP ???}}
  \item{data.long}{Data in long (arm-based) format. \bold{KEEP OR DROP ???}}
  \item{data.design}{List with data in arm-based format (each list
    element corresponds to a single design). \bold{KEEP OR DROP ???}}
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
  Senn S, Gavini F, Magrez D, Scheen A (2013),
  Issues in performing a network meta-analysis.
  \emph{Statistical Methods in Medical Research},
  \bold{22}(2), 169--89
}

\author{Orestis Efthimiou \email{oremiou@gmail.com}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{pairwise}}, \code{\link{netmeta}}}

\examples{
# To be added ...
}
