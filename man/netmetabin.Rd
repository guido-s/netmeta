\name{netmetabin}
\alias{netmetabin}

\title{Network meta-analysis of binary outcome data}

\description{
  Provides three models for the network meta-analysis of binary data
  (Mantel-Haenszel method, based on the non-central hypergeometric
  distribution, and the inverse variance method).
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
  \item{cc.pooled}{A logical indicating whether \code{incr} should be
    used as a continuity correction, when calculating the network
    meta-analysis estimates.}
  \item{incr}{A numerical value which is added to each cell count, i.e.,
    to the numbers of events and non-events, of all treatment arms in
    studies with zero events or non-events in any of the treatment arms
    ("continuity correction").}
  \item{allincr}{A logical indicating whether \code{incr} should be
    added to each cell count of all studies if a continuity correction
    was used for at least one study (only considered if
    \code{method="Inverse"}). If FALSE (default), \code{incr} is used as
    continuity correction only for studies with zero events or zero
    non-events in any of the treatment arms.}
  \item{addincr}{A logical indicating whether \code{incr} should be
    added to each cell count of all studies, irrespective of zero cell
    counts (only considered if \code{method="Inverse"}).}
  \item{allstudies}{A logical indicating whether studies with zero
    events or non-events in all treatment arms should be included in an
    inverse variance meta-analysis (applies only if
    \code{method="Inverse"} and \code{sm} is equal to either \code{"RR"}
    or \code{"OR"}).}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.fixed}{A logical indicating whether a fixed effects (common
    effects) network meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    network meta-analysis should be conducted.}
  \item{prediction}{A logical indicating whether a prediction interval
    should be printed (only considered if \code{method="Inverse"}).}
  \item{level.predict}{The level used to calculate prediction interval
    for a new study (only considered if \code{method="Inverse"}).}
  \item{reference.group}{Reference treatment.}
  \item{baseline.reference}{A logical indicating whether results
    should be expressed as comparisons of other treatments versus the
    reference treatment (default) or vice versa. This argument is only
    considered if \code{reference.group} has been specified.}
  \item{all.treatments}{A logical or \code{"NULL"}. If \code{TRUE},
    matrices with all treatment effects, and confidence limits will be
    printed.}
  \item{seq}{A character or numerical vector specifying the sequence
    of treatments in printouts.}
  \item{tau.preset}{An optional value for manually setting the
    square-root of the between-study variance \eqn{\tau^2} (only
    considered if \code{method="Inverse"}).}
  \item{tol.multiarm}{A numeric for the tolerance for consistency of
    treatment estimates and corresponding variances in multi-arm studies
    which are consistent by design (only considered if
    \code{method="Inverse"}).}
  \item{details.chkmultiarm}{A logical indicating whether treatment
    estimates and / or variances of multi-arm studies with inconsistent
    results or negative multi-arm variances should be printed (only
    considered if \code{method="Inverse"}).}
  \item{sep.trts}{A character used in comparison names as separator
    between treatment labels.}
  \item{nchar.trts}{A numeric defining the minimum number of characters
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
  
  This function implements three models for the network meta-analysis of
  binary data:
  \itemize{
   \item The Mantel-Haenszel network meta-analysis model, as described
     in Efthimiou et al. (2018) (\code{method="MH"});
   \item a network meta-analysis model using the non-central
     hypergeometric distribution with the Breslow approximation, as
     described in Stijnen et al. (2010) (\code{method="NCH"}),
   \item the inverse variance method for network meta-analysis
     (\code{method="Inverse"}), also provided by \code{\link{netmeta}}.
  }

  Comparisons belonging to multi-arm studies are identified by identical
  study labels (argument \code{studlab}). It is therefore important to
  use identical study labels for all comparisons belonging to the same
  multi-arm study.

  Data entry for this function is in \emph{contrast-based} format, that
  is, each line of the data corresponds to a single pairwise comparison
  between two treatments (arguments \code{treat1}, \code{treat2},
  \code{event1}, \code{n1}, \code{event2}, and \code{n2}). If data are
  provided in \emph{arm-based} format, that is, number of events and
  participants are given for each treatment arm separately, function
  \code{\link{pairwise}} can be used to transform the data to
  \emph{contrast-based} format (see help page of function
  \code{\link{pairwise}}).
  
  Note, all pairwise comparisons must be provided for a multi-arm
  study. Consider a multi-arm study of \emph{p} treatments with known
  variances. For this study, the number of events and observations must
  be provided for each treatment, for each of \emph{p}(\emph{p} - 1)/2
  possible comparisons in separate lines in the data. For instance, a
  three-arm study contributes three pairwise comparisons, a four-arm
  study even six pairwise comparisons. Function \code{\link{pairwise}}
  automatically calculates all pairwise comparisons for multi-arm
  studies.
  
  For \code{method = "Inverse"}, both fixed effects and random effects
  models are calculated regardless of values choosen for arguments
  \code{comb.fixed} and \code{comb.random}. Accordingly, the network
  estimates for the random effects model can be extracted from component
  \code{TE.random} of an object of class \code{"netmeta"} even if
  argument \code{comb.random=FALSE}. However, all functions in R package
  \bold{netmeta} will adequately consider the values for
  \code{comb.fixed} and \code{comb.random}. E.g. function
  \code{\link{print.summary.netmeta}} will not print results for the
  random effects model if \code{comb.random=FALSE}.

  For \code{method = "MH"} and \code{method = "NCH"}, only a fixed
  effects model is available.

  By default, treatment names are not abbreviated in printouts. However,
  in order to get more concise printouts, argument \code{nchar.trts} can
  be used to define the minimum number of characters for abbreviated
  treatment names (see \code{\link{abbreviate}}, argument
  \code{minlength}). R function \code{\link{treats}} is utilised
  internally to create abbreviated treatment names.

  Names of treatment comparisons are created by concatenating treatment
  labels of pairwise comparisons using \code{sep.trts} as separator (see
  \code{\link{paste}}). These comparison names are used in the
  covariance matrices \code{Cov.fixed} and \code{Cov.random} and in some
  R functions, e.g, \code{\link{decomp.design}}. By default, a colon is
  used as the separator. If any treatment label contains a colon the
  following characters are used as separator (in consecutive order):
  \code{"-"}, \code{"_"}, \code{"/"}, \code{"+"}, \code{"."},
  \code{"|"}, and \code{"*"}.  If all of these characters are used in
  treatment labels, a corresponding error message is printed asking the
  user to specify a different separator.
}

\value{
  An object of class \code{netmetabin} and \code{netmeta} with
  corresponding \code{print}, \code{summary}, \code{forest}, and
  \code{netrank} functions. The object is a list containing the
  following components:
  \item{studlab, treat1, treat2}{As defined above.}
  \item{n1, n2, event1, event2}{As defined above.}
  \item{TE}{Estimate of treatment effect, i.e. difference between
   first and second treatment (e.g. log odds ratio).}
  \item{seTE}{Standard error of treatment estimate.}   
  \item{k}{Total number of studies.}
  \item{m}{Total number of pairwise comparisons.}
  \item{n}{Total number of treatments.}
  \item{d}{Total number of designs (corresponding to the unique set of
    treatments compared within studies).}
  \item{trts}{Treatments included in network meta-analysis.}
  \item{k.trts}{Number of studies evaluating a treatment.}
  \item{n.trts}{Number of observations receiving a treatment.}
  \item{events.trts}{Number of events observed for a treatment.}
  \item{studies}{Study labels coerced into a factor with its levels
    sorted alphabetically.}
  \item{narms}{Number of arms for each study.}  
  \item{designs}{Unique list of designs present in the network. A design
    corresponds to the set of treatments compared within a study.}
  \item{TE.fixed, seTE.fixed}{\emph{n}x\emph{n} matrix with estimated
    overall treatment effects and standard errors for fixed effects
    model.}  
  \item{lower.fixed, upper.fixed}{\emph{n}x\emph{n} matrices with lower
    and upper confidence interval limits for fixed effects model.}
  \item{zval.fixed, pval.fixed}{\emph{n}x\emph{n}
    matrices with z-value and p-value for test of overall treatment
    effect under fixed effects model.}  
  \item{TE.random, seTE.random}{\emph{n}x\emph{n} matrix with estimated
    overall treatment effects and standard errors for random effects
    model (only available if \code{method="Inverse"}).}  
  \item{lower.random, upper.random}{\emph{n}x\emph{n} matrices with
    lower and upper confidence interval limits for random effects model
    (only available if \code{method="Inverse"}).}  
  \item{zval.random, pval.random}{\emph{n}x\emph{n} matrices with
    z-value and p-value for test of overall treatment effect under
    random effects model (only available if \code{method="Inverse"}).}  
  \item{TE.direct.fixed, seTE.direct.fixed}{\emph{n}x\emph{n} matrix
    with estimated treatment effects and standard errors from direct
    evidence under fixed effects model.}
  \item{lower.direct.fixed, upper.direct.fixed}{\emph{n}x\emph{n}
    matrices with lower and upper confidence interval limits from direct
    evidence under fixed effects model.}
  \item{zval.direct.fixed, pval.direct.fixed}{\emph{n}x\emph{n} matrices
    with z-value and p-value for test of overall treatment effect from
    direct evidence under fixed effects model.}
  \item{TE.direct.random, seTE.direct.random}{\emph{n}x\emph{n} matrix
    with estimated treatment effects and standard errors from direct
    evidence under random effects model (only available if
    \code{method="Inverse"}).}
  \item{lower.direct.random, upper.direct.random}{\emph{n}x\emph{n}
    matrices with lower and upper confidence interval limits from direct
    evidence under random effects model (only available if
    \code{method="Inverse"}).}
  \item{zval.direct.random, pval.direct.random}{\emph{n}x\emph{n}
    matrices with z-value and p-value for test of overall treatment
    effect from direct evidence under random effects model (only
    available if \code{method="Inverse"}).}
  \item{Q}{Overall heterogeneity / inconsistency statistic.}
  \item{df.Q}{Degrees of freedom for test of heterogeneity / inconsistency.}
  \item{pval.Q}{P-value for test of heterogeneity / inconsistency.}
  \item{I2}{I-squared (only available if \code{method="Inverse"}).}
  \item{tau}{Square-root of between-study variance (only available if
    \code{method="Inverse"}).}
  \item{A.matrix}{Adjacency matrix (\emph{n}x\emph{n}).}
  \item{H.matrix}{Hat matrix (\emph{m}x\emph{m})}
  \item{n.matrix}{\emph{n}x\emph{n} matrix with number of observations
    in direct comparisons.}
  \item{events.matrix}{\emph{n}x\emph{n} matrix with number of events
    in direct comparisons.}
  \item{sm, method, level, level.comb}{As defined above.}
  \item{incr, allincr, addincr, allstudies, cc.pooled}{As defined above.}
  \item{comb.fixed, comb.random}{As defined above.}
  \item{prediction, level.predict}{As defined above.}
  \item{reference.group, baseline.reference, all.treatments}{As
    defined above.}
  \item{seq, tau.preset, tol.multiarm, details.chkmultiarm}{As defined
    above.}
  \item{sep.trts, nchar.trts}{As defined above.}
  \item{backtransf, title, warn}{As defined above.}
  \item{data}{Data set (in contrast-based format).}
  \item{data.design}{List with data in arm-based format (each list
    element corresponds to a single design).}
  \item{call}{Function call.}
  \item{version}{Version of R package netmeta used to create object.}
}

\references{
  Efthimiou O, Rücker G, Schwarzer G, Higgins J, Egger M, Salanti G
  (2018):
  A Mantel-Haenszel model for network meta-analysis of rare events.
  \emph{Manuscript submitted for publication}.  
  
  Senn S, Gavini F, Magrez D, Scheen A (2013):
  Issues in performing a network meta-analysis.
  \emph{Statistical Methods in Medical Research},
  \bold{22}, 169--89
  
  Stijnen T, Hamza TH, Ozdemir P (2010):
  Random effects meta-analysis of event outcome in the framework of
  the generalized linear mixed model with applications in sparse
  data.
  \emph{Statistics in Medicine}, \bold{29}, 3046--67
}

\author{Orestis Efthimiou \email{oremiou@gmail.com}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{pairwise}}, \code{\link{netmeta}}}

\examples{
data(Dong2013)

# Only consider first ten studies (to reduce runtime of example)
#
first10 <- subset(Dong2013, id <= 10)

# Transform data from long arm-based format to contrast-based format.
# Argument 'sm' has to be used for odds ratio as summary measure; by
# default the risk ratio is used in the metabin function called
# internally.
p1 <- pairwise(treatment, death, randomized, studlab = id,
               data = first10, sm = "OR")

# Conduct Mantel-Haenszel network meta-analysis
# (without continuity correction)
nb1 <- netmetabin(p1, ref = "plac")
nb1

# Obtain the league table
netleague(nb1)

\dontrun{
# Conduct Mantel-Haenszel network meta-analysis for the whole dataset
#
p2 <- pairwise(treatment, death, randomized, studlab = id,
               data = Dong2013, sm = "OR")
netmetabin(p2, ref = "plac")
  
# Conduct network meta-analysis using the non-central hypergeometric
# model (without continuity correction)
netmetabin(p2, ref = "plac", method = "NCH")

# Conduct Mantel-Haenszel network meta-analysis
# (with continuity correction of 0.5; include all studies)
netmetabin(p2, ref = "plac", cc.pooled = TRUE)


data(Gurusamy2011)

p3 <- pairwise(treatment, death, n, studlab = study,
               data = Gurusamy2011, sm = "OR")

# Conduct Mantel-Haenszel network meta-analysis
# (without continuity correction)
netmetabin(p3, ref = "cont")
}
}
