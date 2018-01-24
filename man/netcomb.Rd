\name{netcomb}

\alias{netcomb}
\alias{discomb}

\title{Additive network meta-analysis for combinations of treatments}

\description{
  Some treatments in a network meta-analysis may be combinations of
  other treatments or have common components. The influence of
  individual components can be evaluated in an additive network
  meta-analysis model assuming that the effect of treatment
  combinations is the sum of the effects of its components. This
  function implements this additive model in a frequentist way.
}

\usage{
netcomb(x,
        inactive = NULL,
        sep.components = "+",
        C.matrix, ...)

discomb(TE, seTE,
        treat1, treat2,
        studlab, data = NULL, subset = NULL,
        ##
        sep.components = "+", 
        C.matrix,
        inactive = NULL,
        ##
        sm,
        level.comb = gs("level.comb"),
        comb.fixed = gs("comb.fixed"),
        comb.random = FALSE | !is.null(tau.preset),
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
        warn = TRUE,
        ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta}.}
  \item{inactive}{A character string defining the inactive treatment
    (see Details).}
  \item{sep.components}{A single character to define separator between
    treatment components.}
  \item{C.matrix}{C matrix (see Details).}
  \item{TE}{Estimate of treatment effect, i.e. difference between
   first and second treatment (e.g. log odds ratio, mean difference,
   or log hazard ratio).}
  \item{seTE}{Standard error of treatment estimate.}
  \item{treat1}{Label/Number for first treatment.}
  \item{treat2}{Label/Number for second treatment.}
  \item{studlab}{An optional - but important! - vector with study
    labels (see \code{\link{netmeta}}).}
  \item{data}{An optional data frame containing the study information.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{sm}{A character string indicating underlying summary measure,
    e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
    \code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
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
  \item{title}{Title of meta-analysis / systematic review.}
  \item{warn}{A logical indicating whether warnings should be printed
    (e.g., if studies are excluded from meta-analysis due to zero
    standard errors).}
  \item{\dots}{Additional arguments.}
}

\details{
Treatments in network meta-analysis (NMA) can be complex
interventions. Some treatments may be combinations of others or have
common components. The standard analysis provided by
\code{\link{netmeta}} is a NMA where all existing (single or combined)
treatments are considered as different nodes in the
network. Exploiting the fact that some treatments are combinations of
common components, an additive NMA model can be used to evaluate the
influence of individual components. This model assumes that the effect
of a treatment combination is the sum of the effects of its components
which implies that common components cancel out in comparisons.

The additive network meta-analysis model has been implemented using
Bayesian methods (Mills et al., 2012; Welton et al., 2013). This
function implements the additive model in a frequentist way (Rücker &
Schwarzer, 2017).

The underlying multivariate model is given by

\deqn{\bold{\delta} = \bold{B} \bold{\theta}, \bold{\theta} = \bold{C} \bold{\beta}}

with
\itemize{
\item[\eqn{\bold{\delta}}] vector of true treatment effects
(differences) from individual studies,
\item[\eqn{\bold{B}}] is a design matrix describing the structure of
the network,
\item[\eqn{\bold{\theta}}] parameter vector that represents the
existing combined treatments,
\item[\eqn{\bold{C}}] matrix describing how the treatments are
composed,
\item[\eqn{\bold{\beta}}] is a parameter vector representing the
treatment components.
}
All parameters are estimated using weighted least squares regression.

Argument \code{inactive} can be used to specify a single component
that does not have any therapeutic value. Accordingly, it is assumed
that the treatment effect of the combination of this component with an
additional treatment component is equal to the treatment effect of the
additional component alone.

Argument \code{sep.components} can be used to specify the separator
between individual components. By default, the matrix \strong{C} is
calculated internally from treatment names. However, it is possible to
specify a different matrix using argument \code{C.matrix}.
}                                                              

\value{
A list containing the following components:

  \item{x, inactive, sep.components, C.matrix}{As defined above.}  
  \item{k}{Total number of studies.}
  \item{m}{Total number of pairwise comparisons.}
  \item{n}{Total number of treatments.}
  \item{c}{Total number of components.}
  \item{comparisons.fixed, comparisons.random}{Lists with components
    studlab, treat1, treat2, TE, seTE, lower, upper, z, p level, df
    (fixed and random effects model).}
  \item{components.fixed, components.random}{Lists with components TE,
    seTE, lower, upper, z, p level, df (fixed and random effects
    model).}
  \item{combinations.fixed, combinations.random}{Lists with components
    TE, seTE, lower, upper, z, p level, df (fixed and random effects
    model).}
  \item{sm}{Summary measure.}
  \item{level.comb}{Level for confidence intervals.}
  \item{comb.fixed}{Fixed effect model}
  \item{comb.random}{Random effects model}
  \item{Q}{Overall heterogeneity / inconsistency statistic.}
  \item{df.Q}{Degrees of freedom for test of heterogeneity /
    inconsistency.}
  \item{pval.Q}{P-value for test of heterogeneity / inconsistency.}   
  \item{Q.comp.fixed, Q.comp.random}{Test statistic for difference in
    goodness of fit between standard model and additive model (fixed and
    random effects model).}
  \item{df.Q.comp.fixed, df.Q.comp.random}{Degrees of freedom for
    difference in goodness of fit between standard model and additive
    model (fixed and random effects model).}
  \item{pval.Q.comp.fixed, pval.Q.comp.random}{P-value for difference
    in goodness of fit between standard model and additive model
    (fixed and random effects model).}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in printouts and forest plots.}
  \item{nchar.trts}{A numeric defining the minium number of characters
    used to create unique treatment and component names.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{call}{Function call.}
  \item{version}{Version of R package netmeta used to create object.}
}

\references{
  König J, Krahn U, Binder H (2013).
  Visualizing the flow of evidence in network meta-analysis and
  characterizing mixed treatment comparisons.
  \emph{Statistics in Medicine}, \bold{32}(30), 5414--29.

  Mills EJ, Thorlund K, Ioannidis JP (2012).
  Calculating additive treatment effects from multiple randomized
  trials provides useful estimates of combination therapies.
  \emph{Journal of Clinical Epidemiology}, \bold{65}(12), 1282--8.
  Doi: 10.1016/j.jclinepi.2012.07.012.

  Rücker G, Schwarzer G (2017).
  Network meta-analysis of combinations of treatments.
  Manuscript.

  Welton NJ, Caldwell DM, Adamopoulos E, Vedhara K (2009).
  Mixed treatment comparison meta-analysis of complex interventions:
  psychological interventions in coronary heart disease.
  \emph{American Journal of Epidemiology}, \bold{169}(9): 1158--65.
  Doi: 10.1093/aje/kwp014.
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido
Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\link{netmeta}, \link{print.netcomb}}

\examples{
data(Linde2016)
#
# Specify, order of treatments
#
trts <- c("TCA", "SSRI", "SNRI", "NRI", "Low-dose SARI", "NaSSa",
          "rMAO-A", "Ind drug", "Hypericum", "Face-to-face CBT",
          "Face-to-face PST", "Face-to-face interpsy", "Face-to-face psychodyn",
          "Other face-to-face", "Remote CBT", "Self-help CBT", "No contact CBT",
          "Face-to-face CBT + SSRI", "Face-to-face interpsy + SSRI",
          "Face-to-face PST + SSRI", "UC", "Placebo")
#
# Note, three treatments are actually combinations of 'SSRI' with
# other components:
#
# "Face-to-face CBT + SSRI",
# "Face-to-face interpsy + SSRI",
# "Face-to-face PST + SSRI"

# Standard Network meta-analysis
#
net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
                data = Linde2016, ref = "placebo",
		seq = trts,
                sm = "OR", comb.fixed = FALSE, comb.random = TRUE)
#
summary(net1)

# Additive model for treatment components (with placebo as inactive
# treatment)
#
nc1 <- netcomb(net1, inactive = "placebo")
#
summary(nc1)
}                                                                 

\keyword{Network meta-analysis}
\keyword{Additive network meta-analysis}
