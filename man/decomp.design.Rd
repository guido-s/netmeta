\name{decomp.design}

\alias{decomp.design}

\title{Design-based decomposition of Cochran's Q in network meta-analysis}

\description{
  This function performs a design-based decomposition of Cochran's Q
  for assessing the homogeneity in the whole network, the homogeneity
  within designs, and the homogeneity/consistency between designs. It
  allows also an assessment of the consistency assumption after
  detaching the effect of single designs.
}

\usage{
decomp.design(x)
}

\arguments{
   \item{x}{An object of class \code{netmeta}.}
}

\details{
  In the context of network meta-analysis and the assessment of the
  homogeneity and consistency assumption, a generalized Cochran's Q
  statistic for multivariate meta-analysis can be used as shown in
  Krahn et al. (2013). This Q statistic can be decomposed in a sum of
  within-design Q statistics and one between-designs Q statistic that
  incorporates the concept of design inconsistency, see Higgins et
  al. (2012).

  Since an inconsistent treatment effect of one design can
  simultaneously inflate several residuals, Krahn et al. (2013)
  suggest for locating the inconsistency in a network to fit a set of
  extended models allowing for example for a deviating effect of each
  study design in turn. The recalculated between-designs Q statistics
  are given in list component \code{Q.detach}. The change of the
  inconsistency contribution of single designs can be investigated in
  more detail by a net heat plot (see function
  \link{netheat}). Designs where only one treatment is involved in
  other designs of the network or where the removal of corresponding
  studies would lead to a splitting of the network do not contribute
  to the inconsistency assessment. These designs are not included in
  \code{Q.detach}.
}

\value{
  A list containing the following components:
  \item{Q.decomp}{Data frame with Q statistics (variable \code{Q}) to
    assess the homogeneity/consistency in the whole network, within
    designs, and between designs. Corresponding degrees of freedom
    (\code{df}) and p-values (\code{pval}) are also given.}
  \item{Q.design}{Data frame with design-specific decomposition of the
    within-designs Q statistic (\code{Q}), corresponding degrees of
    freedom (\code{df}) and p-values (\code{pval}) as well as
    information on the designs (\code{design}) are given.}
  \item{Q.detach}{Data frame with between-designs Q statistics
    (\code{Q}) after detaching of single designs, corresponding
    degrees of freedom (\code{df}) and p-values (\code{pval}) as well
    as information on the detached design (\code{design}) are given.}
  \item{residuals}{Matrix with residuals.}
  \item{Q.inc.design}{A named vector with contribution of single
    designs to the between design Q statistic.}
  \item{call}{Function call.}
  \item{version}{Version of R package netmeta used to create object.}
}

\references{
  Higgins JPT, Jackson D, Barrett JK, Lu G, Ades AE, White IR (2012),
  Consistency and inconsistency in network meta-analysis:
  concepts and models for multi-arm studies.
  \emph{Research Synthesis Methods}, \bold{3}(2), 98--110.

  Krahn U, Binder H, König J (2013),
  A graphical tool for locating inconsistency in network
  meta-analyses.
  \emph{BMC Medical Research Methodology}, \bold{13}, 35.
}

\author{Ulrike Krahn \email{krahnu@uni-mainz.de}}

\seealso{\link{netmeta}, \link{netheat}}

\examples{
##
## Note: Contrary to these data, in Krahn et al. (2013) slightly
## different standard errors related to individual studies are used
## calculated by the function 'metacont' of R package 'meta'.
##
data(Senn2013)

##
## Generation of an object of class 'netmeta' with
## reference treatment 'plac', i.e. placebo
##
net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD", reference="plac")

##
## Decomposition of Cochran's Q 
##
decomp.design(net1) 
}

\keyword{Network meta-analysis}
\keyword{Cochran's Q}
\keyword{Inconsistency}
