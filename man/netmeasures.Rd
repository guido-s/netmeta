\name{netmeasures}

\alias{netmeasures}

\title{Measures for characterizing a network meta-analysis}

\description{
  This function provides measures for quantifying the direct evidence
  proportion, the mean path length and the minimal parallelism (the
  latter on aggregated and study level) as well as the evidence flow
  per design, see König et al. (2013). These measures support the
  critical evaluation of the network meta-analysis results by
  rendering transparent the process of data pooling.
}

\usage{
netmeasures(x)
}

\arguments{
  \item{x}{An object of class \code{netmeta}.}
}

\details{
  The direct evidence proportion gives the contribution of direct
  effect estimates to a network meta-analysis combined for two-arm and
  multi-arm studies.

  Concerning indirectness, comparisons with a mean path length beyond
  two should be interpreted with particular caution, as more than two
  direct comparisons have to be combined serially on average.
  
  Large indices of parallelism, either on study-level or on aggregated
  level, can be considered as supporting the validity of a network
  meta-analysis if there is only a small amount of heterogeneity.

  The network estimates for two treatments are linear combinations of
  direct effect estimates comparing these or other treatments. The
  linear coefficients can be seen as the generalization of weights
  known from classical meta-analysis. These coefficients are given in
  the projection matrix \eqn{H} of the underlying model. For multi-arm
  studies, the coefficients depend on the choice of the study-specific
  baseline treatment, but the flow of evidence can be made explicit as
  shown in König et al. (2013) and is given in \code{H.tilde}.
}

\value{
A list containing the following components:

  \item{proportion}{A named vector of the direct evidence proportion
    per design.}

  \item{meanpath}{A named vector of the mean path length per
    design.}

  \item{minpar}{A named vector of the minimal parallelism on
    aggregated level per design.}

  \item{minpar.study}{A named vector of the minimal parallelism on
    study level per design.}

  \item{H.tilde }{Design-based hat matrix with information on evidence
    flow per design. The number of rows is equal to the number of
    possible pairwise treatment comparisons and the number of columns
    is equal to the number of designs.}
}

\references{
  König J, Krahn U, Binder H (2013).
  Visualizing the flow of evidence in network meta-analysis and
  characterizing mixed treatment comparisons.
  \emph{Statistics in Medicine}, \bold{32}(30), 5414--29.
}

\author{Ulrike Krahn \email{krahnu@uni-mainz.de}, Jochem König \email{koenigjo@uni-mainz.de}}

\seealso{\link{netmeta}}

\examples{
data(Senn2013)

##
## Generation of an object of class 'netmeta' with
## reference treatment 'plac', i.e. placebo
##
net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
        data=Senn2013, sm="MD", reference="plac")

##
## Calculate measures
##        
m <- netmeasures(net1)

##
## Plot of minimal parallelism versus mean path length
##
plot(m$meanpath, m$minpar, pch="",
     xlab="Mean path length", ylab="Minimal parallelism")
text(m$meanpath, m$minpar, names(m$meanpath), cex=0.8)
}

\keyword{Network meta-analysis}
