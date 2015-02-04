\name{parkinson}
\docType{data}
\alias{parkinson}
\title{Network meta-analysis of treatments for Parkinson's disease}

\description{Network meta-analysis comparing the effects of a number
  of treatments for Parkinson's disease.

  The data are the mean lost work-time reduction in patients given
  dopamine agonists as adjunct therapy in Parkinson’s disease.  The
  data are given as sample size, mean and standard deviation in each
  trial arm. Treatments are placebo, coded 1, and four active drugs
  coded 2 to 5. These data are used as an example in the supplemental
  material of Dias et al. (2013).
}

\usage{data(parkinson)}

\format{A data frame with the following columns:
  \describe{
  \item{Study}{Study label}
  \item{Treatment1}{Treatment 1}
  \item{y1}{Treatment effect arm 1}
  \item{sd1}{Standard deviation arm 1}
  \item{n1}{Sample size arm 1}
  \item{Treatment2}{Treatment 2}
  \item{y2}{Treatment effect arm 2}
  \item{sd2}{Standard deviation arm 2}
  \item{n2}{Sample size arm 2}
  \item{Treatment3}{Treatment 3}
  \item{y3}{Treatment effect arm 3}
  \item{sd3}{Standard deviation arm 3}
  \item{n3}{Sample size arm 3}
  }
}

\source{
  Dias S, Sutton AJ, Ades AE and Welton NJ (2013).
  Evidence synthesis for decision making 2: A generalized linear
  modeling framework for pairwise and network meta-analysis of
  randomized controlled trials.
  \emph{Medical Decision Making} \bold{33}, 607--617.
}

\seealso{\code{\link{pairwise}}, \code{\link{metacont}},
\code{\link{netmeta}}, \code{\link{netgraph}}}

\examples{
data(parkinson)

# Transform data from arm-based format to contrast-based format
p1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
               n=list(n1, n2, n3),
               mean=list(y1, y2, y3),
               sd=list(sd1, sd2, sd3),
               data=parkinson, studlab=Study)
p1

# Conduct network meta-analysis
net1 <- netmeta(TE, seTE, treat1, treat2, studlab, data=p1)
net1

# Draw network graphs
netgraph(net1, points=TRUE, cex.points=3, cex=1.5,
         thickness="se.fixed")
netgraph(net1, points=TRUE, cex.points=3, cex = 1.5,
         plastic=TRUE, thickness="se.fixed",
         iterate=TRUE)
netgraph(net1, points=TRUE, cex.points=3, cex = 1.5,
         plastic=TRUE, thickness="se.fixed",
         iterate=TRUE, start="eigen")
}

\keyword{datasets}
