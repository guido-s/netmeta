\name{smokingcessation}
\docType{data}
\alias{smokingcessation}
\title{Network meta-analysis of interventions for smoking cessation}

\description{Network meta-analysis comparing the effects of a number
  of interventions for smoking cessation.

  These data are used as an example in Dias et al. (2013), page 651.
}

\usage{data(smokingcessation)}

\format{A data frame with the following columns:
  \describe{
  \item{event1}{Number of events in arm 1}
  \item{n1}{Number of observations in arm 1}
  \item{event2}{Number of events in arm 2}
  \item{n2}{Number of observations in arm 2}
  \item{event3}{Number of events in arm 3}
  \item{n3}{Number of observations in arm 3}
  \item{treat1}{Treatment 1}
  \item{treat2}{Treatment 2}
  \item{treat3}{Treatment 3}
  }
}

\source{
  Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G and Ades AE (2013).
  Evidence Synthesis for Decision Making 4: Inconsistency in networks
  of evidence based on randomized controlled trials.
  \emph{Medical Decision Making}
  \bold{33}, 641--656.
}

\seealso{\code{\link{pairwise}}, \code{\link{metabin}},
\code{\link{netmeta}}, \code{\link{netgraph}}}

\examples{
data(smokingcessation)

# Transform data from arm-based format to contrast-based format
# Argument 'sm' has to be used for odds ratio as summary measure; by
# default the risk ratio is used in the metabin function called
# internally.
p1 <- pairwise(list(treat1, treat2, treat3),
               event=list(event1, event2, event3),
               n=list(n1, n2, n3),
               data=smokingcessation,
               sm="OR")
p1

# Conduct network meta-analysis
net1 <- netmeta(TE, seTE, treat1, treat2, studlab, data=p1)
net1

# Draw network graph
netgraph(net1, points=TRUE, cex.points=3, cex=1.25)
tname <- c("No intervention","Self-help","Individual counselling","Group
counselling")
netgraph(net1, points=TRUE, cex.points=3, cex=1.25, labels=tname)
}

\keyword{datasets}
