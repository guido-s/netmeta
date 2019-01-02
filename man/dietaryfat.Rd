\name{dietaryfat}
\docType{data}
\alias{dietaryfat}
\title{Network meta-analysis of dietary fat}

\description{Network meta-analysis comparing the effects of two diets to
  control on mortality.

  The data are rates, given as the number of deaths and person-years. 
  These data are used as an example in the supplemental material
  of Dias et al. (2013).
}

\usage{data(dietaryfat)}

\format{A data frame with the following columns:
  \describe{
  \item{treat1}{Treatment 1}
  \item{treat2}{Treatment 2}
  \item{treat3}{Treatment 3}
  \item{years1}{Person years arm 1}
  \item{years2}{Person years arm 2}
  \item{years3}{Person years arm 3}
  \item{d1}{events (deaths) arm 1}
  \item{d2}{events (deaths) arm 2}
  \item{d3}{events (deaths) arm 3}
  \item{ID}{Study ID}
  }
}

\source{
  Dias S, Sutton AJ, Ades AE and Welton NJ (2013):
  Evidence synthesis for decision making 2: A generalized linear
  modeling framework for pairwise and network meta-analysis of
  randomized controlled trials.
  \emph{Medical Decision Making} \bold{33}, 607--17
}

\seealso{\code{\link{pairwise}}, \code{\link{metainc}},
\code{\link{netmeta}}, \code{\link{netgraph}}}

\examples{
data(dietaryfat)

# Transform data from arm-based format to contrast-based format
# Using incidence rate ratios (sm="IRR") as effect measure.
# Note, the argument 'sm' is not necessary as this is the default
# in R function metainc called internally
p1 <- pairwise(list(treat1, treat2, treat3),
               list(d1, d2, d3),
               time=list(years1, years2, years3),
               studlab=ID,
               data=dietaryfat, sm="IRR")
p1

# Conduct network meta-analysis:
net1 <- netmeta(p1)
summary(net1)

# Conduct network meta-analysis using incidence rate differences
# (sm="IRD").
p2 <- pairwise(list(treat1, treat2, treat3),
               list(d1, d2, d3),
               time=list(years1, years2, years3),
               studlab=ID,
               data=dietaryfat, sm="IRD")
net2 <- netmeta(p2)
summary(net2)

# Draw network graph
netgraph(net1, points=TRUE, cex.points=3, cex=1.25)
tname <- c("Control","Diet", "Diet 2")
netgraph(net1, points=TRUE, cex.points=3, cex=1.25, labels=tname)
}

\keyword{datasets}
