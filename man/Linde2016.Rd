\name{Linde2016}
\docType{data}
\alias{Linde2016}
\title{Network meta-analysis of primary care depression treatments}

\description{
Network meta-analysis of 22 treatments (including placebo and usual
care) for the primary care of depression.
}

\usage{data(Linde2016)}

\format{A data frame with the following columns:
  \describe{
  \item{id}{Study ID}
  \item{lnOR}{Response after treatment (log odds ratio)}
  \item{selnOR}{Standard error of log odds ratio}
  \item{treat1}{First treatment}
  \item{treat2}{Second treatment}
  }
}

\source{
  Linde K, Rücker G, Schneider A et al. (2016),
  Questionable assumptions hampered interpretation of a network
  meta-analysis of primary care depression treatments.
  \emph{Journal of Clinical Epidemiology}
  \bold{71}, 86--96
}

\seealso{\code{\link{netmeta}}, \code{\link{netcomb}}}

\examples{
data(Linde2016)

# Standard Network meta-analysis
#
net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
                data = Linde2016, reference.group = "placebo",
                sm = "OR", comb.fixed = FALSE, comb.random = TRUE,
                nchar = 6)
#
summary(net1)
}

\keyword{datasets}
