\name{Senn2013}
\docType{data}
\alias{Senn2013}
\title{Network meta-analysis in diabetes}

\description{Network meta-analysis in diabetes comparing effects of a
  number of drugs on the HbA1c value.

  These data are used as an example in Senn et al. (2013) and have
  been preproced for use in R package netmeta.
}

\usage{data(Senn2013)}

\format{A data frame with the following columns:
  \describe{
  \item{TE}{Treatment effect}
  \item{seTE}{Standard error of treatment effect}
  \item{treat1}{Treatment 1}
  \item{treat2}{Treatment 2}
  \item{studlab}{Study label}
  }
}

\source{
  Senn S, Gavini F, Magrez D, and Scheen A (2013).
  Issues in performing a network meta-analysis.
  \emph{Statistical Methods in Medical Research},
  \bold{22}(2), 169--189.
  First published online 2012 Jan 3.
}

\seealso{\code{\link{netmeta}}}

\examples{
data(Senn2013)

#### Fixed effect model (default)

net1 <- netmeta(TE, seTE, treat1, treat2, studlab, data=Senn2013)
net1
net1$Q.decomp

#### Comparison with reference group

netmeta(TE, seTE, treat1, treat2, studlab, data=Senn2013, 
        reference="placebo")

#### Forest plot

forest(net1, ref="placebo")

#### Random effects model

net2 <- netmeta(TE, seTE, treat1, treat2, studlab, data=Senn2013,
                comb.random = TRUE)
net2
forest(net2, ref="placebo")
}

\keyword{datasets}
