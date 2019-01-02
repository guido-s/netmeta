\name{netdistance}

\alias{netdistance}

\title{Calculate distance matrix for an adjacency matrix}

\description{
  Calculate distance matrix for an adjacency matrix based on distance
  algorithm by Müller et al. (1987).
}

\usage{
netdistance(x)
}

\arguments{
  \item{x}{Either a netmeta object or an adjacency matrix.}
}

\references{
  Müller WR, Szymanski K, Knop JV, and Trinajstic N (1987):
  An algorithm for construction of the molecular distance matrix.
  \emph{Journal of Computational Chemistry},
  \bold{8}, 170--73
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}, \code{\link{netconnection}}}

\examples{
data(Senn2013)

net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data = Senn2013, sm = "MD")

netdistance(net1)
netdistance(net1$A.matrix)
}
