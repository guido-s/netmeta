\name{print.decomp.design}

\alias{print.decomp.design}

\title{Print method for objects of class decomp.design}

\description{
  Print and summary method for objects of class \code{decomp.design}.
}

\usage{
\method{print}{decomp.design}(x, digits=2, ...)
}

\arguments{
  \item{x}{An object of class \code{decomp.design}.}
  \item{digits}{Minimal number of significant digits for Q statistics, see \code{print.default}.}
  \item{\dots}{Additional arguments (ignored at the moment).}
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}, Ulrike Krahn \email{ulrike.krahn@bayer.com}}

\seealso{\code{\link{decomp.design}}}

\examples{
data(Senn2013)

net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD")
print(decomp.design(net1))
}

\keyword{print}
