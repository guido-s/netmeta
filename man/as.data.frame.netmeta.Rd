\name{as.data.frame.netmeta}
\alias{as.data.frame.netmeta}
\title{Additional functions for objects of class netmeta}
\description{
  The \code{as.data.frame} method returns a data frame containing
  information on individual studies, e.g., estimated treatment effect
  and its standard error.
}
\usage{
\method{as.data.frame}{netmeta}(x, row.names=NULL, optional=FALSE, details=FALSE, ...)
}
\arguments{
  \item{x}{An object of class \code{netmeta}.}
  \item{row.names}{\code{NULL} or a character vector giving the row
  names for the data frame.}
  \item{optional}{A logical. If \code{TRUE}, setting row names and
    converting column names (to syntactic names) is optional.}
  \item{details}{A logical. If \code{TRUE}, additional variables of
    less interest are included in data frame.}
  \item{\dots}{Additional arguments.}
}
\value{
  A data frame is returned by the function \code{as.data.frame}.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}
\seealso{\code{\link{netmeta}}}

\examples{
data(Senn2013)

#
# Conduct network meta-analysis
#
net1 <- netmeta(TE, seTE, treat1, treat2, studlab, data=Senn2013, sm="MD")

as.data.frame(net1)

as.data.frame(net1, details=TRUE)
}

\keyword{print}
