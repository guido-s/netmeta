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
    less interested are included in data frame.}
  \item{...}{other arguments}
}
\value{
  A data frame is returned by the function \code{as.data.frame}.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}
\seealso{\code{\link{netmeta}}}

\examples{
data(example1)

example1$se.effect <- sqrt(1/example1$weight)

nm <- netmeta(effect, se.effect, treat1, treat2, author,
              data=example1)

as.data.frame(nm)

as.data.frame(nm, details=TRUE)
}
\keyword{print}
