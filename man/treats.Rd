\name{treats}
\alias{treats}

\title{Abbreviate treatment names}

\description{
  Auxiliary function to create uniquely abbreviated treatment names.
}

\usage{
treats(x, nchar.trts = 8, row = TRUE)
}

\arguments{
  \item{x}{A vector with treatment names or a matrix with treatment
    names as row and / or column names.}
  \item{nchar.trts}{Minimum number of characters to create unique
    treatment names.}
  \item{row}{A logical indicating whether row or column names should
    be used (only considered if argument \code{x} is a matrix).}
}

\details{
  This auxiliary function can be used to create uniquely abbreviated
  treatment names (and is used internally in several R functions for
  this purpose).

  Initially, to construct uniquely abbreviated treatment names,
  \code{\link{substring}} is used to extract the first
  \code{nchar.trts} characters. If these abbreviated treatment names
  are not unique, \code{\link{abbreviate}} with argument
  \code{minlength=nchar.trts} is used.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}, \code{\link{print.netmeta}},
\code{\link{print.summary.netmeta}}}

\examples{
data(Senn2013)
#
net1 <- netmeta(TE, seTE, treat1, treat2,
                studlab, data=Senn2013)

# Use matrix with fixed effect treatment estimates to create unique
# treatment names
#
# Four characters
#
treats(net1$TE.fixed, nchar.trts = 4)

# Two characters
#
treats(net1$TE.fixed, nchar.trts = 2)

# One character
#
treats(net1$TE.fixed, nchar.trts = 1)
}

