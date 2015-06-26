\name{netrank}
\alias{netrank}
\alias{print.netrank}

\title{Frequentist method to rank treatments in network}

\description{
  Ranking treatments in frequentist network meta-analysis without resampling methods.
}

\usage{
netrank(x, small.values="good")

\method{print}{netrank}(x, sort=TRUE,
                        digits=max(4, .Options$digits - 3), ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta} (netrank function) or
    \code{netrank} (print function).}
  \item{small.values}{A character string specifying whether small
    treatment effects indicate a beneficial (\code{"good"}) or harmful
    (\code{"bad"}) effect, can be abbreviated.}
  \item{sort}{A logical indicating whether printout should be sorted
    by decreasing P-score.}
  \item{digits}{Minimal number of significant digits, see \code{print.default}.}
  \item{\dots}{Additional arguments passed on to print.data.frame
    function (used internally).}
}

\details{
  ... Gerta ...
}

\value{
  An object of class \code{c("netrank")} with corresponding
  \code{print} function. The object is a list containing the following
  components:
  \item{Pscore}{A named numeric vector with P-scores.}
  \item{Pmatrix}{Numeric matrix based on pairwise one-sided p-values.}
  \item{small.values, x}{As defined above.}
  \item{version}{Version of R package netmeta used to create object.}
}

\references{
  Rücker G & Schwarzer G (2015),
  Ranking treatments in frequentist network meta-analysis works
  without resampling methods.
  Submitted.

  Salanti G (2012),
  Indirect and mixed-treatment comparison, network, or
  multiple-treatments meta-analysis: Many names, many benefits, many
  concerns for the next generation evidence synthesis tool.
  \emph{Research Synthesis Methods},
  \bold{3}(2), 80--97. 
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}}

\examples{

data(Senn2013)

net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD")

nr1 <- netrank(net1)

nr1

print(nr1, sort=FALSE)
}

\keyword{print}
