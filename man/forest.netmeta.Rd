\name{forest.netmeta}

\alias{forest.netmeta}

\title{Forest plot}

\description{
  Draws a forest plot in the active graphics window (using grid
  graphics system).
}

\usage{
\method{forest}{netmeta}(x,
          pooled=ifelse(x$comb.random, "random", "fixed"),
          reference.group=x$reference.group,
          leftcols="studlab",
          leftlabs="Treatment",
          smlab=NULL,
          ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta}.}
  \item{reference.group}{Reference group.}
  \item{pooled}{A character string indicating whether results for
    fixed effect ( \code{"fixed"}) or random effects model
    (\code{"random"}) should be plotted. Can be abbreviated.}
  \item{leftcols}{A character vector specifying (additional) columns
    to be plotted on the left side of the forest plot or a logical value
    (see Details).}
  \item{leftlabs}{A character vector specifying labels for
    (additional) columns on left side of the forest plot (see
    \code{\link{forest.meta}} help page for details).}
  \item{smlab}{A label printed at top of figure. By default, text for
    either fixed effect or random effects model printed.}
  \item{\dots}{Additional arguments for \code{\link{forest.meta}}
    function.}
}

\details{
  A forest plot, also called confidence interval plot, is drawn in the
  active graphics window. Subgroup analyses are conducted and
  displayed in the plot if \code{byvar} is not missing.

  For more information see help page of \code{\link{forest.meta}}
  function.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{forest.meta}}}

\examples{

data(example1)

example1$se.effect <- sqrt(1/example1$weight)

nm <- netmeta(effect, se.effect, treat1, treat2, author,
              data=example1)

## Fixed effect model (default)
##
forest(nm, reference.group="1")

## Random effects model
##
forest(nm, reference.group="1", pooled="random")

}

\keyword{hplot}
