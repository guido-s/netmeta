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
    fixed effect (\code{"fixed"}) or random effects model
    (\code{"random"}) should be plotted. Can be abbreviated.}
  \item{leftcols}{A character vector specifying (additional) columns
    to be plotted on the left side of the forest plot or a logical
    value (see \code{\link{forest.meta}} help page for details).}
  \item{leftlabs}{A character vector specifying labels for
    (additional) columns on left side of the forest plot (see
    \code{\link{forest.meta}} help page for details).}
  \item{smlab}{A label printed at top of figure. By default, text
    indicating either fixed effect or random effects model is
    printed.}
  \item{\dots}{Additional arguments for \code{\link{forest.meta}}
    function.}
}

\details{
  A forest plot, also called confidence interval plot, is drawn in the
  active graphics window.

  For more information see help page of \code{\link{forest.meta}}
  function.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{forest.meta}}}

\examples{

data(Senn2013)

#### Fixed effect model (default)

net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD")

forest(net1, ref="placebo")

forest(net1, xlim=c(-1.5,1), ref="placebo",
       xlab="HbA1c difference",
       leftcols="studlab", rightcols=NULL,
       leftlabs="Contrast to placebo")

#### Random effects effect model

net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD", comb.random=TRUE)

forest(net2, xlim=c(-1.5,1), ref="placebo",
       xlab="HbA1c difference",
       leftcols="studlab", rightcols=NULL,
       leftlabs="Contrast to placebo")
}

\keyword{hplot}
