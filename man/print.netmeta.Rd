\name{print.netmeta}

\alias{print.netmeta}
\alias{summary.netmeta}
\alias{print.summary.netmeta}

\title{Print and summary method for objects of class netmeta}

\description{
  Print and summary method for objects of class \code{netmeta}.
}

\usage{
\method{print}{netmeta}(x, sortvar, level=x$level, level.comb=x$level.comb,
        comb.fixed=x$comb.fixed, comb.random=x$comb.random,
        reference.group=x$reference.group, all.treatments=x$all.treatments,
        details=TRUE, ma=TRUE, logscale=FALSE,
        digits=max(4, .Options$digits - 3), ...)

\method{summary}{netmeta}(object,
        level=object$level, level.comb=object$level.comb,
        comb.fixed=object$comb.fixed, comb.random=object$comb.random,
        reference.group=object$reference.group, all.treatments=object$all.treatments,
        warn=object$warn, ...)

\method{print}{summary.netmeta}(x, comb.fixed=x$comb.fixed, comb.random=x$comb.random,
        reference.group=x$reference.group, all.treatments=x$all.treatments,
        logscale=FALSE, header=TRUE, digits=max(3, .Options$digits - 3), ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta} or \code{summary.netmeta}.}
  \item{object}{An object of class \code{netmeta}.}
  \item{sortvar}{An optional vector used to sort individual studies
    (must be of same length as \code{x$TE}).}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
  \item{reference.group}{Reference group.}
  \item{all.treatments}{A logical or value \code{"NULL"}. If
        \code{TRUE}, matrices with all treatment effects, and confidence
        limits will be printed.}
  \item{details}{A logical indicating whether further details of
    individual studies should be printed.}
  \item{ma}{A logical indicating whether the summary results of the
    meta-analysis should be printed.}
  \item{logscale}{A logical indicating whether results for summary
    measures 'RR', 'OR', 'HR', or 'PLN' will be printed on logarithmic
    scale.}
  \item{header}{A logical indicating whether information on title of
    meta-analysis, comparison and outcome should be printed at the
    beginning of the printout.}
  \item{digits}{Minimal number of significant digits, see \code{print.default}.}
  \item{warn}{A logical indicating whether the use of
    \code{summary.meta} in connection with \code{metacum} or
    \code{metainf} should result in a warning.}
  \item{\dots}{Additional arguments}
}

\details{
  TODO ...
}

\value{
  A list is returned by the function \code{summary.netmeta} with the
  following elements:
  \item{comparison}{Results for pairwise comparisons (a list with
    elements TE, seTE, lower, upper, z, p, level, df, studlab, treat1,
    treat2).}
  \item{comparison.nma.fixed}{TODO: Results for pairwise comparisons
    based on fixed effect model (a list with elements TE, seTE, lower,
    upper, z, p, level, df, studlab, treat1, treat2, leverage).}
  \item{comparison.nma.random}{TODO: Results for pairwise comparisons
    based on random effects model (a list with elements TE, seTE,
    lower, upper, z, p, level, df, studlab, treat1, treat2).}
  \item{fixed}{TODO: Results for fixed effect model
    (a list with elements TE, seTE, lower, upper, z, p, level, df).}
  \item{random}{TODO: Results for random effects model
    (a list with elements TE, seTE, lower, upper, z, p, level, df).}
  \item{studies}{???}
  \item{narms}{???}
  \item{k}{Total number of studies.}
  \item{m}{Total number of pairwise comparisons.}
  \item{Q}{Overall heterogeneity/inconsistency statistic.}
  \item{df}{Degrees of freedom for test of heterogeneity/inconsistency.}
  \item{tau}{Square-root of between-study variance.}
  \item{I2}{I-squared.}
  \item{sm}{A character string indicating underlying summary measure.}
  \item{ci.lab}{???}
  \item{comb.fixed}{A logical indicating whether result for fixed
    effect meta-analysis should be printed.}
  \item{comb.random}{A logical indicating whether result for random
    effects meta-analysis should be printed.}
  \item{level}{The level used to calculate confidence intervals for
    individual comparisons.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{all.treatments}{A logical or value \code{"NULL"}. If
        \code{TRUE}, matrices with all treatment effects, and confidence
        limits will be printed.}
  \item{reference.group}{Reference group.}
  \item{all.treatments}{A logical or value \code{"NULL"}. If
        \code{TRUE}, matrices with all treatment effects, and confidence
        limits will be printed.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{call}{Function call.}
  \item{version}{Version of R package netmeta used to create object.}
}

\references{
  TODO - any references ???
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}}

\examples{

## TODO: Better example?

data(example1)

example1$se.effect <- sqrt(1/example1$weight)

## Fixed effect model (default)
##
nm <- netmeta(effect, se.effect, treat1, treat2, author,
              data=example1)
summary(nm)
}

\keyword{print}
