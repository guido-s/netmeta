\name{print.netmeta}

\alias{print.netmeta}
\alias{summary.netmeta}
\alias{print.summary.netmeta}

\title{Print and summary method for objects of class netmeta}

\description{
  Print and summary method for objects of class \code{netmeta}.
}

\usage{
\method{print}{netmeta}(x, sortvar,
      comb.fixed=x$comb.fixed, comb.random=x$comb.random,
      prediction=x$prediction,
      reference.group=x$reference.group,
      baseline.reference=x$baseline.reference,
      all.treatments=x$all.treatments,
      details=TRUE, ma=TRUE,
      backtransf=x$backtransf, nchar=x$nchar,
      digits=max(4, .Options$digits - 3), ...)

\method{summary}{netmeta}(object,
        comb.fixed=object$comb.fixed, comb.random=object$comb.random,
        prediction=object$prediction,
        reference.group=object$reference.group,
        baseline.reference=object$baseline.reference,
        all.treatments=object$all.treatments,
        warn=object$warn, ...)

\method{print}{summary.netmeta}(x, comb.fixed=x$comb.fixed, comb.random=x$comb.random,
      prediction=x$prediction,
      reference.group=x$reference.group,
      baseline.reference=x$baseline.reference,
      all.treatments=x$all.treatments,
      backtransf=x$backtransf, nchar=x$nchar,
      header=TRUE, digits=max(3, .Options$digits - 3), ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta} or \code{summary.netmeta}.}
  \item{object}{An object of class \code{netmeta}.}
  \item{sortvar}{An optional vector used to sort individual studies
    (must be of same length as \code{x$TE}).}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
  \item{prediction}{A logical indicating whether prediction intervals
    should be printed.}
  \item{reference.group}{Reference group.}
  \item{baseline.reference}{A logical indicating whether results
    should be expressed as comparisons of other treatments versus the
    reference treatment (default) or vice versa. This argument is only
    considered if \code{reference.group} has been specified.}
  \item{all.treatments}{A logical or value \code{"NULL"}. If
    \code{TRUE}, matrices with all treatment effects, and confidence
    limits will be printed.}
  \item{details}{A logical indicating whether further details for
    individual studies should be printed.}
  \item{ma}{A logical indicating whether summary results of
    meta-analysis should be printed.}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in printouts and forest plots. If
    \code{backtransf=TRUE}, results for \code{sm="OR"} are presented
    as odds ratios rather than log odds ratios, for example.}
  \item{nchar}{A numeric defining the minium number of characters
    used to abbreviate treatment names (see \code{\link{abbreviate}}
    for more information).}
  \item{header}{A logical indicating whether information on title of
    meta-analysis, comparison and outcome should be printed at the
    beginning of the printout.}
  \item{digits}{Minimal number of significant digits, see \code{print.default}.}
  \item{warn}{A logical indicating whether the use of
    \code{summary.meta} in connection with \code{metacum} or
    \code{metainf} should result in a warning.}
  \item{\dots}{Additional arguments.}
}

\value{
  A list is returned by the function \code{summary.netmeta} with the
  following elements:
  \item{comparison}{Results for pairwise comparisons (a list with
    elements TE, seTE, lower, upper, z, p, level, df, studlab, treat1,
    treat2).}
  \item{comparison.nma.fixed}{Results for pairwise comparisons based
    on fixed effect model (a list with elements TE, seTE, lower,
    upper, z, p, level, df, studlab, treat1, treat2, leverage).}
  \item{comparison.nma.random}{Results for pairwise comparisons based
    on random effects model (a list with elements TE, seTE, lower,
    upper, z, p, level, df, studlab, treat1, treat2).}
  \item{fixed}{Results for fixed effect model (a list with elements
    TE, seTE, lower, upper, z, p, level, df).}
  \item{random}{Results for random effects model (a list with elements
    TE, seTE, lower, upper, z, p, level, df).}
  \item{studies}{Study labels coerced into a factor with its levels
    sorted alphabetically.}
  \item{narms}{Number of arms for each study.}
  \item{k}{Total number of studies.}
  \item{m}{Total number of pairwise comparisons.}
  \item{n}{Total number of treatments.}
  \item{Q}{Overall heterogeneity / inconsistency statistic.}
  \item{df.Q}{Degrees of freedom for test of heterogeneity /
    inconsistency.}
  \item{pval.Q}{P-value for test of heterogeneity / inconsistency.}
  \item{I2}{I-squared.}
  \item{tau}{Square-root of between-study variance.}
  \item{Q.heterogeneity}{Overall heterogeneity statistic.}
  \item{df.Q.heterogeneity}{Degrees of freedom for test of overall
    heterogeneity.}
  \item{pval.Q.heterogeneity}{P-value for test of overall heterogeneity.}
  \item{Q.inconsistency}{Overall inconsistency statistic.}
  \item{df.Q.inconsistency}{Degrees of freedom for test of overall
    inconsistency.}
  \item{pval.Q.inconsistency}{P-value for test of overall inconsistency.}
  \item{sm}{A character string indicating underlying summary measure.}
  \item{ci.lab}{Label for confidence interval.}
  \item{comb.fixed}{A logical indicating whether result for fixed
    effect meta-analysis should be printed.}
  \item{comb.random}{A logical indicating whether result for random
    effects meta-analysis should be printed.}
  \item{prediction}{A logical indicating whether prediction intervals
    should be printed.}
  \item{seq}{A character specifying the sequence of treatments.}
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

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}}

\examples{
data(Senn2013)

#
# Fixed effect model (default)
#
net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD")
print(net1, ref="plac", digits=3)
summary(net1)

#
# Random effects model
#
net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD", comb.random=TRUE)
print(net2, ref="plac", digits=3)
summary(net2)
}

\keyword{print}
