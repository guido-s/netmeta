\name{print.netcomb}

\alias{print.netcomb}
\alias{summary.netcomb}
\alias{print.summary.netcomb}

\title{Print and summary method for objects of class netcomb}

\description{
  Print and summary method for objects of class \code{netcomb}.
}

\usage{
\method{print}{netcomb}(x,
      comb.fixed = x$comb.fixed,
      comb.random = x$comb.random,
      backtransf = x$backtransf,
      nchar.trts = x$nchar.trts,
      digits = gs("digits"),
      digits.zval = gs("digits.zval"),
      digits.pval = gs("digits.pval"),
      digits.pval.Q = max(gs("digits.pval.Q"), 2),
      digits.Q = gs("digits.Q"),
      scientific.pval = gs("scientific.pval"),
      big.mark = gs("big.mark"),
      ...)

\method{summary}{netcomb}(object, ...)

\method{print}{summary.netcomb}(x,
      comb.fixed = x$comb.fixed,
      comb.random = x$comb.random,
      backtransf = x$backtransf,
      nchar.trts = x$nchar.trts,
      digits = gs("digits"),
      digits.zval = gs("digits.zval"),
      digits.pval = gs("digits.pval"),
      digits.pval.Q = max(gs("digits.pval.Q"), 2),
      digits.Q = gs("digits.Q"),
      digits.tau2 = gs("digits.tau2"),
      digits.I2 = gs("digits.I2"),
      scientific.pval = gs("scientific.pval"),
      big.mark = gs("big.mark"),
      ...)
}

\arguments{
  \item{x}{An object of class \code{netcomb} or \code{summary.netcomb}.}
  \item{object}{An object of class \code{netcomb}.}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in printouts and forest plots. If
    \code{backtransf=TRUE}, results for \code{sm="OR"} are presented
    as odds ratios rather than log odds ratios, for example.}
  \item{nchar.trts}{A numeric defining the minium number of characters
    used to create unique treatment names (see Details).}
  \item{digits}{Minimal number of significant digits, see
    \code{print.default}.}
  \item{digits.zval}{Minimal number of significant digits for z- or
    t-value, see \code{print.default}.}
  \item{digits.pval}{Minimal number of significant digits for p-value
    of overall treatment effect, see \code{print.default}.}
  \item{digits.pval.Q}{Minimal number of significant digits for
    p-value of heterogeneity tests, see \code{print.default}.}
  \item{digits.Q}{Minimal number of significant digits for
    heterogeneity statistics, see \code{print.default}.}
  \item{digits.tau2}{Minimal number of significant digits for
    between-study variance, see \code{print.default}.}  
  \item{digits.I2}{Minimal number of significant digits for I-squared
    statistic, see \code{print.default}.}
  \item{scientific.pval}{A logical specifying whether p-values should
    be printed in scientific notation, e.g., 1.2345e-01 instead of
    0.12345.}
  \item{big.mark}{A character used as thousands separator.} 
  \item{\dots}{Additional arguments.}
}

\value{
  A list is returned by the function \code{summary.netcomb} with the
  same elements as a \code{\link{netcomb}}linknetcomb object.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netcomb}}}

\examples{
data(Linde2016)

# Standard Network meta-analysis
#
net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
                data = Linde2016, reference.group = "placebo",
                sm = "OR", comb.fixed = FALSE, comb.random = TRUE)

# Additive model for treatment components
#
nc1 <- netcomb(net1)

summary(nc1)
print(summary(nc1), digits = 2, digits.zval = 3)
}

\keyword{print}
