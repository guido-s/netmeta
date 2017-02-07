\name{netsplit}
\alias{netsplit}
\alias{print.netsplit}

\title{Split into direct and indirect evidence}

\description{
  ...
}

\usage{
netsplit(x)

\method{print}{netsplit}(x,
        comb.fixed = x$comb.fixed,
        comb.random = x$comb.random,
        digits = gs("digits"),
        digits.zval = gs("digits.zval"),
        digits.pval = gs("digits.pval"),
        text.NA = ".", backtransf = TRUE,
        showall = FALSE, ...)
}

\arguments{
  \item{x}{netmeta() or netsplit() object.}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
  \item{digits}{Minimal number of significant digits, see \code{print.default}.}
  \item{digits.zval}{Minimal number of significant digits for z- or
    t-value, see \code{print.default}.}
  \item{digits.pval}{Minimal number of significant digits for p-value
    of overall treatment effect, see \code{print.default}.}
  \item{backtransf}{A logical indicating whether printed results
    should be back transformed. For example, if
    \code{backtransf=TRUE}, results for \code{sm="OR"} are printed as
    odds ratios rather than log odds ratios.}
  \item{text.NA}{...}
  \item{showall}{A logical ...}
  \item{...}{Additional arguments (ignored at the moment)}
}

\value{
  An object of class \code{netsplit} with corresponding
  \code{print} function. The object is a list containing the following
  components:
  \item{...}{...}
  \item{call}{Function call.}
  \item{version}{Version of R package netmeta used to create object.}
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}, Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}}

\examples{
data(Woods2010)
#
p1 <- pairwise(treatment, event = r, n = N,
               studlab = author, data = Woods2010, sm = "OR")
#
net1 <- netmeta(TE, seTE, treat1, treat2, studlab, data = p1)
#
oldopts <- options(width = 90)
print(netsplit(net1), digits = 2)
print(netsplit(net1), digits = 2,
      backtransf = FALSE, comb.random = TRUE)
options(oldopts)

data(Senn2013)
#
net2 <- netmeta(TE, seTE, treat1, treat2, studlab, data = Senn2013,
                comb.random = TRUE)
#
print(netsplit(net2), digits = 2)
}
