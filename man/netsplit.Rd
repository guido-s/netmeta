\name{netsplit}
\alias{netsplit}
\alias{print.netsplit}

\title{Split direct and indirect evidence in network meta-analysis}

\description{
  Methods to split network estimates into the contribution of direct and
  indirect evidence and to test for local inconsistency in network
  meta-analysis.
}

\usage{
netsplit(x, method, upper = TRUE,
         reference.group = x$reference.group,
         baseline.reference = x$baseline.reference,
         sep.trts = x$sep.trts, quote.trts = "",
         tol.direct = 0.0005, warn = FALSE)

\method{print}{netsplit}(x,
      comb.fixed = x$comb.fixed,
      comb.random = x$comb.random,
      show = "all",
      overall = TRUE,
      ci = FALSE,
      test = show \%in\% c("all", "with.direct", "both"),
      digits = gs("digits"),
      digits.zval = gs("digits.zval"),
      digits.pval = gs("digits.pval"),
      digits.prop = max(gs("digits.pval")-2, 2),
      text.NA = ".", backtransf = x$backtransf,
      scientific.pval = gs("scientific.pval"),
      big.mark = gs("big.mark"),
      legend = TRUE,
      ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta} or \code{netsplit}.}
  \item{method}{A character string indicating which method to split
    direct and indirect evidence is to be used. Either
    \code{"Back-calculation"} or \code{"SIDDE"}, can be abbreviated. See
    Details.}
  \item{upper}{A logical indicating whether treatment comparisons
    should be selected from the lower or upper triangle of the
    treatment effect matrices (see list elements \code{TE.fixed} and
    \code{TE.random} in the \code{netmeta} object).}
  \item{reference.group}{Reference group.}
  \item{baseline.reference}{A logical indicating whether results
    should be expressed as comparisons of other treatments versus the
    reference treatment or vice versa. This argument is only
    considered if \code{reference.group} is not equal to \code{""}.} 
  \item{sep.trts}{A character string used in comparison names as
    separator between treatment labels, e.g., " vs ".}
  \item{quote.trts}{A character used to print around treatment
    labels.}
  \item{tol.direct}{A numeric defining the maximum deviation of the
    direct evidence proportion from 0 or 1 to classify a comparison as
    providing only indirect or direct evidence, respectively.}
  \item{comb.fixed}{A logical indicating whether results for fixed
    effect model should be printed.}
  \item{comb.random}{A logical indicating whether results for random
    effects model should be printed.}
  \item{show}{A character string indicating which comparisons should
    be printed (see Details).}
  \item{overall}{A logical indicating whether estimates from network
    meta-analyis should be printed in addition to direct and indirect
    estimates.}
  \item{ci}{A logical indicating whether confidence intervals should
    be printed in addition to treatment estimates.}  
  \item{test}{A logical indicating whether results of a test comparing
    direct and indirect estimates should be printed.}
  \item{digits}{Minimal number of significant digits, see
    \code{print.default}.}
  \item{digits.zval}{Minimal number of significant digits for z-value
    of test of agreement between direct and indirect evidence, see
    \code{print.default}.}
  \item{digits.pval}{Minimal number of significant digits for p-value
    of test of agreement between direct and indirect evidence, see
    \code{print.default}.}
  \item{digits.prop}{Minimal number of significant digits for direct
    evidence proportions, see \code{print.default}.}
  \item{text.NA}{A character string specifying text printed for
    missing values.}
  \item{backtransf}{A logical indicating whether printed results
    should be back transformed. For example, if
    \code{backtransf=TRUE}, results for \code{sm="OR"} are printed as
    odds ratios rather than log odds ratios.}
  \item{scientific.pval}{A logical specifying whether p-values should
    be printed in scientific notation, e.g., 1.2345e-01 instead of
    0.12345.}
  \item{big.mark}{A character used as thousands separator.}
  \item{legend}{A logical indicating whether a legend show be
    printed.}
  \item{warn}{A logical indicating whether warnings should be printed.}
  \item{...}{Additional arguments (ignored at the moment)}
}

\details{
  A comparison of direct and indirect treatment estimates can serve as
  check for consistency of network meta-analysis (Dias et al.,
  2010).

  This function provides two methods to derive indirect estimates:
  \itemize{
   \item Back-calculation method deriving indirect estimates from direct
     pairwise comparisons and network estimates. The direct evidence
     proportion as described in König et al. (2013) is used in the
     calculation of the indirect evidence
   \item Separate Indirect from Direct Design Evidence (SIDDE)
  }

  Note, for the back-calculation method, indirect treatment estimates
  are already calculated in \code{\link{netmeta}} and this function
  combines and prints these estimates in a user-friendly way.

  Argument \code{show} determines which comparisons are printed:
  \tabular{ll}{
  \dQuote{all} \tab All comparisons \cr
  \dQuote{both} \tab Only comparisons contributing both direct and indirect evidence \cr
  \dQuote{with.direct} \tab Comparisons providing direct evidence \cr
  \dQuote{direct.only} \tab Comparisons providing only direct evidence \cr
  \dQuote{indirect.only} \tab Comparisons providing only indirect evidence \cr
  }
}

\value{
  An object of class \code{netsplit} with corresponding
  \code{print} function. The object is a list containing the following
  components:
  \item{comb.fixed, comb.random}{As defined above.}
  \item{comparison}{A vector with treatment comparisons.}
  \item{prop.fixed, prop.random}{A vector with direct evidence
    proportions (fixed effect / random effects model).}
  \item{fixed, random}{Results of network meta-analysis (fixed effect
    / random effects model), i.e., data frame with columns comparison,
    TE, seTE, lower, upper, z, and p.}  
  \item{direct.fixed, direct.random}{Network meta-analysis results
    based on direct evidence (fixed effect / random effects model),
    i.e., data frame with columns comparison, TE, seTE, lower, upper,
    z, and p.}
  \item{indirect.fixed, indirect.random}{Network meta-analysis results
    based on indirect evidence (fixed effect / random effects model),
    i.e., data frame with columns comparison, TE, seTE, lower, upper,
    z, and p.}
  \item{compare.fixed, compare.random}{Comparison of direct and
    indirect evidence in network meta-analysis (fixed effect / random
    effects model), i.e., data frame with columns comparison, TE,
    seTE, lower, upper, z, and p.}
  \item{sm}{A character string indicating underlying summary measure}  
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{version}{Version of R package netmeta used to create object.}
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}, Gerta Rücker
\email{ruecker@imbi.uni-freiburg.de}}

\references{
  Dias S, Welton NJ, Caldwell DM, Ades AE (2010),
  Checking consistency in mixed treatment comparison meta-analysis.
  \emph{Statistics in Medicine},
  \bold{29}, 932--44

  König J, Krahn U, Binder H (2013),
  Visualizing the flow of evidence in network meta-analysis and
  characterizing mixed treatment comparisons.
  \emph{Statistics in Medicine},
  \bold{32}(30), 5414--29
  
  Puhan MA, Schünemann HJ, Murad MH, et al. (2014),  
  A GRADE working group approach for rating the quality of treatment
  effect estimates from network meta-analysis.
  \emph{British Medical Journal},
  \bold{349}, g5630
}

\seealso{\code{\link{netmeta}}, \code{\link{netmeasures}}}

\examples{
data(Woods2010)
#
p1 <- pairwise(treatment, event=r, n=N,
               studlab=author, data=Woods2010, sm="OR")
#
net1 <- netmeta(p1)
#
print(netsplit(net1), digits=2)
print(netsplit(net1), digits=2,
      backtransf=FALSE, comb.random=TRUE)

data(Senn2013)
#
net2 <- netmeta(TE, seTE, treat1.long, treat2.long,
                studlab, data=Senn2013,
                comb.random=TRUE)
#
print(netsplit(net2), digits=2)
# Layout of Puhan et al. (2014), Table 1
print(netsplit(net2), digits=2, ci=TRUE, test=FALSE)
}
