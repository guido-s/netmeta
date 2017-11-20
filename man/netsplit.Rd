\name{netsplit}
\alias{netsplit}
\alias{print.netsplit}

\title{Split direct and indirect evidence in network meta-analysis}

\description{
  Back-calculation method to split contribution of direct and indirect
  evidence in network meta-analysis.
}

\usage{
netsplit(x, lower = FALSE,
         reference.group = x$reference.group,
         baseline.reference = x$baseline.reference,
         sep.trts = x$sep.trts, quote = "")

\method{print}{netsplit}(x,
      comb.fixed = x$comb.fixed,
      comb.random = x$comb.random,
      showall = TRUE,
      overall = TRUE,
      ci = FALSE,
      test = TRUE,
      digits = gs("digits"),
      digits.zval = gs("digits.zval"),
      digits.pval = gs("digits.pval"),
      text.NA = ".", backtransf = x$backtransf,
      ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta} or \code{netsplit}.}
  \item{lower}{A logical indicating whether treatment comparisons
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
  \item{quote}{A character used to print around treatment labels.}
  \item{comb.fixed}{A logical indicating whether results for fixed
    effect model should be printed.}
  \item{comb.random}{A logical indicating whether results for random
    effects model should be printed.}
  \item{showall}{A logical indicating whether all comparisons
    (default) or only comparisons contributing both direct and
    indirect evidence should be printed.} 
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
  \item{backtransf}{A logical indicating whether printed results
    should be back transformed. For example, if
    \code{backtransf=TRUE}, results for \code{sm="OR"} are printed as
    odds ratios rather than log odds ratios.}
  \item{text.NA}{A character string specifying text printed for
    missing values.}
  \item{...}{Additional arguments (ignored at the moment)}
}

\details{
  Direct and indirect treatment estimates are calculated in
  \code{\link{netmeta}}. This function combines and prints these
  estimates in a user-friendly way.

  A comparison of direct and indirect treatment estimates can serve as
  check for consistency of network meta-analysis (Dias et al.,
  2010). This function provides a back-calculation method to derive
  indirect estimates from direct pairwise comparisons and network
  estimates. The direct evidence proportion as describe in König et
  al. (2013) is used to calculate the indirect evidence.
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
    i.e., data frame with columsn comparison, TE, seTE, lower, upper,
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
  Dias S, Welton NJ, Caldwell DM, Ades AE (2010).
  Checking consistency in mixed treatment comparison meta-analysis.
  \emph{Statistics in Medicine}, \bold{29}, 932--44.

  König J, Krahn U, Binder H (2013).
  Visualizing the flow of evidence in network meta-analysis and
  characterizing mixed treatment comparisons.
  \emph{Statistics in Medicine}, \bold{32}(30), 5414--29.
  
  Puhan MA, Schünemann HJ, Murad MH, et al. (2014).  
  A GRADE working group approach for rating the quality of treatment
  effect estimates from network meta-analysis.
  \emph{British Medical Journal}, \bold{349}, g5630
}

\seealso{\code{\link{netmeta}}, \code{\link{netmeasures}}}

\examples{
data(Woods2010)
#
p1 <- pairwise(treatment, event = r, n = N,
               studlab = author, data = Woods2010, sm = "OR")
#
net1 <- netmeta(p1)
#
print(netsplit(net1), digits = 2)
print(netsplit(net1), digits = 2,
      backtransf = FALSE, comb.random = TRUE)

data(Senn2013)
#
net2 <- netmeta(TE, seTE, treat1.long, treat2.long,
                studlab, data = Senn2013,
                comb.random = TRUE)
#
print(netsplit(net2), digits = 2)
# Layout of Puhan et al. (2014), Table 1
print(netsplit(net2), digits = 2, ci = TRUE, test = FALSE)
}
