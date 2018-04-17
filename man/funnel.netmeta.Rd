\name{funnel.netmeta}

\alias{funnel.netmeta}

\title{\sQuote{Comparison-adjusted} funnel plot}

\description{
  Draw a \sQuote{comparison-adjusted} funnel plot to assess funnel
  plot asymmetry in network meta-analysis.
}

\usage{
\method{funnel}{netmeta}(x,
       order,
       pooled = ifelse(x$comb.random, "random", "fixed"),
       ##
       xlab,
       level = x$level,
       ##
       pch,
       col = "black",
       ##
       legend = TRUE,
       linreg = FALSE,
       rank = FALSE,
       mm = FALSE,
       ##
       pos.legend = "topright",
       pos.tests = "topleft",
       ##
       text.linreg = "(Egger)",
       text.rank = "(Begg-Mazumdar)",
       text.mm = "(Thompson-Sharp)",
       ##
       sep.trts = x$sep.trts,
       nchar.trts = x$nchar.trts,
       ##
       backtransf = x$backtransf,
       digits.pval = gs("digits.pval"),
       ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta}.}
  \item{order}{A mandatory character or numerical vector specifying
    the order of treatments (see Details).}
  \item{pooled}{A character string indicating whether results for the
    fixed effect (\code{"fixed"}) or random effects model
    (\code{"random"}) should be plotted. Can be abbreviated.}
  \item{xlab}{A label for the x-axis.}
  \item{level}{The confidence level utilised in the plot. For the
    funnel plot, confidence limits are not drawn if
    \code{yaxis="size"}.}
  \item{pch}{The plotting symbol(s) used for individual studies within
    direct comparisons.}
  \item{col}{The colour(s) used for individual studies within direct
    comparisons.}
  \item{legend}{A logical indicating whether a legend with information
    on direct comparisons should be added to the plot.}
  \item{linreg}{A logical indicating whether result of linear
    regression test for funnel plot asymmetry should be added to
    plot.}
  \item{rank}{A logical indicating whether result of rank test for
    funnel plot asymmetry should be added to plot.}
  \item{mm}{A logical indicating whether result of linear regression
    test for funnel plot asymmetry allowing for between-study
    heterogeneity should be added to plot.}
  \item{pos.legend}{The position of the legend describing plotting
    symbols and colours for direct comparisons.}
  \item{pos.tests}{The position of results for test(s) of funnel plot
    asymmetry.}
  \item{text.linreg}{A character string used in the plot to label the
    linear regression test for funnel plot asymmetry.}
  \item{text.rank}{A character string used in the plot to label the
    rank test for funnel plot asymmetry.}
  \item{text.mm}{A character string used in the plot to label the
    linear regression test for funnel plot asymmetry allowing for
    between-study heterogeneity.}
  \item{sep.trts}{A character used in comparison names as separator
    between treatment labels.}
  \item{nchar.trts}{A numeric defining the minium number of characters
    used to create unique treatment names (see
    \code{\link{netmeta}}).}
  \item{backtransf}{A logical indicating whether results for relative
    summary measures (argument \code{sm} equal to \code{"OR"},
    \code{"RR"}, \code{"HR"}, or \code{"IRR"}) should be back
    transformed in funnel plots. If \code{backtransf=TRUE}, results
    for \code{sm="OR"} are printed as odds ratios rather than log odds
    ratios, for example.}
  \item{digits.pval}{Minimal number of significant digits for p-value
    of test(s) for funnel plot asymmetry.}
  \item{\dots}{Additional graphical arguments passed as arguments to
    \code{\link{funnel.meta}}.}
}

\details{
  A \sQuote{comparison-adjusted} funnel plot (Chaimani et al., 2013)
  is drawn in the active graphics window.

  Argument \code{order} is mandatory to determine the order of
  treatments (Chaimani et al., 2013):
  
  \emph{\dQuote{Before using this plot, investigators should order the
    treatments in a meaningful way and make assumptions about how
    small studies differ from large ones. For example, if they
    anticipate that newer treatments are favored in small trials, then
    they could name the treatments from oldest to newest so that all
    comparisons refer to \sQuote{old versus new intervention}. Other
    possibilities include defining the comparisons so that all refer
    to an active treatment versus placebo or sponsored versus
    non-sponsored intervention.}}

  The treatments can be either in increasing or decreasing order. 
  
  In the funnel plot, if \code{yaxis} is \code{"se"}, the standard error
  of the treatment estimates is plotted on the y-axis which is likely to
  be the best choice (Sterne & Egger, 2001). Other possible choices for
  \code{yaxis} are \code{"invvar"} (inverse of the variance),
  \code{"invse"} (inverse of the standard error), and \code{"size"}
  (study size).
}

\value{
  A data frame with the following columns:

  \item{studlab}{Study label.}
  \item{treat1}{Label/Number for first treatment.}
  \item{treat2}{Label/Number for second treatment.}
  \item{comparison}{Treatment comparison.}
  \item{TE}{Estimate of treatment effect, e.g., log odds ratio.}
  \item{TE.direct}{Pooled estimate from direct evidence.}
  \item{TE.adj}{\sQuote{Comparison-adjusted} treatment effect (TE - TE.direct).}
  \item{seTE}{Standard error of treatment estimate.}
  \item{pch}{Plotting symbol(s).}
  \item{col}{Colour of plotting symbol(s).}
}

\references{
  Chaimani A, Higgins JP, Mavridis D, Spyridonos P, Salanti G (2013),
  Graphical tools for network meta-analysis in STATA.
  PLOS ONE, \bold{8}, e76654.

  Sterne JAC & Egger M (2001),
  Funnel plots for detecting bias in meta-analysis: Guidelines on choice
  of axis. \emph{Journal of Clinical Epidemiology}, \bold{54},
  1046--1055.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}, \code{\link{funnel.meta}}, \code{\link{metabias}}}

\examples{
data(Senn2013)
net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data = Senn2013, sm = "MD")

#
# 'Comparison-adjusted' funnel plot not created as argument 'order' is
# missing
#
funnel(net1)

#
# (Non-sensical) alphabetic order of treatments with placebo as last
# treatment
#
ord <- c("a", "b", "me", "mi", "pi", "r", "si", "su", "v", "pl")
funnel(net1, order = ord)

#
# Add results for tests of funnel plot asymmetry and use different
# plotting symbols and colours
#
funnel(net1, order = ord,
       pch = rep(c(15:18, 1), 3), col = 1:3,
       linreg = TRUE, rank = TRUE, mm = TRUE, digits.pval = 2)

#
# Same results for tests of funnel plot asymmetry using reversed order
# of treatments
#
funnel(net1, order = rev(ord),
       pch = rep(c(15:18, 1), 3), col = 1:3,
       linreg = TRUE, rank = TRUE, mm = TRUE, digits.pval = 2)

#
# Calculate tests for funnel plot asymmetry
#
f1 <- funnel(net1, order = ord,
             pch = rep(c(15:18, 1), 3), col = 1:3,
             linreg = TRUE, rank = TRUE, mm = TRUE)
#
metabias(metagen(TE.adj, seTE, data = f1))
metabias(metagen(TE.adj, seTE, data = f1), method = "rank")
metabias(metagen(TE.adj, seTE, data = f1), method = "mm")
}

\keyword{hplot}
