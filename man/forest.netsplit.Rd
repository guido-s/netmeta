\name{forest.netsplit}

\alias{forest.netsplit}

\title{Forest plot for direct and indirect evidence}

\description{
  Forest plot to show direct and indirect evidence in network
  meta-analysis. Furthermore, estimates from network meta-analysis as
  well as prediction intervals can be printed.
}

\usage{
\method{forest}{netsplit}(x,
       pooled=ifelse(x$comb.random, "random", "fixed"),
       showall = FALSE,
       subgroup = "comparison",
       overall = TRUE, direct = TRUE, indirect = TRUE,
       prediction = x$prediction,
       text.overall = "Network estimate",
       text.direct = "Direct estimate",
       text.indirect = "Indirect estimate",
       text.predict = "Prediction interval",
       type.overall, type.direct, type.indirect,
       col.square = "gray", col.square.lines = col.square,
       col.inside = "white",
       col.diamond = "gray", col.diamond.lines = "black",
       col.predict = "red", col.predict.lines = "black",
       equal.size = FALSE,
       leftcols, leftlabs,
       rightcols = c("effect", "ci"), rightlabs = NULL,
       digits = gs("digits.forest"),
       digits.prop = max(gs("digits.pval") - 2, 2),
       backtransf = x$backtransf,
       lab.NA = "", smlab, ...)
}

\arguments{
  \item{x}{An object of class \code{netsplit}.}
  \item{pooled}{A character string indicating whether results for
    the fixed effect (\code{"fixed"}) or random effects model
    (\code{"random"}) should be plotted. Can be abbreviated.}
  \item{showall}{A logical indicating whether all comparisons
    (default) or only comparisons contributing both direct and
    indirect evidence should be shown in forest plot.}
  \item{overall}{A logical indicating whether network meta-analysis
    estimates should be printed.}
  \item{direct}{A logical indicating whether direct estimates should
    be printed.}
  \item{indirect}{A logical indicating whether indirect estimates
    should be printed.}
  \item{prediction}{A logical indicating whether prediction intervals
    should be printed.}
  \item{subgroup}{A character string indicating which layout should be
    used in forest plot: subgroups by comparisons
    (\code{"comparison"}) or subgroups by estimates
    (\code{"estimate"}). Can be abbreviated.}
  \item{text.overall}{A character string used in the plot to label the
    network estimates.}    
  \item{text.direct}{A character string used in the plot to label the
    direct estimates.}
  \item{text.indirect}{A character string used in the plot to label
    the indirect estimates.}
  \item{text.predict}{A character string used in the plot to label the
    prediction interval.}
  \item{type.overall}{A character string specifying how to plot
    treatment effects and confidence intervals for the overall network
    evidence.}
  \item{type.direct}{A character string specifying how to plot
    treatment effects and confidence intervals for the direct
    evidence.}
  \item{type.indirect}{A character string specifying how to plot
    treatment effects and confidence intervals for the indirect
    evidence.}
  \item{col.square}{The colour for squares.}
  \item{col.square.lines}{The colour for the outer lines of squares.}
  \item{col.inside}{The colour for results and confidence limits if
    confidence limits are completely within squares squares.}
  \item{col.diamond}{The colour of diamonds.}
  \item{col.diamond.lines}{The colour of the outer lines of diamonds.}
  \item{col.predict}{Background colour of prediction intervals.}
  \item{col.predict.lines}{Colour of outer lines of prediction intervals.}
  \item{equal.size}{A logical indicating whether all squares should be
    of equal size. Otherwise, the square size is proportional to the
    precision of estimates.}
  \item{leftcols}{A character vector specifying columns to be plotted
    on the left side of the forest plot (see Details).}
  \item{leftlabs}{A character vector specifying labels for columns on
    left side of the forest plot.}
  \item{rightcols}{A character vector specifying columns to be plotted
    on the right side of the forest plot (see Details).} 
  \item{rightlabs}{A character vector specifying labels for columns on
    right side of the forest plot.}
  \item{digits}{Minimal number of significant digits for treatment
    effects and confidence intervals, see \code{print.default}.}  
  \item{digits.prop}{Minimal number of significant digits for the
    direct evidence proportion.}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in forest plots. If \code{backtransf=TRUE},
    results for \code{sm="OR"} are presented as odds ratios rather
    than log odds ratios, for example.}
  \item{lab.NA}{A character string to label missing values.}
  \item{smlab}{A label printed at top of figure. By default, text
    indicating either fixed effect or random effects model is
    printed.}
  \item{\dots}{Additional arguments for \code{\link{forest.meta}}
    function.}
}

\details{
  A forest plot, also called confidence interval plot, is drawn in the
  active graphics window.

  The arguments \code{leftcols} and \code{rightcols} can be used to
  specify columns which are plotted on the left and right side of the
  forest plot, respectively. If argument \code{rightcols} is
  \code{FALSE}, no columns will be plotted on the right side.
  
  If direct estimates are included in the forest plot
  (\code{direct=TRUE}, default), the following columns will be printed
  on the left side of the forest plot: the comparisons (column
  \code{"studlab"} in \code{\link{forest.meta}}), number of pairwise
  comparisons (\code{"k"}), and direct evidence proportion
  (\code{"k"}).

  If direct estimates are not included in the forest plot
  (\code{direct=FALSE}), only the comparisons (\code{"studlab"}) are
  printed on the left side of the forest plot.
     
  For more information see help page of \code{\link{forest.meta}}
  function.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{forest.meta}}}

\examples{
data(Senn2013)

\dontrun{
#
# Fixed effect model (default)
#
net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD")

forest(net1, ref="plac")

forest(net1, xlim=c(-1.5,1), ref="plac",
       xlab="HbA1c difference", rightcols=FALSE)
}

#
# Random effects effect model
#
net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD", comb.random=TRUE)

forest(net2, xlim=c(-1.5,1), ref="plac",
       xlab="HbA1c difference")

#
# Add column with P-Scores on right side of forest plot
#
forest(net2, xlim=c(-1.5,1), ref="plac",
       xlab="HbA1c difference",
       rightcols=c("effect", "ci", "Pscore"),
       rightlabs="P-Score",
       just.addcols="right")

\dontrun{
#
# Add column with P-Scores on left side of forest plot
#
forest(net2, xlim=c(-1.5,1), ref="plac",
       xlab="HbA1c difference",
       leftcols=c("studlab", "Pscore"),
       leftlabs=c("Treatment", "P-Score"),
       just.addcols="right")

#
# Sort forest plot by descending P-Score
#
forest(net2, xlim=c(-1.5,1), ref="plac",
       xlab="HbA1c difference",
       rightcols=c("effect", "ci", "Pscore"),
       rightlabs="P-Score",
       just.addcols="right",
       sortvar=-Pscore)

#
# Drop reference group and sort by and print number of studies
# with direct treatment comparisons
#
forest(net2, xlim=c(-1.5,1), ref="plac",
       xlab="HbA1c difference",
       leftcols=c("studlab", "k"),
       leftlabs=c("Contrast\nto Placebo", "Direct\nComparisons"),
       sortvar=-k,
       drop=TRUE,
       smlab="Random Effects Model")
}
}

\keyword{hplot}
