\name{plot.netposet}

\alias{plot.netposet}

\title{Scatter plot or biplot showing partially order of treatment ranks}

\description{
  This function generates a scatter plot or biplot of P-scores with an
  overlay describing partial order of treatment ranks.
}

\usage{
\method{plot}{netposet}(x,
     plottype = "scatter",
     pooled=ifelse(x$comb.random, "random", "fixed"),
     dim = "2d",
     sel.x = 1, sel.y = 2, sel.z = 3,
     cex = 1, col = "black",
     cex.text = cex, col.text = col,
     adj.x = 0, adj.y = 1,
     offset.x = 0.005, offset.y = -0.005,
     pch = NULL, cex.points = cex, col.points = col,
     col.lines = "black", lty.lines = 1, lwd.lines = 1,
     arrows = FALSE,
     length = 0.05,
     grid = TRUE,
     col.grid = "gray", lty.grid = 2, lwd.grid = 1,
     ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta} (mandatory).}
  \item{plottype}{A character string indicating whether a scatter plot
    or biplot should be produced, either \code{"scatter"} or
    \code{"biplot"}. Can be abbreviated.}
  \item{pooled}{A character string indicating whether scatter plot
    should be drawn for fixed effect (\code{"fixed"}) or random
    effects model (\code{"random"}). Can be abbreviated.}
  \item{dim}{A character string indicating whether a 2- or
    3-dimensional plot should be produced, either \code{"2d"} or
    \code{"3d"}. Can be abbreviated.}
  \item{sel.x}{A numeric specifying number of outcome to use for the
    x-axis in a scatterplot (argument is not considered for a biplot).}
  \item{sel.y}{A numeric specifying number of outcome to use for the
    y-axis in a scatterplot (argument is not considered for a biplot).}
  \item{sel.z}{A numeric specifying number of outcome to use for the
    z-axis in a scatterplot (argument is not considered for a biplot).}
  \item{cex}{The magnification to be used for treatment labels and points.}
  \item{col}{Colour(s) of treatment labels and points.}
  \item{cex.text}{The magnification to be used for treatment labels.}
  \item{col.text}{Colour(s) of treatment labels.}
  \item{adj.x}{Value(s) in [0, 1] to specify adjustment of treatment
    labels on x-axis (only considered in 2-D plots); see
    \code{\link[graphics]{text}}.}
  \item{adj.y}{Value(s) in [0, 1] to specify adjustment of treatment
    labels on y-axis (only considered in 2-D plots); see
    \code{\link[graphics]{text}}.}
  \item{offset.x}{Offset(s) of treatment labels on x-axis (only
    considered in 2-D plots).}
  \item{offset.y}{Offset(s) of treatment labels on y-axis (only
    considered in 2-D plots).}
  \item{pch}{Plot symbol(s) for points; no points printed if equal to \code{NULL}.}
  \item{cex.points}{Magnification(s) to be used for points.}
  \item{col.points}{Colour(s) of points.}
  \item{col.lines}{Line colour.}
  \item{lty.lines}{Line type.}
  \item{lwd.lines}{Line width.}
  \item{arrows}{A logical indicating whether arrows should be printed
    (only considered in 2-D plots).}
  \item{length}{Length of arrows; see \code{\link[graphics]{arrows}}.}  
  \item{grid}{A logical indicating whether grid lines should be added to plot.}
  \item{col.grid}{Colour of grid lines.}
  \item{lty.grid}{Line type of grid lines.}
  \item{lwd.grid}{Line width of grid lines.}
  \item{\dots}{Additional graphical arguments.}
}

\details{
  By default (arguments \code{plottype="scatter"} and
  \code{dim="2d"}), a scatter plot is created showing P-scores (see
  \code{\link{netrank}}) for the first two outcomes considered in the
  generation of a partially ordered set of treatment ranks (using
  \code{\link{netposet}}). In addition to the P-scores, the partially
  order of treatment ranks is shown as lines connecting treatments
  which is analogous to a Hasse diagram. If argument
  \code{dim="3d"}), a 3-D scatter plot is generated showing P-scores
  for the first three outcomes.

  To overcome the restriction of two or three dimension, a biplot
  (Gabriel, 1971) can be generated using argument
  \code{plottype="biplot"}. This is essentially a scatter plot using
  the first two (\code{dim="2d"}) or three (\code{dim="3d"})
  components in a principal components analysis (using
  \code{\link[stats]{prcomp}}).

  Arguments \code{sel.x} and \code{sel.y} can be used to select
  different outcomes to show on x- and y-axis in a 2-D scatter plot;
  argument \code{sel.z} can be used accordingly in a 3-D scatter plot.
  These arguments are ignored for a biplot.

  Note, in order to generate 3-D plots (argument \code{dim="3d"}), R
  package \bold{rgl} is necessary. Note, under macOS the X.Org X
  Window System must be available (see \url{https://www.xquartz.org}).
}

\references{
  Carlsen L, Bruggemann R (2014),
  Partial order methodology: a valuable tool in chemometrics.
  \emph{Journal of Chemometrics},
  \bold{28} 226--34, DOI:10.1002/cem.2569

  Gabriel KR (1971),
  The biplot graphic display of matrices with application to principal
  component analysis.
  \emph{Biometrika},
  \bold{58} 453--67
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}, \code{\link{netrank}}, \code{\link{netposet}}, \code{\link{hasse}}}

\examples{
# Use depression dataset
#
data(Linde2015)
#
# Define order of treatments
#
trts <- c("TCA", "SSRI", "SNRI", "NRI",
          "Low-dose SARI", "NaSSa", "rMAO-A", "Hypericum",
          "Placebo")
#
# Outcome labels
#
outcomes <- c("Early response", "Early remission")
#
# (1) Early response
#
p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
               event = list(resp1, resp2, resp3),
               n = list(n1, n2, n3),
               studlab = id, data = Linde2015, sm = "OR")
#
net1 <- netmeta(p1,
                comb.fixed = FALSE, comb.random = TRUE,
                seq = trts, ref = "Placebo")
#
# (2) Early remission
#
p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
               event = list(remi1, remi2, remi3),
               n = list(n1, n2, n3),
               studlab = id, data = Linde2015, sm = "OR")
#
net2 <- netmeta(p2,
                comb.fixed = FALSE, comb.random = TRUE,
                seq = trts, ref = "Placebo")
#
# Partial order of treatment rankings
#
po <- netposet(netrank(net1, small.values = "bad"),
               netrank(net2, small.values = "bad"),
               outcomes = outcomes)
#
# Scatter plot
#
plot(po)
#
# Biplot
#
plot(po, "biplot")
}

\keyword{Scatterplot of partially ordered rankings}
\keyword{poset}
\keyword{hplot}
