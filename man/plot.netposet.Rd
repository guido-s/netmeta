\name{plot.netposet}

\alias{plot.netposet}

\title{Scatter plot or bi-plot of partially order of treatment ranks}

\description{
  This function generates a scatter plot or bi-plot of a partial order
  of treatment ranks.
}

\usage{
\method{plot}{netposet}(x,
     pooled=ifelse(x$comb.random, "random", "fixed"),
     sel.x = 1, sel.y = 2, sel.z = 3,
     dim = "2d", plottype = "scatter",
     cex = 1, col = "black",
     adj.x = 0, adj.y = 1,
     offset.x = 0.005, offset.y = -0.005,
     arrows = FALSE,
     col.lines = "black", lty.lines = 1, lwd.lines = 1,
     length = 0.05,
     grid = TRUE,
     col.grid = "gray", lty.grid = 2, lwd.grid = 1,
     ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta} (mandatory).}
  \item{pooled}{A character string indicating whether scatter plot
    should be drawn for fixed effect (\code{"fixed"}) or random
    effects model (\code{"random"}). Can be abbreviated.}
  \item{sel.x}{.}
  \item{sel.y}{.}
  \item{sel.z}{.}
  \item{dim}{A character string indicating whether a 2- or
    3-dimensional plot should be produced, either \code{"2d"} or
    \code{"3d"}.}
  \item{plottype}{A character string indicating whether a scatter plot
    or bi-plot should be produced, either \code{"scatter"} or
    \code{"biplot"}.}
  \item{cex}{The magnification to be used for treatment labels.}
  \item{col}{A vector with with colour of treatment labels.}
  \item{adj.x}{Value(s) in [0, 1] to specify adjustment of treatment
    labels on x-axis (only considered in 2-D plots); see
    \code{\link[graphics]{text}}.}
  \item{adj.y}{Value(s) in [0, 1] to specify adjustment of treatment
    labels on y-axis (only considered in 2-D plots); see
    \code{\link[graphics]{text}}.}
  \item{offset.x}{Offset of treatment labels on x-axis (only
    considered in 2-D plots).}
  \item{offset.y}{Offset of treatment labels on x-axis (only
    considered in 2-D plots).}
  \item{arrows}{A logical indicating whether arrows should be printed
    (only considered in 2-D plots).}
  \item{col.lines}{Line colour.}
  \item{lty.lines}{Line type.}
  \item{lwd.lines}{Line width.}
  \item{length}{Length of arrows; see \code{\link[graphics]{arrows}}.}
  \item{grid}{A logical indicating whether a grid lines should be added to plot.}
  \item{col.grid}{Colour of grid lines.}
  \item{lty.grid}{Line type of grid lines.}
  \item{lwd.grid}{Line width of grid lines.}
  \item{\dots}{Additional graphical arguments.}
}

\details{
  Scatter plot ...

  In order to generate 3-D plots (argument \code{dim = "3d"}), R
  package \bold{rgl} is necessary. Note, under macOS the X.Org X
  Window System must be available (see \url{https://www.xquartz.org}).
}

\references{
  Carlsen L, Bruggemann R (2014),
  Partial order methodology: a valuable tool in chemometrics.
  \emph{Journal of Chemometrics},
  \bold{28} 226--34, DOI:10.1002/cem.2569
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}, \code{\link{netposet}}, \code{\link{hasse}}}

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
# Partial order of treatment rankings (all five outcomes)
#
po <- netposet(netrank(net1, small.values = "bad"),
               netrank(net2, small.values = "bad"),
               outcomes = outcomes)
#
# Scatter plot
#
plot(po)
}

\keyword{Scatterplot of partially ordered rankings}
\keyword{poset}
\keyword{hplot}
