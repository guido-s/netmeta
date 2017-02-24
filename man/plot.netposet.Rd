\name{plot.netposet}

\alias{plot.netposet}

\title{Scatter plot of partially order of treatment ranks}

\description{
  This function generates a scatter plot of a partial order of
  treatment ranks.
}

\usage{
\method{plot}{netposet}(x,
     sel.x = 1, sel.y = 2, sel.z = 3,
     dim = "2d",
     cex = 1.25, col = "black",
     adj.x = 0, adj.y = 1,
     offset.x = 0.005, offset.y = -0.005,
     arrows = TRUE,
     col.lines = "black", lty.lines = 1, lwd.lines = 1,
     length = 0.05,
     grid = TRUE,
     col.grid = "gray", lty.grid = 2, lwd.grid = 1,
     ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta} (mandatory).}
  \item{sel.x}{.}
  \item{sel.y}{.}
  \item{sel.z}{.}
  \item{dim}{A character string indicating whether a 2- or
    3-dimensional plot should be produced, either \code{"2d"} or
    \code{"3d"}.}
  \item{cex}{The magnification to be used for treatment labels.}
  \item{col}{A vector with with colour of treatment labels.}
  \item{adj.x}{Value(s) in [0, 1] to specify adjustment of treatment
    labels on x-axis; see \code{\link[graphics]{text}}.}
  \item{adj.y}{Value(s) in [0, 1] to specify adjustment of treatment
    labels on y-axis; see \code{\link[graphics]{text}}.}
  \item{offset.x}{Offset of treatment labels on x-axis.}
  \item{offset.y}{Offset of treatment labels on x-axis.}
  \item{arrows}{A logical indicating whether arrows should be printed.}
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
...
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}, \code{\link{netposet}}, \code{\link{hasse}}}

\examples{
## Ribassin-Majed et al. (2017), J. Clin. Oncol.
##
outcomes <- c("OS", "PFS", "LC", "DC")
treatments <- c("RT", "IC-RT", "IC-CRT", "CRT",
                "CRT-AC", "RT-AC", "IC-RT-AC")
##
## P-scores (from Table 1)
##
pscore.os  <- c(15, 33, 63, 70, 96, 28, 45) / 100
pscore.pfs <- c( 4, 46, 79, 52, 94, 36, 39) / 100
pscore.lc  <- c( 9, 27, 47, 37, 82, 58, 90) / 100
pscore.dc  <- c(16, 76, 95, 48, 72, 32, 10) / 100

pscore.matrix <- data.frame(pscore.os, pscore.pfs, pscore.lc, pscore.dc)
rownames(pscore.matrix) <- treatments
colnames(pscore.matrix) <- outcomes
pscore.matrix


po <- netposet(pscore.matrix)
po12 <- netposet(pscore.matrix[, 1:2])

oldpar <- par(pty = "s")
plot(po12)
par(oldpar)
}

\keyword{Scatterplot of partially ordered rankings}
\keyword{poset}
\keyword{hplot}
