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
     col = "darkblue", col.lines = 8,
     cex = 1, pch = 3, lty = 1, lwd = 2,
     tcex = 0.8,
     ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta} (mandatory).}
  \item{sel.x}{.}
  \item{sel.y}{.}
  \item{sel.z}{.}
  \item{dim}{.}
  \item{col}{.}
  \item{col.lines}{.}
  \item{cex}{.}
  \item{pch}{.}
  \item{lty}{.}
  \item{lwd}{.}
  \item{tcex}{.}
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
