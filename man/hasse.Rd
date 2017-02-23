\name{hasse}

\alias{hasse}

\title{Hasse diagram}

\description{
  This function generates a Hasse diagram for a partial order of
  treatment ranks in a network meta-analysis.
}

\usage{
hasse(x, newpage = TRUE)
}

\arguments{
  \item{x}{An object of class \code{netposet} (mandatory).}
  \item{newpage}{A logical value indicating whether a new figure should
    be printed in an existing graphics window.}
}

\details{
  This R function is a wrapper function for R function
  \code{\link[hasseDiagram]{hasse}} in R package \bold{hasseDiagram}
  (Krzysztof Ciomek, \url{https://github.com/kciomek/hasseDiagram}),
  i.e., function \code{hasse} can only be used if R package
  \bold{hasseDiagram} is installed.
}

\references{
  ...
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}, \code{\link{netposet}}}

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
po
po12

hasse(po)
hasse(po12)
}

\keyword{Hasse diagram}
\keyword{hplot}
