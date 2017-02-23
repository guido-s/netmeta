\name{netposet}
\alias{netposet}
\alias{print.netposet}

\title{Partial order of treatments in network meta-analysis}

\description{
  Partial order of treatments in network meta-analysis. ... called \emph{poset} ...
}

\usage{
netposet(..., outcomes, small.values)

\method{print}{netposet}(x, ...)
}

\arguments{
  \item{\dots}{See details.}
  \item{outcomes}{A character vector with outcome names.}
  \item{small.values}{A character vector specifying whether small
    treatment effects indicate a beneficial (\code{"good"}) or harmful
    (\code{"bad"}) effect, can be abbreviated.}
  \item{x}{An object of class \code{netposet}.}
}

\details{
 ... called \emph{poset} ...

  In function \code{netposet}, argument \code{\dots} can be either a
  ranking matrix or an arbitrary number of \code{netmeta} or
  \code{netrank} objects. In function \code{print.netposet}, this
  argument is passed on to the printing function.

  Argument \code{outcomes} can be used to label outcomes. If argument
  \code{outcomes} is missing,
  \itemize{
  \item column names of the ranking matrix are used as outcome labels,
  \item capital letters 'A', 'B', \dots are used as outcome labels and
    a corresponding warning is printed.
  }
  
  If argument \code{\dots} consists of \code{netrank} objects,
  argument \code{small.values} can be used to specify for each outcome
  whether small values are good or bad; see \code{\link{netrank}}.
}

\value{
  An object of class \code{netposet} with corresponding \code{print}
  and \code{plot} function. The object is a list containing the
  following components:
  \item{outcomes}{As defined above.}
  \item{treatments}{A character vector with treatment names.}
  \item{M0}{Hasse matrix skipping unnecessary paths.}
  \item{M}{"Full" Hasse matrix.}
  \item{PO}{Matrix with information about partial ordering.}
  \item{mean.rank}{???}
  \item{pscore.matrix}{Ranking matrix.}
  \item{version}{Version of R package netmeta used to create object.}
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}, \code{\link{hasse}}, \code{\link{plot.netposet}}}

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

oldpar <- par(pty = "s")
plot(po12)
par(oldpar)
}

\keyword{Partial order of rankings}
\keyword{poset}
