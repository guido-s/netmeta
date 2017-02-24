\name{netposet}
\alias{netposet}
\alias{print.netposet}

\title{Partial order of treatments in network meta-analysis}

\description{
  Partial order of treatments in network meta-analysis. The set of
  treatments in a network is called a partially ordered set (in short,
  a \emph{poset}), if different outcomes provide different treatment
  ranking lists.
}

\usage{
netposet(..., outcomes, treatments, small.values)

\method{print}{netposet}(x, ...)
}

\arguments{
  \item{\dots}{See details.}
  \item{outcomes}{A character vector with outcome names.}
  \item{treatments}{A character vector with treatment names.}
  \item{small.values}{See details.}
  \item{x}{An object of class \code{netposet}.}
}

\details{
  In network meta-analysis, frequently different outcomes are
  considered which may each provide a different ordering of
  treatments. The concept of a partially ordered set (in short, a
  \emph{poset}) of treatments can be used to gain further insights in
  situations with apparently conflicting orderings.

  In function \code{netposet}, argument \code{\dots} can be any of the
  following:
  \itemize{
  \item single ranking matrix with each column providing P-scores
    (Rücker and Schwarzer 2015) or SUCRA values (Salanti et al. 2011)
    for an outcome and rows corresponding to treatments;
  \item arbitrary number of \code{netrank} objects providing P-scores;
  \item arbitrary number of \code{netmeta} objects.
  }
  Note, albeit in general a ranking matrix is not constrained to have
  values between 0 and 1, \code{netposet} stops with an error in this
  case as this function expects a matrix with P-scores or SUCRA
  values.

  Argument \code{outcomes} can be used to label outcomes. If argument
  \code{outcomes} is missing,
  \itemize{
  \item column names of the ranking matrix are used as outcome labels
    (if first argument is a ranking matrix and column names are
    available);
  \item capital letters 'A', 'B', \dots are used as outcome labels and
    a corresponding warning is printed.
  }
  Argument \code{treatments} can be used to provide treatment labels
  if the first argument is a ranking matrix. If argument
  \code{treatment} is missing,
  \itemize{
  \item row names of the ranking matrix are used as treatment labels
    (if available);
  \item letters 'a', 'b', \dots are used as treatment labels and a
    corresponding warning is printed.
  }
  If argument \code{\dots} consists of \code{netmeta} objects,
  \code{netrank} is called internally to calculate P-scores. In this
  case, argument \code{small.values} can be used to specify for each
  outcome whether small values are good or bad; see
  \code{\link{netrank}}. This argument is ignored for a ranking matrix
  and \code{netrank} objects.
  
  In function \code{print.netposet}, argument \code{\dots} is passed on
  to the printing function.
}

\value{
  An object of class \code{netposet} with corresponding \code{print},
  \code{plot}, and \code{hasse} function. The object is a list
  containing the following components:
  \item{P.matrix}{Ranking matrix with rows corresponding to treatments
    and columns corresponding to outcomes.}
  \item{M0.matrix}{Hasse matrix skipping unnecessary paths.}
  \item{M.matrix}{"Full" Hasse matrix.}
  \item{O.matrix}{Matrix with information about partial ordering.}
  \item{version}{Version of R package netmeta used to create object.}
}

\references{
  Rücker G & Schwarzer G (2015),
  Ranking treatments in frequentist network meta-analysis works
  without resampling methods.
  \emph{BMC Medical Research Methodology},
  \bold{15}, 58, DOI:10.1186/s12874-015-0060-8 .

  Salanti G, Ades AE, Ioannidis JP (2011). 
  Graphical methods and numerical summaries for presenting results
  from multiple-treatment meta-analysis: an overview and tutorial.
  \emph{Journal of Clinical Epidemiology},
  \bold{64}(2), 163--171.
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
