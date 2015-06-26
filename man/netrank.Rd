\name{netrank}
\alias{netrank}
\alias{print.netrank}

\title{Frequentist method to rank treatments in network}

\description{
  Ranking treatments in frequentist network meta-analysis without resampling methods.
}

\usage{
netrank(x, small.values="good")

\method{print}{netrank}(x, sort=TRUE, digits=max(4, .Options$digits - 3), ...)            
}

\arguments{
  \item{x}{An object of class \code{netmeta} (netrank function) or
    \code{netrank} (print function).}
  \item{small.values}{A character string specifying whether small
    treatment effects indicate a beneficial (\code{"good"}) or harmful
    (\code{"bad"}) effect, can be abbreviated.}
  \item{sort}{A logical indicating whether printout should be sorted
    by decreasing P-score.}
  \item{digits}{Minimal number of significant digits, see
    \code{\link{print.default}}.}    
  \item{\dots}{Additional arguments passed on to
    \code{\link{print.data.frame}} function (used internally).}
}

\details{
  Treatments are ranked based on a network meta-analysis. Ranking is
  performed by P-scores. P-scores are based solely on the point
  estimates and standard errors of the network estimates. They measure
  the extent of certainty that a treatment is better than another
  treatment, averaged over all competing treatments (Rücker and
  Schwarzer 2015).

  The P-score of treatment \emph{i} is defined as the mean of all 1 -
  P[\emph{j}] where P[\emph{j}] denotes the one-sided P-value of
  accepting the alternative hypothesis that treatment \emph{i} is
  better than one of the competing treatments \emph{j}. Thus, if
  treatment \emph{i} is better than many other treatments, many of
  these P-values will be small and the P-score will be large. Vice
  versa, if treatment \emph{i} is worse than most other treatments,
  the P-score is small.

  The P-score of treatment \emph{i} can be interpreted as the mean
  extent of certainty that treatment \emph{i} is better than another
  treatment. This interpretation is comparable to that of the Surface
  Under the Cumulative RAnking curve (SUCRA) which is the rank of
  treatment \emph{i} within the range of treatments, measured on a
  scale from 0 (worst) to 1 (best) (Salanti et al. 2011).
}

\value{
  An object of class \code{c("netrank")} with corresponding
  \code{print} function. The object is a list containing the following
  components:
  \item{Pscore}{A named numeric vector with P-scores.}
  \item{Pmatrix}{Numeric matrix based on pairwise one-sided p-values.}
  \item{small.values, x}{As defined above.}
  \item{version}{Version of R package netmeta used to create object.}
}

\references{
  Rücker G & Schwarzer G (2015),
  Ranking treatments in frequentist network meta-analysis works
  without resampling methods.
  Submitted.

  Salanti G, Ades AE, Ioannidis JP (2011). 
  Graphical methods and numerical summaries for presenting results
  from multiple-treatment meta-analysis: an overview and tutorial.
  \emph{Journal of Clinical Epidemiology},
  \bold{64}(2), 163--171.
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}}

\examples{

data(Senn2013)

net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD")

nr1 <- netrank(net1)

nr1

print(nr1, sort=FALSE)
}

\keyword{print}
