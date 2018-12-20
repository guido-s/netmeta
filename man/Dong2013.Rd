\name{Dong2013}
\docType{data}
\alias{Dong2013}
\title{Network meta-analysis for chronic obstructive pulmonary disease}

\description{  
  Network meta-analysis comparing inhaled medications in patients with
  chronic obstructive pulmonary disease.
}

\usage{data(Dong2013)}

\format{A data frame with the following columns:
  \describe{
  \item{id}{Study ID}
  \item{treatment}{Treatment}
  \item{death}{Mortality}
  \item{randomized}{Number of individuals in treatment arm}
  }
}

\source{
  Dong Y-H, Lin H-H, Shau W-Y, Wu Y-C, Chang C-H, Lai M-S (2013).  
  Comparative safety of inhaled medications in patients with chronic
  obstructive pulmonary disease: systematic review and mixed treatment
  comparison meta-analysis of randomised controlled trials.
  \emph{Thorax};\bold{68}:48--56.
}

\seealso{\code{\link{pairwise}}, \code{\link{metabin}},
\code{\link{netmetabin}}}

\examples{
data(Dong2013)

# Transform data from long arm-based format to contrast-based format.
# Argument 'sm' has to be used for odds ratio as summary measure; by
# default the risk ratio is used in the metabin function called
# internally.
p1 <- pairwise(treatment, death, randomized, studlab = id,
               data = Dong2013, sm = "OR")

# Conduct Mantel-Haenszel network meta-analysis
netmetabin(p1, ref = "plac")
}

\keyword{datasets}
