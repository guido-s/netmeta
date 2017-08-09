\name{Woods2010}
\docType{data}
\alias{Woods2010}
\title{Count statistics of survival data}

\description{
Count mortality statistics in randomised controlled trials of
treatments for chronic obstructive pulmonary disease (Woods
et al. (2010), Table 1).
}

\usage{data(Woods2010)}

\format{A data frame with the following columns:
  \describe{
  \item{author}{First author / study name}
  \item{treatment}{Treatment}
  \item{r}{Number of deaths in treatment arm}
  \item{N}{Number of patients in treatment arm}
  }
}

\source{
  Woods BS, Hawkins N, Scott DA (2010).
  Network meta-analysis on the log-hazard scale, combining count and
  hazard ratio statistics accounting for multi-arm trials: A
  tutorial.  
  \emph{BMC Medical Research Methodology}
  \bold{10}, 54.
}

\seealso{\code{\link{pairwise}}, \code{\link{metabin}}, \code{\link{netmeta}}}

\examples{
data(Woods2010)

# Transform data from long arm-based format to contrast-based format
# Argument 'sm' has to be used for odds ratio as summary measure; by
# default the risk ratio is used in the metabin function called
# internally.
p1 <- pairwise(treatment, event = r, n = N,
               studlab = author, data = Woods2010, sm = "OR")
p1

# Conduct network meta-analysis
net1 <- netmeta(p1)
net1

\dontrun{
# Show forest plot
forest(net1, ref = "Placebo", drop = TRUE,
       leftlabs = "Contrast to Placebo")
}
}

\keyword{datasets}
