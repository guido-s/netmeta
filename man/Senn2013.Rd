\name{Senn2013}
\docType{data}
\alias{Senn2013}
\title{Network meta-analysis in diabetes}

\description{Network meta-analysis in diabetes comparing effects of a
  number of drugs on the HbA1c value.

  These data are used as an example in Senn et al. (2013) and have
  been preprocessed for use in R package netmeta.
}

\usage{data(Senn2013)}

\format{A data frame with the following columns:
  \describe{
  \item{TE}{Treatment effect}
  \item{seTE}{Standard error of treatment effect}
  \item{treat1}{Treatment 1}
  \item{treat2}{Treatment 2}
  \item{treat1.long}{Treatment 1 (full treatment names)}
  \item{treat2.long}{Treatment 2 (full treatment names)}
  \item{studlab}{Study label}
  }
}

\details{
  Treatment labels provided by columns \code{treat1} and \code{treat2}
  have been abbreviated:
  \itemize{
   \item acar = Acarbose
   \item benf = Benfluorex
   \item metf = Metformin
   \item migl = Miglitol
   \item piog = Pioglitazone
   \item plac = Placebo
   \item rosi = Rosiglitazone
   \item sita = Sitagliptin
   \item sulf = Sulfonylurea
   \item vild = Vildagliptin
  }
  Full treatment names are available in columns \code{treat1.long} and
  \code{treat2.long}.
}

\source{
  Senn S, Gavini F, Magrez D, and Scheen A (2013).
  Issues in performing a network meta-analysis.
  \emph{Statistical Methods in Medical Research},
  \bold{22}(2), 169--189.
  First published online 2012 Jan 3.
}

\seealso{\code{\link{netmeta}}}

\examples{
data(Senn2013)

#
# Fixed effect model (default)
#
net1 <- netmeta(TE, seTE, treat1, treat2,
                studlab, data=Senn2013)
net1
net1$Q.decomp

#
# Forest plot
#
forest(net1, ref="plac")

\dontrun{
#
# Comparison with reference group
#
netmeta(TE, seTE, treat1.long, treat2.long,
        studlab, data=Senn2013, 
        reference="plac")

#
# Random effects model
#
net2 <- netmeta(TE, seTE, treat1, treat2, studlab, data=Senn2013,
                comb.random = TRUE)
net2
forest(net2, ref="plac")
}
}

\keyword{datasets}
