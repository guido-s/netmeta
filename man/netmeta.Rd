\name{netmeta}
\alias{netmeta}

\title{Network meta-analysis using graph-theoretical method}

\description{
  Network meta-analysis using graph-theoretical method.
}

\usage{netmeta(TE, seTE, treat1, treat2, studlab, data=NULL,
         sm="", level=0.95, level.comb=0.95,
         comb.fixed=TRUE, comb.random=FALSE, reference.group="",
         all.treatments=NULL, title="", warn=TRUE)
}
\arguments{
  \item{TE}{Estimate of treatment effect.}
  \item{seTE}{Standard error of treatment estimate.}
  \item{treat1}{Label/Number for first treatment.}
  \item{treat2}{Label/Number for second treatment.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information.}
  \item{sm}{A character string indicating underlying summary measure,
    e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"AS"},
    \code{"MD"}, \code{"SMD"}.}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
    reference.group
    all.treatments
  \item{reference.group}{Reference group.}
  \item{all.treatments}{A logical or value \code{"NULL"}. If
        \code{TRUE}, matrices with all treatment effects, and confidence
        limits will be printed.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{warn}{A logical indicating whether warnings should be printed
    (e.g., if studies are excluded from meta-analysis due to zero
    standard errors).}
}

\details{Network meta-analysis using graph-theoretical method.}

\value{...}

\references{
  Rücker G (2012),
  Network meta-analysis, electrical networks and graph theory.
  \emph{Research Synthesis Methods},
  \bold{3}, 312--324.
}

\author{Gerta Ruecker \email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metagen}}}

\examples{
data(example1)

example1$se.effect <- sqrt(1/example1$weight)

## Fixed effect model (default)
##
nm <- netmeta(effect, se.effect, treat1, treat2, author,
              data=example1)
nm

## Random effects model
##
nm <- netmeta(effect, se.effect, treat1, treat2, author,
              data=example1, comb.fixed=FALSE, comb.random=TRUE)
nm

## Fixed effect and random effects model
##
nm <- netmeta(effect, se.effect, treat1, treat2, author,
              data=example1, comb.random=TRUE)
nm


data(thromb2)

nm2 <- netmeta(logor, selogor, treat1, treat2,
               data=thromb2, sm="logOR")
as.data.frame(nm2)


data(thromb)

nm3 <- netmeta(logor, selogor, treat1, treat2,
               data=thromb, sm="logOR")
as.data.frame(nm3)

## Fixed effect and random effects model
##
nm3 <- netmeta(logor, selogor, treat1, treat2,
               data=thromb, sm="logOR", comb.random=TRUE)
nm3
summary(nm3)

}
