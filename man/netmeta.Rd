\name{netmeta}
\alias{netmeta}

\title{Network meta-analysis using graph-theoretical method}

\description{
  TODO: Short description of method.

  Network meta-analysis using graph-theoretical method.
}

\usage{netmeta(TE, seTE, treat1, treat2, studlab, data=NULL, subset=NULL,
         sm="", level=0.95, level.comb=0.95,
         comb.fixed=TRUE, comb.random=FALSE, reference.group="",
         all.treatments=NULL, tau.preset=NULL, title="", warn=TRUE)
}
\arguments{
  \item{TE}{Estimate of treatment effect.}
  \item{seTE}{Standard error of treatment estimate.}
  \item{treat1}{Label/Number for first treatment.}
  \item{treat2}{Label/Number for second treatment.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{sm}{A character string indicating underlying summary measure,
    e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"AS"},
    \code{"MD"}, \code{"SMD"}.}
  \item{level}{The level used to calculate confidence intervals for
    individual comparisons.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
  \item{reference.group}{Reference group.}
  \item{all.treatments}{A logical or value \code{"NULL"}. If
        \code{TRUE}, matrices with all treatment effects, and confidence
        limits will be printed.}
  \item{tau.preset}{Prespecified value for between-study variance tau-squared.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{warn}{A logical indicating whether warnings should be printed
    (e.g., if studies are excluded from meta-analysis due to zero
    standard errors).}
}

\details{TODO: Some / more details on method}

\value{
  TODO: Check and add short description for list elements.

  An object of class \code{c("netmeta")} with corresponding
  \code{print}, \code{summary}, \code{forest} function. The object is
  a list containing the following components:
  \item{TE, seTE, studlab, treat1, treat2}{}
  \item{sm, level, level.comb,}{}
  \item{comb.fixed, comb.random,}{}
  \item{reference.group, all.treatments,}{}
  \item{tau.preset, title, warn}{As defined above.}
  \item{studies}{???}
  \item{narms}{???}
  \item{TE.nma.fixed}{???}
  \item{seTE.nma.fixed}{???}
  \item{lower.nma.fixed}{???}
  \item{upper.nma.fixed}{???}
  \item{leverage.fixed}{???}
  \item{w.fixed}{Weight of individual studies (in fixed effect
    model).}
  \item{TE.fixed}{Matrix with estimated overall treatment effects (fixed effect model).}
  \item{seTE.fixed}{Matrix with standard errors (fixed effect model).}
  \item{lower.fixed, upper.fixed}{Matrices with lower and upper
    confidence interval limits (fixed effect model).}
  \item{zval.fixed, pval.fixed}{Matrices with z-value and p-value for
    test of overall treatment effect (fixed effect model).}
  \item{Q.fixed}{???}
  \item{TE.random}{Matrix with estimated overall treatment effects (random effects model).}
  \item{seTE.random}{Matrix with standard errors (random effects model).}
  \item{lower.random, upper.random}{Matrices with lower and upper
    confidence interval limits (random effects model).}
  \item{zval.random, pval.random}{Matrices with z-value and p-value for
    test of overall treatment effect (random effects model).}
  \item{k}{Total number of studies.}
  \item{m}{Total number of pairwise comparisons.}
  \item{Q}{Overall heterogeneity/inconsistency statistic.}
  \item{df}{Degrees of freedom for test of heterogeneity/inconsistency.}
  \item{pval.Q}{P-value for test of heterogeneity/inconsistency.}
  \item{I2}{I-squared.}
  \item{tau}{Square-root of between-study variance.}
  \item{Q.heterogeneity}{Overall heterogeneity statistic.}
  \item{Q.inconsistency}{Overall inconsistency statistic.}
  \item{A.matrix}{Matrix with ???.}
  \item{L.matrix}{Matrix with ???.}
  \item{Lplus.matrix}{Matrix with ???.}
  \item{Q.matrix}{Matrix with ???.}
  \item{G.matrix}{Matrix with ???.}
  \item{H.matrix}{Matrix with ???.}
  \item{Q.decomp}{Matrix with ???.}
  \item{call}{Function call.}
  \item{version}{Version of R package netmeta used to create object.}
}

\references{
 TODO: Additional references?

  Rücker G (2012),
  Network meta-analysis, electrical networks and graph theory.
  \emph{Research Synthesis Methods},
  \bold{3}, 312--324.
}

\author{Gerta Ruecker \email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metagen}}}

\examples{

## TODO: Better example?

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
