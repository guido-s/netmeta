\name{netmetabin}
\alias{netmetabin}

\title{Network meta-analysis of binary outcome data}

\description{
  Network meta-analysis of binary outcome data using the
  Mantel-Haenszel method or the non-central hypergeometric
  distribution distribution.
}

\usage{
netmetabin(event1, n1, event2, n2,
           treat1, treat2, studlab,
           data = NULL, subset = NULL,
           sm,
           method = "MH",
           incr = if (method == "Inverse") 0.5 else 0,
           allincr = FALSE, addincr = FALSE,
           allstudies = FALSE,
           level = 0.95, level.comb = 0.95,
           comb.fixed = gs("comb.fixed"),
           comb.random = method == "Inverse" &
             (gs("comb.random") | !is.null(tau.preset)),
           ##
           prediction = FALSE,
           level.predict = 0.95,
           ##
           reference.group = "",
           baseline.reference = TRUE,
           all.treatments = NULL,
           seq = NULL,
           ##
           tau.preset = NULL,
           ##
           tol.multiarm = 0.0005,
           details.chkmultiarm = FALSE,
           ##
           sep.trts = ":",
           nchar.trts = 666,
           ##
           backtransf = gs("backtransf"),
           ##
           title = "",
           keepdata = gs("keepdata"),
           warn = TRUE,
           ...)
}

\arguments{
  \item{event1}{...}
  \item{n1}{...}
  \item{event2}{...}
  \item{n2}{...}
  \item{treat1}{...}
  \item{treat2}{...}
  \item{studlab}{...}
  \item{data}{...}
  \item{subset}{...}
  \item{sm}{...}
  \item{method}{...}
  \item{incr}{...}
  \item{allincr}{...}
  \item{addincr}{...}
  \item{allstudies}{...}
  \item{level}{...}
  \item{level.comb}{...}
  \item{comb.fixed}{...}
  \item{comb.random}{...}
  \item{prediction}{...}
  \item{level.predict}{...}
  \item{reference.group}{...}
  \item{baseline.reference}{...}
  \item{all.treatments}{...}
  \item{seq}{...}
  \item{tau.preset}{...}
  \item{tol.multiarm}{...}
  \item{details.chkmultiarm}{...}
  \item{sep.trts}{...}
  \item{nchar.trts}{...}
  \item{backtransf}{...}
  \item{title}{...}
  \item{keepdata}{...}
  \item{warn}{...}
  \item{\dots}{...}
}

\details{
  ... more details ...

  Comparisons belonging to multi-arm studies are identified by
  identical study labels (argument \code{studlab}). It is therefore
  important to use identical study labels for all comparisons
  belonging to the same multi-arm study, e.g., study label
  "Willms1999" for the three-arm study in Senn et al. (2013).

  Data entry for this function is in \emph{contrast-based} format,
  that is, each line corresponds to a single pairwise comparison
  between two treatments (arguments \code{treat1}, \code{treat2},
  \code{event1}, \code{n1}, \code{event2}, and \code{n2}). If data are
  provided in \emph{arm-based} format, that is, number of events and
  participants are given for each treatment arm separately, function
  \code{\link{pairwise}} can be used to transform the data to
  \emph{contrast-based} format (see help page of function
  \code{\link{pairwise}}).
  
  Note, all pairwise comparisons must be provided for a multi-arm
  study. Consider a multi-arm study of \emph{p} treatments with known
  variances. For this study, treatment effects and standard errors
  must be provided for each of \emph{p}(\emph{p} - 1)/2 possible
  comparisons. For instance, a three-arm study contributes three
  pairwise comparisons, a four-arm study even six pairwise
  comparisons. Function \code{\link{pairwise}} automatically
  calculates all pairwise comparisons for multi-arm studies.
  
  For \code{method = "Inverse"}, both fixed effects and random effects
  models are calculated regardless of values choosen for arguments
  \code{comb.fixed} and \code{comb.random}. Accordingly, the network
  estimates for the random effects model can be extracted from
  component \code{TE.random} of an object of class \code{"netmeta"}
  even if argument \code{comb.random=FALSE}. However, all functions in
  R package \bold{netmeta} will adequately consider the values for
  \code{comb.fixed} and \code{comb.random}. E.g. function
  \code{\link{print.summary.netmeta}} will not print results for the
  random effects model if \code{comb.random=FALSE}. By default
  (\code{comb.fixed=TRUE}), a fixed effects network meta-analysis is
  reported. Use argument \code{comb.random=TRUE} to report results of
  a random effects model.

  For \code{method = "MH"} and \code{method = "NCH"}, only a fixed
  effects model is available.

  By default, treatment names are not abbreviated in
  printouts. However, in order to get more concise printouts, argument
  \code{nchar.trts} can be used to define the minimum number of
  characters for abbreviated treatment names. R function
  \code{\link{treats}} is utilised internally to create abbreviated
  treatment names.

  Names of treatment comparisons are created by concatenating
  treatment labels of pairwise comparisons using \code{sep.trts} as
  separator (see \code{\link{paste}}). These comparison names are used
  in the covariance matrices \code{Cov.fixed} and \code{Cov.random}
  and in some R functions, e.g, \code{\link{decomp.design}}. By
  default, a colon is used as the separator. If any treatment label
  contains a colon the following characters are used as separator (in
  consecutive order): \code{"-"}, \code{"_"}, \code{"/"}, \code{"+"},
  \code{"."}, \code{"|"}, and \code{"*"}.  If all of these characters
  are used in treatment labels, a corresponding error message is
  printed asking the user to specify a different separator.
}

\value{
  An object of class \code{netmetabin} and \code{netmeta} with
  corresponding \code{print}, \code{summary}, \code{forest}, and
  \code{netrank} function. The object is a list containing the
  following components:
  \item{studlab, treat1, treat2}{As defined above.}
  \item{n1, n2, event1, event2}{As defined above.}
  \item{k}{Total number of studies.}
  \item{m}{Total number of pairwise comparisons.}
  \item{n}{Total number of treatments.}
  \item{d}{Total number of designs.}
  \item{trts}{Treatments included in network meta-analysis.}
  \item{studies}{Study labels coerced into a factor with its levels
    sorted alphabetically.}
  \item{TE.fixed, TE.random}{\emph{n}x\emph{n} matrix with estimated overall
    treatment effects (fixed effect / random effects model).}
  \item{seTE.fixed, seTE.random}{\emph{n}x\emph{n} matrix with standard errors
    (fixed effect / random effects model).}
  \item{lower.fixed, upper.fixed, lower.random, upper.random}{\emph{n}x\emph{n}
    matrices with lower and upper confidence interval limits (fixed
    effect / random effects model).}
  \item{zval.fixed, pval.fixed, zval.random, pval.random}{\emph{n}x\emph{n}
    matrices with z-value and p-value for test of overall treatment
    effect (fixed effect / random effects model).}
  \item{TE.direct.fixed, TE.direct.random}{\emph{n}x\emph{n} matrix
    with estimated treatment effects from direct evidence (fixed
    effect / random effects model).}
  \item{seTE.direct.fixed, seTE.direct.random}{\emph{n}x\emph{n}
    matrix with estimated standard errors from direct evidence (fixed
    effect / random effects model).}
  \item{lower.direct.fixed, upper.direct.fixed, lower.direct.random,
    upper.direct.random}{\emph{n}x\emph{n} matrices with lower and
    upper confidence interval limits from direct evidence (fixed
    effect / random effects model).}
  \item{zval.direct.fixed, pval.direct.fixed, zval.direct.random,
    pval.direct.random}{\emph{n}x\emph{n} matrices with z-value and
    p-value for test of overall treatment effect from direct evidence
    (fixed effect / random effects model).}
  \item{Q}{Overall heterogeneity / inconsistency statistic.}
  \item{df.Q}{Degrees of freedom for test of heterogeneity / inconsistency.}
  \item{pval.Q}{P-value for test of heterogeneity / inconsistency.}
  \item{y}{??? REMOVE or DESCRIBE ???}
  \item{V}{??? REMOVE or DESCRIBE ???}
  \item{H}{??? REMOVE or DESCRIBE ???}
  \item{X}{??? REMOVE or DESCRIBE ???}
  \item{A.matrix}{Adjacency matrix (\emph{n}x\emph{n}).}
  \item{data}{Data (contrast-based format).}
  \item{data.long}{Data (arm-based format).}
  \item{sm, level, level.comb}{As defined above.}
  \item{comb.fixed, comb.random}{As defined above.}
  \item{prediction, level.predict}{As defined above.}
  \item{reference.group, baseline.reference, all.treatments}{As
    defined above.}
  \item{seq, tau.preset, tol.multiarm, details.chkmultiarm}{As defined
    above.}
  \item{sep.trts, nchar.trts}{As defined above.}
  \item{backtransf, title, warn}{As defined above.}
  \item{call}{Function call.}
  \item{version}{Version of R package netmeta used to create object.}
}

\references{
  Senn S, Gavini F, Magrez D, and Scheen A (2013),
  Issues in performing a network meta-analysis.
  \emph{Statistical Methods in Medical Research},
  \bold{22}(2), 169--189.
  First published online 2012 Jan 3.
}

\author{Orestis Efthimiou \email{oremiou@gmail.com}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{pairwise}}, \code{\link{netmeta}}}

\examples{
# To be added ...
}
