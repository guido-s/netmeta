\name{pairwise}
\alias{pairwise}

\title{Transform meta-analysis data from two arm-based formats into
contrast-based format}

\description{This function transforms data that are given in wide or
long arm-based format (e.g. input format for WinBUGS) to a
contrast-based format that is needed as input to R function
\code{\link{netmeta}}. The function can transform data with binary,
continuous, or generic outcomes as well as incidence rates from
arm-based to contrast-based format.
}

\usage{
pairwise(treat, event, n, mean, sd, TE, seTE, time,
         data=NULL, studlab,
         incr=0.5, allincr=FALSE, addincr=FALSE, allstudies=FALSE,
         keep=NULL, func="mode", ties.method="random",
         ...)
}

\arguments{
  \item{treat}{A list or vector with treatment information for
    individual treatment arms (see Details).}
  \item{event}{A list or vector with information on number of events
    for individual treatment arms (see Details).}
  \item{n}{A list or vector with information on number of observations
    for individual treatment arms (see Details).}
  \item{mean}{A list or vector with estimated means for individual
    treatment arms (see Details).}
  \item{sd}{A list or vector with information on the standard
    deviation for individual treatment arms (see Details).}
  \item{TE}{A list or vector with estimated treatment effects for
    individual treatment arms (see Details).}
  \item{seTE}{A list or vector with standard errors of estimated
    treatment effect for individual treatment arms (see Details).}
  \item{time}{A list or vector with information on person time at risk
    for individual treatment arms (see Details).}
  \item{data}{An optional data frame containing the study
    information.}
  \item{studlab}{A vector with study labels (optional).}
  \item{incr}{A numerical value which is added to each cell frequency
    for studies with a zero cell count.}
  \item{allincr}{A logical indicating if \code{incr} is added to each
    cell frequency of all studies if at least one study has a zero
    cell count. If FALSE (default), \code{incr} is added only to each
    cell frequency of studies with a zero cell count.}
  \item{addincr}{A logical indicating if \code{incr} is added to each
    cell frequency of all studies irrespective of zero cell counts.}
  \item{allstudies}{A logical indicating if studies with zero or all
    events in two treatment arms are to be included in the
    meta-analysis (applies only if \code{sm} is equal to \code{"RR"}
    or \code{"OR"}).}
  \item{keep}{A vector with names of variables to keep in dataset
    (only relevant for data in long arm-based format, see Details).}
  \item{func}{A vector with function names to summarize values within
    individual studies; see Details.}
  \item{ties.method}{A character string describing how ties are
    handled if \code{func="mode"}; see Details.}
  \item{\dots}{Additional arguments passed-through to the functions to
    calculate effects.}
}

\details{
  R function \code{\link{netmeta}} expects data in a
  \bold{contrast-based format}, where each row corresponds to a
  comparison of two treatments and contains a measure of the treatment
  effect comparing two treatments with standard error, labels for the
  two treatments and an optional study label. In contrast-based
  format, a three-arm study contributes three rows with treatment
  comparison and corresponding standard error for pairwise comparison
  \emph{A} vs \emph{B}, \emph{A} vs \emph{C}, and \emph{B} vs \emph{C}
  whereas a four-arm study contributes six rows / pairwise
  comparisons: \emph{A} vs \emph{B}, \emph{A} vs \emph{C}, \dots,
  \emph{C} vs \emph{D}.

  Other programs for network meta-analysis in WinBUGS and Stata
  require data in an \emph{arm-based} format, i.e. treatment estimate
  for each treatment arm instead of a difference of two treatments. A
  common \bold{(wide) arm-based format} consists of one data row per
  study, containing treatment and other necessary information for all
  study arms. For example, a four-arm study contributes one row with
  four treatment estimates and corresponding standard errors for
  treatments \emph{A}, \emph{B}, \emph{C}, and \emph{D}. Another
  possible arm-based format is a long format where each row
  corresponds to a single study arm. Accordingly, in the \bold{long
  arm-based format} a study contributes as many rows as treatments
  considered in the study.
  
  The pairwise function transforms data given in (wide or long)
  arm-based format into the contrast-based format which consists of
  \emph{pairwise} comparisons and is needed as input to the
  \code{\link{netmeta}} function.

  The pairwise function can transform data with binary outcomes (using
  the \code{\link{metabin}} function from R package meta), continuous
  outcomes (\code{\link{metacont}} function), incidence rates
  (\code{\link{metainc}} function), and generic outcomes
  (\code{\link{metagen}} function). Depending on the outcome, the
  following arguments are mandatory:  
  \itemize{
   \item treat, event, n (see \code{\link{metabin}})
   \item treat, n, mean, sd (see \code{\link{metacont}})
   \item treat, event, time (see \code{\link{metainc}})
   \item treat, TE, seTE (see \code{\link{metagen}})
  }
  
  Argument \code{treat} is mandatory to identify the individual
  treatments. The other arguments contain outcome specific data. These
  arguments must be either lists (wide arm-based format, i.e., one row
  per study) or vectors (long arm-based format, i.e. multiple rows per
  study) of the same length.

  For the wide arm-based format, each list consists of as many vectors
  of the same length as the multi-arm study with the largest number of
  treatments. If a single multi-arm study has five arms, five vectors
  have to be provided for each lists. Two-arm studies have entries
  with \code{NA} for the third and subsequent vectors. Each list entry
  is a vector with information for each individual study; i.e. the
  length of this vector corresponds to the total number of studies
  incorporated in the network meta-analysis. Typically, list elements
  are part of a data frame (argument \code{data}, optional); see
  Examples. An optional vector with study labels can be provided which
  can be part of the data frame.

  In the long arm-based format, argument \code{studlab} is mandatory
  to identify rows contributing to individual studies.

  Additional arguments for meta-analysis functions can be provided
  using argument \code{'\dots'}. The following is a list of some
  important arguments:
  \tabular{cll}{
   Argument \tab Description \tab R function \cr
   \code{sm} \tab Summary measure \tab \code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metainc}}, \code{\link{metagen}} \cr
   \code{method} \tab Meta-analysis method \tab \code{\link{metabin}}, \code{\link{metainc}} \cr
   \code{method.tau} \tab Estimation of between-study variance \tab \code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metainc}}, \code{\link{metagen}} \cr
   \code{method.smd} \tab Standardised mean difference \tab \code{\link{metacont}}
   }
  More information on these as well as other arguments is given in the
  help pages of R functions \code{\link{metabin}},
  \code{\link{metacont}}, \code{\link{metainc}}, and
  \code{\link{metagen}}, respectively.
    
  The value of pairwise is a data frame with as many rows as there are
  pairwise comparisons. For each study with \emph{p} treatments,
  \emph{p*(p-1)/2} contrasts are generated. Each row contains the
  treatment effect (\code{TE}), its standard error (\code{seTE}), the
  treatments compared ((\code{treat1}), (\code{treat2})) and the study
  label ((\code{studlab})). Further columns are added according to
  type of data.

  If data are provided in the wide arm-based format, all variables
  from the original dataset are also part of the output dataset. If
  data are provided in the long arm-based format, no variable will be
  added to the output dataset by default. However, using argument
  \code{keep}, a vector with variable names can be provided that
  should be added to the output dataset. In this case, for each study,
  unique values will be calculated for the additional variable(s)
  based on the argument \code{func} (which must be of length 1 or the
  same length as argument \code{keep}): "mode" (most common value),
  "min" (minimum value), "max", "mean", "median", and "sum". In order
  to determine the most common value, the argument \code{ties.method}
  can be used in the case of ties with "first" meaning that the first
  / smallest value will be selected; similar for "last" (last /
  largest value) and "random" (random selection).
  
  Any variable names in the original dataset that are identical to the
  main variables (i.e., "TE", "seTE", ...) will be renamed to variable
  name with ending ".orig".
}

\value{
  A data frame with the following columns
  \item{TE}{Treatment estimate comparing treatment 'treat1' and 'treat2'.}
  \item{seTE}{Standard error of treatment estimate.}
  \item{studlab}{Study labels.}
  \item{treat1}{First treatment in comparison.}
  \item{treat2}{Second treatment in comparison.}
  \item{event1}{Number of events for first treatment arm (for metabin and metainc).}
  \item{event2}{Number of events for second treatment arm (for metabin and metainc).}
  \item{n1}{Number of observations for first treatment arm (for metabin and metacont).}
  \item{n2}{Number of observations for second treatment arm (for metabin and metacont).}
  \item{mean1}{Estimated mean for first treatment arm (for metacont).}
  \item{mean2}{Estimated mean for second treatment arm (for metacont).}
  \item{sd1}{Standard deviation for first treatment arm (for metacont).}
  \item{sd2}{Standard deviation for second treatment arm (for metacont).}
  \item{TE1}{Estimated treatment effect for first treatment arm (for metagen).}
  \item{TE2}{Estimated treatment effect for second treatment arm (for metagen).}
  \item{seTE1}{Standard error of estimated treatment effect for first treatment arm (for metagen).}
  \item{seTE2}{Standard error of estimated treatment effect for second treatment arm (for metagen).}
  \item{time1}{Person time at risk for first treatment arm (for metainc).}
  \item{time2}{Person time at risk for second treatment arm (for metainc).}

  Additional columns are available in the data frame depending on the
  specific arm-based format and argument \code{keep}; see Details.
}

\author{Gerta Rücker\email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}, \code{\link{metacont}},
\code{\link{metagen}}, \code{\link{metabin}}, \code{\link{metainc}},
\code{\link{netgraph}}}

\examples{
#
# Example using continuous outcomes (internal call of function metacont)
#
data(parkinson)
# Transform data from arm-based format to contrast-based format
p1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
               n=list(n1, n2, n3),
               mean=list(y1, y2, y3),
               sd=list(sd1, sd2, sd3),
               data=parkinson, studlab=Study)
p1

# Conduct network meta-analysis
net1 <- netmeta(p1)
net1

# Draw network graphs
netgraph(net1, points=TRUE, cex.points=3, cex=1.5,
         thickness="se.fixed")
netgraph(net1, points=TRUE, cex.points=3, cex = 1.5,
         plastic=TRUE, thickness="se.fixed",
         iterate=TRUE)
netgraph(net1, points=TRUE, cex.points=3, cex = 1.5,
         plastic=TRUE, thickness="se.fixed",
         iterate=TRUE, start="eigen")


#
# Example using generic outcomes (internal call of function metagen)
#
# Calculate standard error for means y1, y2, y3
parkinson$se1 <- with(parkinson, sqrt(sd1^2/n1))
parkinson$se2 <- with(parkinson, sqrt(sd2^2/n2))
parkinson$se3 <- with(parkinson, sqrt(sd3^2/n3))
# Transform data from arm-based format to contrast-based format using
# means and standard errors (note, argument 'sm' has to be used to
# specify that argument 'TE' is a mean difference)
p2 <- pairwise(list(Treatment1, Treatment2, Treatment3),
               TE=list(y1, y2, y3),
               seTE=list(se1, se2, se3),
               data=parkinson, studlab=Study,
               sm="MD")
p2

# Compare pairwise objects p1 (based on continuous outcomes) and p2
# (based on generic outcomes)
all.equal(p1[, c("TE", "seTE", "studlab", "treat1", "treat2")],
          p2[, c("TE", "seTE", "studlab", "treat1", "treat2")])

# Same result as network meta-analysis based on continuous outcomes
# (object net1)
\dontrun{net2 <- netmeta(p2)
net2}


#
# Example with binary data
#
data(smokingcessation)
# Transform data from arm-based format to contrast-based format
# (interal call of metabin function). Argument 'sm' has to be used for
# odds ratio as risk ratio (sm="RR") is default of metabin function.
p3 <- pairwise(list(treat1, treat2, treat3),
               list(event1, event2, event3),
               list(n1, n2, n3),
               data=smokingcessation,
               sm="OR")
p3

# Conduct network meta-analysis
net3 <- netmeta(p3)
net3

#
# Example with incidence rates
#
data(dietaryfat)

# Transform data from arm-based format to contrast-based format
p4 <- pairwise(list(treat1, treat2, treat3),
               list(d1, d2, d3),
               time=list(years1, years2, years3),
               studlab=ID,
               data=dietaryfat)
p4

# Conduct network meta-analysis using incidence rate ratios (sm="IRR").
# Note, the argument 'sm' is not necessary as this is the default in R
# function metainc called internally
net4 <- netmeta(p4, sm="IRR")
summary(net4)

#
# Example with long data format
#
data(Woods2010)

# Transform data from long arm-based format to contrast-based format
# Argument 'sm' has to be used for odds ratio as summary measure; by
# default the risk ratio is used in the metabin function called
# internally.
p5 <- pairwise(treatment, event=r, n=N,
               studlab=author, data=Woods2010, sm="OR")
p5

# Conduct network meta-analysis
net5 <- netmeta(p5)
net5
}
