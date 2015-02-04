\name{pairwise}
\alias{pairwise}

\title{Transform meta-analysis data from arm-based format into
contrast-based format}

\description{This function transforms data that are given in an
arm-based format (e.g. input for WinBUGS is of this format) to a
contrast-based format that is needed as input to R function
\code{\link{netmeta}}. The function can transform data with binary,
continuous, or generic outcomes as well as incidence rates from
arm-based to contrast-based format.
}

\usage{
pairwise(treat, event, n, mean, sd, TE, seTE, time,
         data=NULL, studlab, ...)
}

\arguments{
  \item{treat}{A list with treatment information for individual
    treatment arms.}
  \item{event}{A list with information on number of events for
    individual treatment arms.}
  \item{n}{A list with information on number of observations for
    individual treatment arms.}
  \item{mean}{A list with estimated means for individual treatment
    arms.}
  \item{sd}{A list with information on the standard deviation for
    individual treatment arms.}
  \item{TE}{A list with estimated treatment effects for individual
    treatment arms.}
  \item{seTE}{A list with standard errors of estimated treatment
    effect for individual treatment arms.}
  \item{time}{A list with information on person time at risk for
    individual treatment arms.}
  \item{data}{An optional data frame containing the study
    information.}
  \item{studlab}{A vector with study labels (optional).}
  \item{\dots}{Additional arguments passed-through to the functions to
    calculate effects.}
}

\details{
  R function \code{\link{netmeta}} expects data in a
  \emph{contrast-based format}, where each row corresponds to a
  comparison of two treatments and contains a measure of the treatment
  effect comparing two treatments with standard error, labels for the
  two treatments and an optional study label. Other programs,
  e.g. WinBUGS, require data in an \emph{arm-based format},
  i.e. treatment estimate for each treatment arm instead of a
  difference of the two or more treatments. For example, in
  contrast-based format, a three-arm study contributes three rows with
  treatment comparison and corresponding standard error for comparison
  \emph{A} vs \emph{B}, \emph{A} vs \emph{C}, and \emph{B} vs
  \emph{C}, whereas, in arm-based format, this study contributes one
  row with three treatment estimates and corresponding standard error
  for treatment \emph{A}, \emph{B}, and \emph{C}.

  The pairwise function transforms data given in an arm-based format
  into the contrast-based format which is needed as input to R
  function \code{\link{netmeta}}.

  The pairwise functions expects a number of lists as mandatory input
  depending on the type of data (see next paragrah). List elements can
  be part of a data frame (argument \code{data}, optional). Each list
  entry is a vector with information for each individual study; the
  length of this vector corresponds to the total number of studies
  incorporated in the network meta-analysis. An optional vector with
  study labels can be provided; again this vector can be part of the
  data frame provided with argument \code{data}.

  Pairwise can transform data with binary outcomes (using the
  \code{\link{metabin}} function from R package meta), continuous
  outcomes (\code{\link{metacont}} function), incidence rates
  (\code{\link{metainc}} function), and generic outcomes
  (\code{\link{metagen}} function). The following lists are mandatory
  for these different outcome types:
  \itemize{
   \item treat, event, n (see \code{\link{metabin}})
   \item treat, n, mean, sd (see \code{\link{metacont}})
   \item treat, event, time (see \code{\link{metainc}})
   \item treat, TE, seTE (see \code{\link{metagen}})
  }
  Additional arguments for these functions, e.g. argument \code{sm} to
  specify the summary measure, can be provided using argument
  \code{\dots}.
 
  The \code{treat} list is mandatory to identify the individual
  treatments. Each list consists of as many vectors as the multi-arm
  study with the largest number of treatments. If a single multi-arm
  study has five arms, five vectors have to be provided for each
  lists. Two-arm studies have entries with \code{NA} for the third and
  subsequent vectors.
  
  The value of pairwise is a data frame with as many rows as there are
  pairwise comparisons. For each study with \emph{p} treatments,
  \emph{p*(p-1)/2} contrasts are generated. Each row contains the
  treatment effect, its standard error, the treatments compared and
  the study label.  Further columns are added according to type of
  data.
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
net1 <- netmeta(TE, seTE, treat1, treat2, studlab, data=p1)
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
\dontrun{net2 <- netmeta(TE, seTE, treat1, treat2, studlab, data=p2)
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
net3 <- netmeta(TE, seTE, treat1, treat2, studlab, data=p3)
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
net4 <- netmeta(TE, seTE, treat1, treat2, studlab, data=p4, sm="IRR")
summary(net4)
}
