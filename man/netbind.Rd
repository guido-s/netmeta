\name{netbind}
\alias{netbind}

\title{Combine network meta-analysis objects}

\description{  
This function can be used to combine network meta-analysis objects which
is especially useful to generate a forest plot with results of several
network meta-analyses.
}

\usage{
netbind(..., name, comb.fixed, comb.random,
        col.study = "black", col.inside = "white",
        col.square = "gray", col.square.lines = col.square,
        backtransf, reference.group, baseline.reference)
}

\arguments{
  \item{...}{Any number of meta-analysis objects (see Details).}
  \item{name}{An optional character vector providing descriptive names
    for the network meta-analysis objects.}
  \item{comb.fixed}{A logical indicating whether results for the fixed
    effects model should be reported.}
  \item{comb.random}{A logical indicating whether results for the random
    effects model should be reported.}
  \item{col.study}{The colour for network estimates and confidence
    limits.}
  \item{col.inside}{The colour for network estimates and confidence
    limits if confidence limits are completely within squares.}
  \item{col.square}{The colour for squares.}
  \item{col.square.lines}{The colour for the outer lines of squares.}
  \item{backtransf}{A logical indicating whether results should be back
    transformed. If \code{backtransf=TRUE} (default), results for
    \code{sm="OR"} are printed as odds ratios rather than log odds
    ratios, for example.}
  \item{reference.group}{Reference group.}
  \item{baseline.reference}{A logical indicating whether results
    should be expressed as comparisons of other treatments versus the
    reference treatment (default) or vice versa. This argument is only
    considered if \code{reference.group} has been specified.}
}

\details{
  This function can be used to combine network meta-analysis objects.
  and is, for
example, useful to generate a forest plot with results of subgroup
analyses.
}

\value{
  An object of class \code{"netbind"} with corresponding \code{forest}
  function. The object is a list containing the following components:
  \item{fixed}{A data frame with results for the fixed effects model.}
  \item{random}{A data frame with results for the random effects model.}
  \item{sm}{Summary measure used in network meta-analyses.}
  \item{level.comb}{Level for confidence intervals.}
  \item{comb.fixed, comb.random, backtransf}{As defined above.}
  \item{reference.group, baseline.reference}{As defined above.}
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}, \code{\link{forest.netbind}}}

\examples{
data(Linde2016)

# Only consider studies including Face-to-face PST (to reduce runtime
# of example)
#
face <- subset(Linde2016, id \%in\% c(16, 24, 49, 118))

# Standard random effects NMA model (with placebo as reference treatment)
#
net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
                data = face, reference.group = "placebo",
                sm = "OR", comb.fixed = FALSE, comb.random = TRUE)

# Additive CNMA model with placebo as inactive component and reference
#
nc1 <- netcomb(net1, inactive = "placebo")

# Combine results of standard NMA and CNMA
#
nb1 <- netbind(nc1, net1,
               name = c("Additive CNMA", "Standard NMA"),
               col.study = c("red", "black"),
               col.square = c("red", "black"))
forest(nb1,
       col.by = "black", addrow.subgroups = FALSE,
       fontsize = 10, spacing = 0.7, squaresize = 0.9,
       label.left = "Favours Placebo",
       label.right = "Favours other")
}
