\name{hasse}

\alias{hasse}

\title{Hasse diagram}

\description{
  This function generates a Hasse diagram for a partial order of
  treatment ranks in a network meta-analysis.
}

\usage{
hasse(x, newpage = TRUE)
}

\arguments{
  \item{x}{An object of class \code{netposet} (mandatory).}  
  \item{newpage}{A logical value indicating whether a new figure
    should be printed in an existing graphics window. Otherwise, the
    Hasse diagram is added to the existing figure.}
}

\details{
  Generate a Hasse diagram for a partial order of treatment ranks in a
  network meta-analysis (Carlsen and Bruggemann, 2014).
  
  This R function is a wrapper function for R function
  \code{\link[hasseDiagram]{hasse}} in R package \bold{hasseDiagram}
  (Krzysztof Ciomek, \url{https://github.com/kciomek/hasseDiagram}),
  i.e., function \code{hasse} can only be used if R package
  \bold{hasseDiagram} is installed.
}

\references{
  Carlsen L, Bruggemann R (2014),
  Partial order methodology: a valuable tool in chemometrics.
  \emph{Journal of Chemometrics},
  \bold{28} 226--34, DOI:10.1002/cem.2569
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}, \code{\link{netposet}}}

\examples{
# Use depression dataset
#
data(Linde2015)
#
# Define order of treatments
#
trts <- c("TCA", "SSRI", "SNRI", "NRI",
          "Low-dose SARI", "NaSSa", "rMAO-A", "Hypericum",
          "Placebo")
#
# Outcome labels
#
outcomes <- c("Early response", "Early remission")
#
# (1) Early response
#
p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
               event = list(resp1, resp2, resp3),
               n = list(n1, n2, n3),
               studlab = id, data = Linde2015, sm = "OR")
#
net1 <- netmeta(p1,
                comb.fixed = FALSE, comb.random = TRUE,
                seq = trts, ref = "Placebo")
#
# (2) Early remission
#
p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
               event = list(remi1, remi2, remi3),
               n = list(n1, n2, n3),
               studlab = id, data = Linde2015, sm = "OR")
#
net2 <- netmeta(p2,
                comb.fixed = FALSE, comb.random = TRUE,
                seq = trts, ref = "Placebo")
#
# Partial order of treatment rankings (all five outcomes)
#
po <- netposet(netrank(net1, small.values = "bad"),
               netrank(net2, small.values = "bad"),
               outcomes = outcomes)
#
# Hasse diagram
#
hasse(po)
}

\keyword{Hasse diagram}
\keyword{hplot}
