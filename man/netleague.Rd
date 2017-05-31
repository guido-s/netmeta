\name{netleague}
\alias{netleague}

\title{Print league table for network meta-analysis results}

\description{
  A league table is a square matrix showing all pairwise comparisons
  in a network meta-analysis. Typically, both treatment estimates and
  confidence intervals are shown.
}

\usage{
netleague(x, y,
          comb.fixed = x$comb.fixed, comb.random = x$comb.random,
          seq = x$seq, ci = TRUE, backtransf = TRUE,
          digits = gs("digits"))
}

\arguments{
  \item{x}{An object of class \code{netmeta} (mandatory).}
  \item{y}{An object of class \code{netmeta} (optional).}
  \item{comb.fixed}{A logical indicating whether a league table for
    fixed effect meta-analyses should be printed.}
  \item{comb.random}{A logical indicating whether a league table for
    random effects meta-analyses should be printed.}
  \item{seq}{A character or numerical vector specifying the sequence
    of treatments in rows and columns of a league table.}
  \item{ci}{A logical indicating whether confidence intervals should
    be shown.}
  \item{backtransf}{A logical indicating whether printed results
    should be back transformed. If \code{backtransf=TRUE}, results for
    \code{sm="OR"} are printed as odds ratios rather than log odds
    ratios, for example.}
  \item{digits}{Minimal number of significant digits, see
    \code{print.default}.}
}

\details{
  If argument \code{y} is not provided, the league table contains the
  same information in the lower and upper triangle, i.e., treatment
  comparisons and confidence intervals for network meta-analysis
  object \code{x}.

  If argument \code{y} is provided, the league table contains
  information on treatment comparisons from network meta-analysis
  object \code{x} in the lower triangle and from network meta-analysis
  object \code{y} in the upper triangle.

  R function \code{\link{netrank}} can be used to change the order of
  rows and columns in the league table (see examples).
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}, Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}, \code{\link{netposet}}, \code{\link{netrank}}}

\examples{
# Network meta-analysis of count mortality statistics
#
data(Woods2010)

p0 <- pairwise(treatment, event = r, n = N,
               studlab = author, data = Woods2010, sm = "OR")
net0 <- netmeta(p0)

cilayout(bracket = "(", separator = " - ")
oldopts <- options(width = 100)

# League table for fixed effect model
#
netleague(net0, digits = 2)

# League table for fixed effect and random effects model
#
netleague(net0, comb.random = TRUE, digits = 2)

# Change order of treatments according to treatment ranking
#
netleague(net0, comb.random = TRUE, digits = 2,
          seq = netrank(net0))
#
print(netrank(net0), comb.random = TRUE)


# Use depression dataset
#
data(Linde2015)
cilayout()
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

options(width = 200)
netleague(net1, digits = 2)

netleague(net1, digits = 2, ci = FALSE)
netleague(net2, digits = 2, ci = FALSE)
netleague(net1, net2, digits = 2, ci = FALSE)

netleague(net1, net2, seq = netrank(net1, small = "bad"), ci = FALSE)
netleague(net1, net2, seq = netrank(net2, small = "bad"), ci = FALSE)

print(netrank(net1, small = "bad"), comb.random = TRUE)
print(netrank(net2, small = "bad"), comb.random = TRUE)

options(oldopts)


\dontrun{
# Generate a partial order of treatment rankings 
#
np <- netposet(net1, net2, outcomes = outcomes, small.values = rep("bad",2))
hasse(np)
plot(np)
}
}
