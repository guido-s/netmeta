\name{netconnection}
\alias{netconnection}
\alias{print.netconnection}

\title{Get information on network connectivity (number of subnetworks, distance matrix)}

\description{
  To determine the network structure and to test whether a given
  network is fully connected. Network information is provided as a
  triple of vectors \code{treat1}, \code{treat2}, and \code{studlab}
  where each row corresponds to an existing pairwise treatment
  comparison (\code{treat1}, \code{treat2}) in a study
  (\code{studlab}). The function calculates the number of subnetworks
  (connectivity components; value of 1 corresponds to a fully
  connected network) and the distance matrix (in block-diagonal form
  in the case of subnetworks).
}

\usage{
netconnection(treat1, treat2, studlab, data = NULL, subset = NULL,
              title = "", warn = FALSE)

\method{print}{netconnection}(x, digits = max(4, .Options$digits - 3), ...)
}

\arguments{
  \item{treat1}{Label/Number for first treatment.}
  \item{treat2}{Label/Number for second treatment.}
  \item{studlab}{An optional - but important! - vector with study labels (see Details).}
  \item{data}{An optional data frame containing the study information.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{warn}{A logical indicating whether warnings should be printed.}
  \item{x}{An object of class \code{netconnection}.}
  \item{digits}{Minimal number of significant digits, see
    \code{\link{print.default}}.}
  \item{...}{Additional arguments (ignored at the moment)}
}

\value{
  An object of class \code{netconnection} with corresponding
  \code{print} function. The object is a list containing the following
  components:

  \item{treat1, treat2, studlab, title, warn}{As defined above.}
  \item{k}{Total number of studies.}
  \item{m}{Total number of pairwise comparisons.}
  \item{n}{Total number of treatments.}
  \item{n.subnets}{Number of subnetworks; equal to 1 for a fully connected network.}
  \item{D.matrix}{Distance matrix.}
  \item{A.matrix}{Adjacency matrix.}
  \item{L.matrix}{Laplace matrix.}
  \item{call}{Function call.}
  \item{version}{Version of R package netmeta used to create object.}
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}, \code{\link{netdistance}}}

\examples{
data(Senn2013)

nc1 <- netconnection(treat1, treat2, studlab, data = Senn2013)
nc1

# Extract number of (sub)networks
nc1$n.subnets

# Extract distance matrix
nc1$D.matrix

# Conduct network meta-analysis (results not shown)
net1 <- netmeta(TE, seTE, treat1, treat2, studlab, data = Senn2013)


#
# Artificial example with two subnetworks
#
t1 <- c("G", "B", "B", "D", "A", "F")
t2 <- c("B", "C", "E", "E", "H", "A")
#
nc2 <- netconnection(t1, t2)
nc2

# Number of subnetworks
nc2$n.subnets

# Extract distance matrix
nc2$D.matrix

# Conduct network meta-analysis
# (results in an error message due to unconnected network)
try(net2 <- netmeta(1:6, 1:6, t1, t2, 1:6))

# Conduct network meta-analysis on first subnetwork
net2.1 <- netmeta(1:6, 1:6, t1, t2, 1:6,
                  subset = (t1 \%in\% c("A", "F", "H") & t2 \%in\% c("A", "F", "H")))

# Conduct network meta-analysis on first subnetwork
net2.2 <- netmeta(1:6, 1:6, t1, t2, 1:6,
                  subset = !(t1 \%in\% c("A", "F", "H") & t2 \%in\% c("A", "F", "H")))

summary(net2.1)
summary(net2.2)
}
