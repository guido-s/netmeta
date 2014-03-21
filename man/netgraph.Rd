\name{netgraph}

\alias{netgraph}

\title{Network graph}

\description{
  This function generates a graph of the evidence network.
}

\usage{
netgraph(x, seq=x$seq, labels=NULL, plastic=TRUE, thickness=TRUE,
         col="slateblue", lwd=1,
         highlight=NULL, col.highlight="red2", lwd.highlight=2, ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta}.}
  \item{seq}{A character or numerical vector specifying the sequence
    of treatments arrangement anticlockwise.}
  \item{labels}{An optional vector with treatment labels.}
  \item{plastic}{A logical indicating whether the appearance of the
    comparisons should be in 3D look.}
  \item{thickness}{A logical indicating if the thickness of the
    comparisons should be proportional to the inverse standard errors
    or not if \code{plastic=FALSE}.}
  \item{col}{Color for comparisons if \code{plastic=FALSE}.}
  \item{lwd}{A numeric for scaling the line width of comparisons.}
  \item{highlight}{A character identifying a comparison that should be
    marked in the network graph,
    e.g. \code{highlight="treat1:treat2"}.}
  \item{col.highlight}{Color for highlighting the comparison given by
    \code{highlight}.}
  \item{lwd.highlight}{A numeric for the line width for highlighting
    the comparison given by \code{highlight}.}
  \item{\dots}{Additional graphical arguments (ignored at the moment).}
}

\details{
  The nodes in the graph are associated with the treatments and are
  placed on the circumcircle. The comparisons display the observed
  treatment comparisons. The thickness of the comparisons can be
  chosen to be proportional to the inverse standard error of the
  treatment effects, aggregated over all studies including the two
  respective treatments.
}

\references{
  Krahn U, Binder H, König J (2013),
  A graphical tool for locating inconsistency in network
  meta-analyses.
  \emph{BMC Medical Research Methodology}, \bold{13}, 35.
}

\author{Ulrike Krahn \email{krahnu@uni-mainz.de}, Jochem König \email{koenigjo@uni-mainz.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\link{netmeta}}

\examples{
data(Senn2013)

##
## Generation of an object of class 'netmeta' with reference treatment 'plac'
##
net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
        data=Senn2013, sm="MD", reference="plac")

##
## Network graph with specified order of the treatments and one
## highlighted comparison
##
trts <- c("plac", "benf", "migl", "acar", "sulf",
          "metf", "rosi", "piog", "sita", "vild")
netgraph(net1, highlight="acar:metf", seq=trts)
##
## Same network graph using argument 'seq' in netmeta function
##
net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD", reference="plac",
                seq=trts)
netgraph(net2, highlight="acar:metf")


##
## Network graph without 3D look of the comparisons and one highlighted comparison
##
netgraph(net1, plastic=FALSE, highlight="acar:metf")

##
## Network graph without 3D look and comparisons with same thickness
##
netgraph(net1, plastic=FALSE, thickness=FALSE)

##
## Network graph with changed labels and specified order of the
## treatments
##
netgraph(net1, seq=c(1, 3, 5, 2, 9, 4, 7, 6, 8, 10),
         labels=LETTERS[1:10])
}

\keyword{Network meta-analysis}
\keyword{Network graph}
\keyword{hplot}
