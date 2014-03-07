\name{netheat}

\alias{netheat}

\title{Net heat plot}

\description{
  This function creates a net heat plot, a graphical tool for locating
  inconsistency in network meta-analyses.
}

\usage{
netheat(x, ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta}.}
  \item{\dots}{Additional arguments.}
}

\details{
  The net heat plot is a matrix visualization proposed by Krahn et
  al. (2013) that highlights hot spots of inconsistency between
  specific direct evidence in the whole network and renders
  transparent possible drivers.

  In this plot, the area of a gray square displays the contribution of
  the direct estimate of one design in the column to a network
  estimate in a row. In combination, the colors show the detailed
  change in inconsistency when relaxing the assumption of consistency
  for the effects of single designs.  The colors on the diagonal
  represent the inconsistency contribution of the corresponding
  design.  The colors on the off-diagonal are associated with the
  change in inconsistency between direct and indirect evidence in a
  network estimate in the row after relaxing the consistency
  assumption for the effect of one design in the column. Cool colors
  indicate an increase and warm colors a decrease: the stronger the
  intensity of the color, the greater the difference between the
  inconsistency before and after the detachment. So, a blue colored
  element indicates that the evidence of the design in the column
  supports the evidence in the row.  A clustering procedure is applied
  to the heat matrix in order to find warm colored hot spots of
  inconsistency. In the case that the colors of a column corresponding
  to design \eqn{d} are identical to the colors on the diagonal, the
  detaching of the effect of design \eqn{d} dissolves the total
  inconsistency in the network.

  The pairwise contrasts corresponding to designs of three- or
  multi-arm studies are marked by '_' following the treatments of the
  design.

  Designs where only one treatment is involved in other designs of the
  network or where the removal of corresponding studies would lead to
  a splitting of the network do not contribute to the inconsistency
  assessment and are not incorporated into the net heat plot.
}

\references{
  Krahn U, Binder H, König J (2013),
  A graphical tool for locating inconsistency in network
  meta-analyses.
  \emph{BMC Medical Research Methodology}, \bold{13}, 35.
}

\author{Ulrike Krahn \email{krahnu@uni-mainz.de}}

\seealso{\link{netmeta}}

\examples{
data(Senn2013)

##
## Generation of an object of class 'netmeta' with
## reference treatment 'plac', i.e. placebo
##
net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
        data=Senn2013, sm="MD", reference="plac")
        
##
## Generate net heat plot 
##
netheat(net1) 
}

\keyword{Network meta-analysis}
\keyword{Inconsistency}
