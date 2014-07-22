\name{netheat}

\alias{netheat}

\title{Net heat plot}

\description{
  This function creates a net heat plot, a graphical tool for locating
  inconsistency in network meta-analyses.
}

\usage{
netheat(x, random=FALSE, tau.preset=NULL, ...)                                
}

\arguments{
  \item{x}{An object of class \code{netmeta}.}
  \item{random}{A logical indicating whether the net heat plot should be based  
    on a random effects model.}                                                 
  \item{tau.preset}{An optional value for the square-root of the
    between-study variance \eqn{tau^2} for a random effects model on
    which the net heat plot will be based.}
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
                                                                                
  In the case of \code{random=TRUE}, the net heat plot is based on a random     
  effects model generalised for multivariate meta-analysis in which the         
  between-study variance \eqn{tau^2} is estimated by the method of moments      
  (see Jackson et al., 2012) and embedded in a full design-by-treatment         
  interaction model (see Higgins et al., 2012).                                 
}

\references{
  Krahn U, Binder H, König J (2013),
  A graphical tool for locating inconsistency in network
  meta-analyses.
  \emph{BMC Medical Research Methodology}, \bold{13}, 35.
  
  Jackson D, White IR and Riley RD (2012),                                       
  Quantifying the impact of between-study heterogeneity in multivariate         
  meta-analyses.                                                                
  \emph{Statistics in Medicine}, \bold{31}(29), 3805--3820.                     
                                                                                
  Higgins JPT, Jackson D, Barrett JK, Lu G, Ades AE, White IR (2012),           
  Consistency and inconsistency in network meta-analysis: concepts              
  and models for multi-arm studies.                                             
  \emph{Research Synthesis Methods}, \bold{3}, 98--110.                         
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
## Generate a net heat plot based on a fixed effects model
##
netheat(net1) 

##                                                                              
## Generate a net heat plot based on a random effects model                     
##                                                                              
netheat(net1, random=TRUE)                                                      
}

\keyword{Network meta-analysis}
\keyword{Inconsistency}
