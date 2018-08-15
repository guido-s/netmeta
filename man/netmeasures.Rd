\name{netmeasures}

\alias{netmeasures}

\title{Measures for characterizing a network meta-analysis}

\description{
  This function provides measures for quantifying the direct evidence
  proportion, the mean path length and the minimal parallelism (the
  latter on aggregated and study level) of mixed treatment comparisons          
  (network estimates) as well as the evidence flow per design, see König 
  et al. (2013). These measures support the critical evaluation of the network 
  meta-analysis results by rendering transparent the process of data pooling.
}

\usage{
netmeasures(x,
            random = x$comb.random | !missing(tau.preset),
            tau.preset = x$tau.preset, warn = TRUE)
}

\arguments{
  \item{x}{An object of class \code{netmeta}.}
  \item{random}{A logical indicating whether random effects model
    should be used to calculate network measures.}
  \item{tau.preset}{An optional value for the square-root of the
    between-study variance \eqn{\tau^2}.}
   \item{warn}{A logical indicating whether warnings should be
     printed.}
}

\details{
  The direct evidence proportion gives the absolute contribution of direct      
  effect estimates combined for two-arm and multi-arm studies to one            
  network estimate.

  Concerning indirectness, comparisons with a mean path length beyond
  two should be interpreted with particular caution, as more than two
  direct comparisons have to be combined serially on average.
  
  Large indices of parallelism, either on study-level or on aggregated
  level, can be considered as supporting the validity of a network
  meta-analysis if there is only a small amount of heterogeneity.

  The network estimates for two treatments are linear combinations of
  direct effect estimates comparing these or other treatments. The
  linear coefficients can be seen as the generalization of weights
  known from classical meta-analysis. These coefficients are given in
  the projection matrix \eqn{H} of the underlying model. For multi-arm
  studies, the coefficients depend on the choice of the study-specific
  baseline treatment, but the absolute flow of evidence can be made explicit   
  for each design as shown in König et al. (2013) and is given in \code{H.tilde}.
  
  All measures are calculated based on the fixed effects meta-analysis
  by default. In the case that in function \code{netmeta} the argument
  \code{comb.random=TRUE}, all measures are calculated for a random
  effects model. The value of the square-root of the between-study
  variance \eqn{tau^2} can also be prespecified by argument
  \code{tau.preset} in function \code{netmeta}.
}                                                              

\value{
A list containing the following components:

  \item{random, tau.preset}{As defined above.}
  
  \item{proportion}{A named vector of the direct evidence proportion            
    of each network estimate.}

  \item{meanpath}{A named vector of the mean path length of each network        
  estimate.}

  \item{minpar}{A named vector of the minimal parallelism on aggregated level   
  of each network estimate.}

  \item{minpar.study}{A named vector of the minimal parallelism on              
    study level of each network estimate.}

  \item{H.tilde }{Design-based hat matrix with information on absolute evidence 
    flow per design. The number of rows is equal to the number of
    possible pairwise treatment comparisons and the number of columns
    is equal to the number of designs.}
}

\references{
  König J, Krahn U, Binder H (2013),
  Visualizing the flow of evidence in network meta-analysis and
  characterizing mixed treatment comparisons.
  \emph{Statistics in Medicine},
  \bold{32}(30), 5414--29
}

\author{Ulrike Krahn \email{ulrike.krahn@bayer.com}, Jochem König \email{koenigjo@uni-mainz.de}}

\seealso{\link{netmeta}}

\examples{
data(Senn2013)

#
# Generation of an object of class 'netmeta' with
# reference treatment 'plac', i.e. placebo based 
# on a fixed effects model
#
net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
        data=Senn2013, sm="MD", reference="plac")

#
# Calculate measures based on a fixed effects model
#        
nm1 <- netmeasures(net1)

#
# Plot of minimal parallelism versus mean path length
#
plot(nm1$meanpath, nm1$minpar, pch="",
     xlab="Mean path length", ylab="Minimal parallelism")
text(nm1$meanpath, nm1$minpar, names(nm1$meanpath), cex=0.8)

# Generation of an object of class 'netmeta' with
# reference treatment 'plac' based on a random
# effects model
#
net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data=Senn2013, sm="MD", reference="plac", comb.random=TRUE)

#
# Calculate measures based on a random effects model
#                          
nm2 <- netmeasures(net2)
}                                                                              

\keyword{Network meta-analysis}
