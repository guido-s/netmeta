\name{forest.netbind}

\alias{forest.netbind}

\title{Forest plot for direct and indirect evidence}

\description{
  Forest plot to show direct and indirect evidence in network
  meta-analysis. Furthermore, estimates from network meta-analysis as
  well as prediction intervals can be printed.
}

\usage{
\method{forest}{netbind}(x,
       pooled = ifelse(x$comb.random, "random", "fixed"),
       equal.size = FALSE,
       leftcols = "studlab", leftlabs = "Treatment",
       rightcols = c("effect", "ci"), rightlabs = NULL,
       digits = gs("digits.forest"),
       digits.prop = max(gs("digits.pval") - 2, 2),
       backtransf = x$backtransf,
       lab.NA = "", smlab, ...)
}

\arguments{
  \item{x}{An object of class \code{netbind}.}
  \item{pooled}{A character string indicating whether results for
    the fixed effects (\code{"fixed"}) or random effects model
    (\code{"random"}) should be plotted. Can be abbreviated.}
  \item{equal.size}{A logical indicating whether all squares should be
    of equal size. Otherwise, the square size is proportional to the
    precision of estimates.}
  \item{leftcols}{A character vector specifying columns to be plotted
    on the left side of the forest plot (see Details).}
  \item{leftlabs}{A character vector specifying labels for columns on
    left side of the forest plot.}
  \item{rightcols}{A character vector specifying columns to be plotted
    on the right side of the forest plot (see Details).} 
  \item{rightlabs}{A character vector specifying labels for columns on
    right side of the forest plot.}
  \item{digits}{Minimal number of significant digits for treatment
    effects and confidence intervals, see \code{print.default}.}  
  \item{digits.prop}{Minimal number of significant digits for the
    direct evidence proportion.}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in forest plots. If \code{backtransf=TRUE},
    results for \code{sm="OR"} are presented as odds ratios rather
    than log odds ratios, for example.}
  \item{lab.NA}{A character string to label missing values.}
  \item{smlab}{A label printed at top of figure. By default, text
    indicating either fixed effect or random effects model is
    printed.}
  \item{\dots}{Additional arguments for \code{\link{forest.meta}}
    function.}
}

\details{
  A forest plot, also called confidence interval plot, is drawn in the
  active graphics window.

  The arguments \code{leftcols} and \code{rightcols} can be used to
  specify columns which are plotted on the left and right side of the
  forest plot, respectively. If argument \code{rightcols} is
  \code{FALSE}, no columns will be plotted on the right side.
     
  For more information see help page of \code{\link{forest.meta}}
  function.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{forest.meta}}}

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


\keyword{hplot}
