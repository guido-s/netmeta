#' Generic function for heat plots
#' 
#' @description
#' Generic function for heat plots
#' 
#' @param x An R object.
#' @param \dots Additional arguments.
#' 
#' @details
#' 
#' For more details, look at the following function to generate heat
#' plots:
#' \itemize{
#' \item \code{\link{heatplot.netmeta}}
#' \item \code{\link[crossnma]{heatplot.crossnma}} (if installed)
#' }
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de
#' }
#' 
#' @keywords hplot
#'
#' @examples
#' # Only consider first five studies (to reduce runtime of example)
#' #
#' studies <- unique(Senn2013$study)
#' Senn2013.5 <- subset(Senn2013, study %in% studies[1:5])
#' 
#' # Transform data from long arm-based to contrast-based format
#' #
#' pw <- pairwise(studlab = study, treat = treatment,
#'   n = n, mean = mean, sd = sd, data = Senn2013.5,
#'   varnames = c("MD", "seMD"))
#' 
#' # Conduct random effects network meta-analysis with
#' # placebo as reference treatment
#' #
#' nma <- netmeta(pw, common = FALSE, reference = "plac")
#' 
#' # Generate a heat plot (with abbreviated treatment labels)
#' #
#' heatplot(nma, nchar.trts = 4) 
#' 
#' @keywords hplot
#' 
#' @rdname heatplot
#' @export heatplot

heatplot <- function(x, ...)
  UseMethod("heatplot")
