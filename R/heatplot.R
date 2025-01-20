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
#' data(Senn2013)
#' 
#' # Only consider first five studies (to reduce runtime of example)
#' #
#' studies <- unique(Senn2013$studlab)
#' Senn2013.5 <- subset(Senn2013, studlab %in% studies[1:5])
#' 
#' # Conduct random effects network meta-analysis with
#' # placebo as reference treatment
#' #
#' net1 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'   data = Senn2013.5, sm = "MD", common = FALSE, reference = "plac")
#'       
#' # Generate a heat plot (with abbreviated treatment labels)
#' #
#' heatplot(net1, nchar.trts = 4) 
#' 
#' @keywords hplot
#' 
#' @rdname heatplot
#' @export heatplot

heatplot <- function(x, ...)
  UseMethod("heatplot")
