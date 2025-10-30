#' Heat plot for network path
#'
#' @description
#' BRIEF DESCRIPTION (1-2 sentence(s)).
#'
#' @param x An object of class \code{netmeta}.
#' @param \dots Additional arguments (ignored).
#' 
#' @return
#' Heat plot ...
#' 
#' @keywords hplot
#' 
#' @author Noosheen R. Tahmasebi
#'   \email{noosheen.rajabzadehtahmasebi@@uniklinik-freiburg.de},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @examples
#' \dontrun{
#' data(Senn2013)
#' nma1 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'   data = Senn2013, sm = "MD", random = FALSE, nchar.trts = 4)
#' hm <- hatmatrix(nma1, method = "Davies", type = "full")
#' 
#' np1 <- netpath(nma1, hm, node1 = "Placebo", node2 = "Sulfonylurea")
#' heatplot(np1)
#' }
#'
#' @method heatplot netpath
#' @export

heatplot.netpath <- function(x, ...) {
  chkclass(x, "netpath")
  #
  xhat <- x$theta_p
  Sigma <- x$S
  index <- x$results$path_index
  #
  create_standardized_heatmap(xhat, Sigma, index)  
}
