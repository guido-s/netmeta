#' Hasse diagram
#' 
#' @description
#' This function generates a Hasse diagram for a partial order of
#' treatment ranks in a network meta-analysis.
#' 
#' @aliases hasse hasse.netposet
#' 
#' @param x An object of class \code{netposet} (mandatory).
#' @param pooled A character string indicating whether Hasse diagram
#'   show be drawn for common (\code{"common"}) or random effects model
#'   (\code{"random"}). Can be abbreviated.
#' @param newpage A logical value indicating whether a new figure
#'   should be printed in an existing graphics window. Otherwise, the
#'   Hasse diagram is added to the existing figure.
#' @param shape Shape of node borders, either \code{"roundrect"},
#'   \code{"rect"}, or \code{"none"}, can be abbreviated.
#' @param col.lines Line colour.
#' @param col.nodes Colour for treatment node borders.
#' @param lwd Width of lines and node borders.
#' @param \dots Additional arguments (ignored).
#' 
#' @details
#' Generate a Hasse diagram (Carlsen & Bruggemann, 2014) for a partial
#' order of treatment ranks in a network meta-analysis (Rücker &
#' Schwarzer, 2017).
#' 
#' This R function is a wrapper function for a modified version of R function
#' \code{hasse} in R package \bold{hasseDiagram}
#' (Krzysztof Ciomek, \url{https://github.com/kciomek/hasseDiagram}) which is
#' available under the MIT license.
#' 
#' @author Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}, Guido
#'   Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netposet}},
#'   \code{\link{netrank}}, \code{\link{plot.netrank}},
#'   \code{\link[metabook]{Linde2015}}
#' 
#' @references
#' Carlsen L, Bruggemann R (2014):
#' Partial order methodology: a valuable tool in chemometrics.
#' \emph{Journal of Chemometrics},
#' \bold{28}, 226--34
#' 
#' Rücker G, Schwarzer G (2017):
#' Resolve conflicting rankings of outcomes in network meta-analysis:
#' Partial ordering of treatments.
#' \emph{Research Synthesis Methods},
#' \bold{8}, 526--36
#' 
#' @name hasse.netposet
#' 
#' @keywords plot
#' 
#' @examples
#' \dontrun{
#' # Only run example if R package 'Rgraphviz' from Bioconductor is available
#' #
#' if (requireNamespace("Rgraphviz", quietly = TRUE)) {
#' # Define order of treatments in depression dataset Linde2015
#' #
#' trts <- c("TCA", "SSRI", "SNRI", "NRI",
#'   "Low-dose SARI", "NaSSa", "rMAO-A", "Hypericum", "Placebo")
#' 
#' # Outcome labels
#' #
#' outcomes <- c("Early response", "Early remission")
#' 
#' # (1) Early response
#' #
#' pw1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3),
#'   n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#' #
#' nma1 <- netmeta(pw1, common = FALSE,
#'   seq = trts, ref = "Placebo", small.values = "undesirable")
#' 
#' # (2) Early remission
#' #
#' pw2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3),
#'   n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#' #
#' nma2 <- netmeta(pw2, common = FALSE,
#'   seq = trts, ref = "Placebo", small.values = "undesirable")
#' 
#' # Partial order of treatment rankings
#' #
#' po <- netposet(netrank(nma1), netrank(nma2), outcomes = outcomes)
#' 
#' # Hasse diagram
#' #
#' hasse(po)
#' }
#' }
#' 
#' @method hasse netposet
#' @export

hasse.netposet <- function(x,
                           pooled = ifelse(x$random, "random", "common"),
                           newpage = TRUE,
                           shape = "roundrect",
                           col.lines = "black", col.nodes = "black",
                           lwd = 1,
                           ...) {
  
  if (!is_installed_package("Rgraphviz", stop = FALSE)) {
    message(paste0("Package 'Rgraphviz' missing.",
                   "\n  ",
                   "Please use the following R commands for installation:",
                   "\n  ",
                   "install.packages(\"BiocManager\")",
                   "\n  ",
                   "BiocManager::install(\"Rgraphviz\")"))
    #
    return(invisible(NULL))
  }
  
  chkclass(x, "netposet")
  x <- updateversion(x)
  #
  pooled <- setchar(pooled, c("common", "random", "fixed"))
  pooled[pooled == "fixed"] <- "common"
  #
  chklogical(newpage)
  #
  shape <- setchar(shape, c("roundrect", "rect", "none"))
  #
  chkcolor(col.lines, length = 1)
  chkcolor(col.nodes, length = 1)
  #
  chknumeric(lwd, min = 0, zero = TRUE, length = 1)
  
  if (pooled == "common")
    M <- x$M.common
  else
    M <- x$M.random
  #
  M <- matrix(as.logical(M), dim(M), dimnames = dimnames(M))
  #
  pars <- list(newpage = newpage,  shape = shape,
               col.lines = col.lines, col.nodes = col.nodes,
               lwd = lwd)
  #
  hasse_hasseDiagram(M, pars)
  
  invisible()
}


#' @rdname hasse.netposet
#' @export

hasse <- function(x, ...)
  UseMethod("hasse")
