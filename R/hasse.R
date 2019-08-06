#' Hasse diagram
#' 
#' @description
#' This function generates a Hasse diagram for a partial order of
#' treatment ranks in a network meta-analysis.
#' 
#' @param x An object of class \code{netposet} (mandatory).
#' @param pooled A character string indicating whether Hasse diagram
#'   show be drawn for fixed effect (\code{"fixed"}) or random effects
#'   model (\code{"random"}). Can be abbreviated.
#' @param newpage A logical value indicating whether a new figure
#'   should be printed in an existing graphics window. Otherwise, the
#'   Hasse diagram is added to the existing figure.
#' 
#' @details
#' Generate a Hasse diagram (Carlsen & Bruggemann, 2014) for a partial
#' order of treatment ranks in a network meta-analysis (Rücker &
#' Schwarzer, 2017).
#' 
#' This R function is a wrapper function for R function
#' \code{\link[hasseDiagram]{hasse}} in R package \bold{hasseDiagram}
#' (Krzysztof Ciomek, \url{https://github.com/kciomek/hasseDiagram}),
#' i.e., function \code{hasse} can only be used if R package
#' \bold{hasseDiagram} is installed.
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netposet}},
#'   \code{\link{netrank}}, \code{\link{plot.netrank}}
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
#' @keywords plot
#' 
#' @examples
#' # Use depression dataset
#' #
#' data(Linde2015)
#' 
#' # Define order of treatments
#' #
#' trts <- c("TCA", "SSRI", "SNRI", "NRI",
#'           "Low-dose SARI", "NaSSa", "rMAO-A", "Hypericum",
#'           "Placebo")
#' 
#' # Outcome labels
#' #
#' outcomes <- c("Early response", "Early remission")
#' 
#' # (1) Early response
#' #
#' p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'                event = list(resp1, resp2, resp3),
#'                n = list(n1, n2, n3),
#'                studlab = id, data = Linde2015, sm = "OR")
#' #
#' net1 <- netmeta(p1, comb.fixed = FALSE,
#'                 seq = trts, ref = "Placebo")
#' 
#' # (2) Early remission
#' #
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'                event = list(remi1, remi2, remi3),
#'                n = list(n1, n2, n3),
#'                studlab = id, data = Linde2015, sm = "OR")
#' #
#' net2 <- netmeta(p2, comb.fixed = FALSE,
#'                 seq = trts, ref = "Placebo")
#' 
#' # Partial order of treatment rankings
#' #
#' po <- netposet(netrank(net1, small.values = "bad"),
#'                netrank(net2, small.values = "bad"),
#'                outcomes = outcomes)
#' 
#' # Hasse diagram
#' #
#' hasse(po)
#' 
#' @export hasse


hasse <- function(x,
                  pooled = ifelse(x$comb.random, "random", "fixed"),
                  newpage = TRUE) {
  
  meta:::chkclass(x, "netposet")
  ##
  pooled <- meta:::setchar(pooled, c("fixed", "random"))
  
  
  if (!meta:::is.installed.package("hasseDiagram", stop = FALSE))
    stop(paste("Library 'hasseDiagram' missing.",
               "\n  ",
               "Please use the following R commands for installation:",
               "\n  ",
               "install.packages(\"BiocManager\")",
               "\n  ",
               "BiocManager::install()",
               "\n  ",
               "BiocManager::install(\"Rgraphviz\")",
               "\n  ",
               "install.packages(\"hasseDiagram\")",
               sep = ""),
         call. = FALSE)
  
  if (pooled == "fixed")
    M <- x$M.fixed
  else
    M <- x$M.random
  ##
  hasseDiagram::hasse(matrix(as.logical(M), dim(M), dimnames = dimnames(M)),
                      parameters = list(newpage = newpage))
  
  
  invisible(NULL)
}
