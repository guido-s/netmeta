#' Path-based test of inconsistency in network meta-analysis
#' 
#' @description
#' Performs a path-based test of inconsistency for a specific comparison
#' (treatment pair) in a network meta-analysis, based on the decomposition of
#' the hat matrix into independent paths connecting the two nodes.
#'
#' @param x A \code{netmeta} object.
#' @param random A logical indicating whether the path algorithm is
#'   based on a random effects model.
#' @param node1 First node.
#' @param node2 Second node.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
#' @param sep.trts A character used in comparison names as separator
#'   between treatment labels.
#' @param digits.Q Minimal number of significant digits for
#'   heterogeneity statistics, see \code{print.default}.
#' @param digits.pval.Q Minimal number of significant digits for
#'   p-value of heterogeneity tests, see \code{print.default}.
#' @param details.methods A logical specifying whether details on statistical
#'   methods should be printed.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param \dots Additional arguments (igored).
#'
#' @details
#' This function implements a path-based approach to assess inconsistency in a
#' network meta-analysis (Tahmasebi et al., 2025). Starting from the hat matrix
#' of the network (as calculated by \code{\link{hatmatrix}} with
#' \code{method = "Davies"} and \code{type = "full"}), the direct and indirect
#' evidence contributing to the comparison between \code{node1} and \code{node2}
#' is decomposed into a set of independent evidence paths. A depth-first search
#' algorithm is then used to identify these paths. A test statistic Q is
#' calculated to test whether the estimates derived from these independent paths
#' are consistent with each other. Under the null hypothesis of consistency, Q
#' approximately follows a chi-squared distribution with degrees of freedom
#' equal to the number of independent paths minus one.
#'
#' Depending on the argument \code{random}, the network estimates are based on
#' either the common effects model (\code{random = FALSE}) or the random effects
#' model (\code{random = TRUE}).
#' 
#' In order to get more concise printouts, argument \code{nchar.trts} can be
#' used to define the minimum number of characters for abbreviated treatment
#' names (see \code{\link{abbreviate}}, argument \code{minlength}). R function
#' \code{\link{treats}} is utilised internally to create abbreviated
#' treatment names.
#'
#' @return
#' A list of class "netpath" containing the following elements:
#' \item{results}{Data frame containing results of the Q tests}
#' \item{A.matrix}{Path-adjacency matrix}
#' \item{Sigma}{Standardized matrix derived from the linearly independent paths}
#' \item{theta_p}{Theta p}
#' \item{random, node1, node2}{As defined above.}
#' \item{nchar.trts}{As defined above.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create object.}
#'
#' @author Noosheen R. Tahmasebi
#'   \email{noosheen.rajabzadehtahmasebi@@uniklinik-freiburg.de},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{netmeta}}, \code{\link{heatplot.netpath}}
#'
#' @references
#' Tahmasebi NR, Davies AL, Papakonstantinou T, Rücker G,
#' Nikolakopoulou A (2025):
#' Path-based approach for detecting and assessing inconsistency in network
#' meta-analysis: A novel method.
#' \emph{arXiv}, \doi{https://doi.org/10.48550/arXiv.2506.20364}
#'
#' @examples
#' \dontrun{
#' # Transform data from long arm-based to contrast-based format
#' #
#' pw <- pairwise(studlab = study, treat = treatment,
#'   n = n, mean = mean, sd = sd, data = Senn2013,
#'   varnames = c("MD", "seMD"))
#'
#' # Conduct common effects network meta-analysis
#' #
#' nma <- netmeta(pw, random = FALSE, nchar.trts = 4)
#' 
#' np <- netpath(nma, node1 = "Placebo", node2 = "Sulfonylurea")
#' np
#' }
#'
#' @export netpath

netpath <- function(x, random = x$random, node1, node2,
                    nchar.trts = x$nchar.trts) {
  
  chkclass(x, "netmeta")
  x <- updateversion(x)
  #
  chklogical(random)
  #
  node1 <- setchar(node1, x$trts)
  node2 <- setchar(node2, x$trts)
  #
  nchar.trts <- replaceNULL(nchar.trts, 666)
  chknumeric(nchar.trts, min = 1, length = 1)
  #
  if (random)
    hm <- hatmatrix(x, method = "Davies", type = "full")$random
  else
    hm <- hatmatrix(x, method = "Davies", type = "full")$common
  
  # Run path inconsistency analysis
  #
  res <- run_path_inconsistency(x, hm, node1, node2, x$sep.trts)
  #
  res$random <- random
  res$node1 <- node1
  res$node2 <- node2
  #
  res$trts <- x$trts
  res$nchar.trts <- nchar.trts
  res$sep.trts <- x$sep.trts
  #
  res$call <- match.call()
  res$version = packageDescription("netmeta")$Version
  #
  class(res) <- "netpath"
  #
  res
}


#' @rdname netpath
#' @method print netpath
#' @export

print.netpath <- function(x,
                          #
                          nchar.trts = x$nchar.trts,
                          sep.trts = x$sep.trts,
                          #
                          digits.Q = gs("digits.Q"),
                          digits.pval.Q = gs("digits.pval.Q"),
                          #
                          details.methods = gs("details"),
                          legend = gs("legend"),
                          ...) {
  chkclass(x, "netpath")
  #
  chknumeric(nchar.trts, min = 1, length = 1)
  #
  missing.sep.trts <- missing(sep.trts)
  sep.trts <- replaceNULL(sep.trts, ":")
  chkchar(sep.trts, length = 1)
  sep.trts <- setsep(x$trts, sep.trts, missing = missing.sep.trts)
  #
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.pval.Q, min = 1, length = 1)
  #
  chklogical(details.methods)
  chklogical(legend)
  
  # Get rid of warning 'Undefined global functions or variables'
  path_index <- note <- comparison <- Q <- pval <- NULL
  
  comps <- unique(x$results$comparison)
  comps.abbr <- comps(comps, x$trts, x$sep.trts, nchar.trts)
  #
  trts <- unlist(strsplit(comps, x$sep.trts))
  trts.abbr <- unlist(strsplit(comps.abbr, x$sep.trts))
  #
  if (x$sep.trts != sep.trts)
    comps.abbr <- sub(x$sep.trts, sep.trts, comps.abbr)
  #
  dat <- x$results %>% select(-path_index, -note) %>% unique()
  #
  rownames(dat) <- comps.abbr
  #
  dat %<>% rename(N.paths = comparison) %>%
    mutate(N.paths = nrow(x$results),
           Q = formatN(Q, digits = digits.Q),
           pval = formatPT(pval, digits = digits.pval.Q))
  
  cat("Path-based inconsistency test\n\n")
  #
  print(dat)
  
  if (details.methods) {
    txt.details <- "\nDetails on statistical methods:"
    #
    txt.details <- paste0(txt.details,
                          paste0("\n- ", if (x$random) "Random" else "Common",
                                 " effects network meta-analysis"))
    #
    if (!is.null(x$results$note) || all(x$results$note == ""))
      txt.details <- paste0(txt.details, paste0("\n- ", unique(x$results$note)))
    #
    cat(paste0(txt.details, "\n"))
  }
  #
  # Add legend with abbreviated treatment labels
  #
  legendabbr(trts, trts.abbr, legend)
  #
  invisible(NULL)
}
