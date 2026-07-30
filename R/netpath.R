#' TITLE
#' 
#' @description
#' BRIEF DESCRIPTION (1-2 sentence(s)).
#' 
#' @param x A \code{netmeta} object.
#' @param random A logical indicating whether the path algorithm is
#'   based on a random effects model.
#' @param node1 First node.
#' @param node2 Second node.
#' 
#' @return
#' A netpath object.
#' 
#' @author Noosheen R. Tahmasebi
#'   \email{noosheen.rajabzadehtahmasebi@@uniklinik-freiburg.de},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}
#' 
#' @examples
#' \dontrun{
#' data(Senn2013)
#' nma1 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'   data = Senn2013, sm = "MD", random = FALSE, nchar.trts = 4)
#' 
#' np1 <- netpath(nma1, node1 = "Placebo", node2 = "Sulfonylurea")
#' np1
#' }
#' 
#' @export netpath

netpath <- function(x, random = x$random, node1, node2) {
  chkclass(x, "netmeta")
  #
  chklogical(random)
  
  if (random)
    hm <- hatmatrix(x, method = "Davies", type = "full")$random
  else
    hm <- hatmatrix(x, method = "Davies", type = "full")$common
  
  # Run path inconsistency analysis
  #
  res <- run_path_inconsistency(x, hm, node1, node2)
  #
  res
}


#' @method print netpath
#' @export

print.netpath <- function(x, ...) {
  chkclass(x, "netpath")
  #
  cat("---------- Path-based inconsistency test ----------\n")
  cat("Comparison:", unique(x$results$comparison), "\n")
  cat("Number of independent paths:", x$results$df + 1, "\n")
  cat("Quadratic form Q:", x$results$Q, "\n")
  cat("Degrees of freedom:", x$results$df, "\n")
  cat("p-value:", x$results$pval, "\n")
  if (!is.null(x$results$note))
    cat(x$results$note, "\n")
  cat("---------------------------------------------------\n")
  #
  invisible(NULL)
}
