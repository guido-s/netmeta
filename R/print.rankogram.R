#' Print method for rankograms
#' 
#' @description
#' Print method for objects of class \code{rankogram}.
#' 
#' @param x An R object of class \code{rankogram}.
#' @param common A logical indicating to print ranking probabilities
#'   and SUCRAs for the common effects model.
#' @param random A logical indicating to print ranking probabilities
#'   and SUCRAs for the random effects model.
#' @param cumulative.rankprob A logical indicating whether cumulative
#'   ranking probabilities should be printed.
#' @param sort A logical indicating whether treatments should be
#'   sorted by decreasing SUCRAs.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments for printing.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{rankogram}}, \code{\link{plot.rankogram}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP (2011):
#' Graphical methods and numerical summaries for presenting results
#' from multiple-treatment meta-analysis: an overview and tutorial.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{64}, 163--71
#'
#' @examples
#' data(Woods2010)
#' p1 <- pairwise(treatment, event = r, n = N, studlab = author,
#'   data = Woods2010, sm = "OR")
#' net1 <- netmeta(p1, small.values = "desirable")
#'
#' ran1 <- rankogram(net1, nsim = 100)
#' ran1
#' print(ran1, cumulative.rankprob = TRUE)
#'
#' @rdname print.rankogram
#' @method print rankogram
#' @export

print.rankogram <- function(x,
                            common = x$common,
                            random = x$random,
                            cumulative.rankprob = x$cumulative.rankprob,
                            sort = TRUE,
                            nchar.trts = x$nchar.trts,
                            digits = gs("digits.prop"),
                            legend = gs("legend"),
                            warn.deprecated = gs("warn.deprecated"),
                            ...) {
  
  #
  #
  # (1) Check for rankogram object and upgrade object
  #
  #
  
  chkclass(x, "rankogram")
  x <- updateversion(x)
  
  
  #
  #
  # (2) Check other arguments
  #
  #
  
  chklogical(cumulative.rankprob)
  chklogical(sort)
  #
  chknumeric(nchar.trts, length = 1)
  #
  chknumeric(digits, length = 1)
  chklogical(legend)
  #
  # Check for deprecated arguments in '...'
  #
  args  <- list(...)
  chklogical(warn.deprecated)
  #
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  #
  random <-
    deprecated(random, missing(random), args, "comb.random", warn.deprecated)
  chklogical(random)
  
  
  #
  #
  # (3) Print results
  #
  #
  
  if (common | random)
    cat(if (cumulative.rankprob)
      "Cumulative ranking probabilities" else "Rankogram",
      " (based on ", x$nsim, " simulation",
      if (x$nsim > 1) "s", ")\n\n",
      sep = "")
  #
  if (common) {
    if (cumulative.rankprob)
      rank.common <- x$cumrank.matrix.common
    else
      rank.common <- x$ranking.matrix.common
    rownames(rank.common) <- treats(rank.common, nchar.trts)
    #
    if (sort)
      rank.common <-
      rank.common[rev(order(x$ranking.common)), , drop = FALSE]
    #
    cat("Common effects model: \n\n")
    prmatrix(formatN(rank.common, digits), quote = FALSE, right = TRUE, ...)
    if (random)
      cat("\n")
  }
  #
  if (random) {
    if (cumulative.rankprob)
      rank.random <- x$cumrank.matrix.random
    else
      rank.random <- x$ranking.matrix.random
    rownames(rank.random) <-
      treats(rank.random, nchar.trts)
    #
    if (sort)
      rank.random <-
      rank.random[rev(order(x$ranking.random)), , drop = FALSE]
    #
    cat("Random effects model: \n\n")
    prmatrix(formatN(rank.random, digits), quote = FALSE, right = TRUE, ...)
  }
  #
  # Add legend with abbreviated treatment labels
  #
  if ((common | random) & legend) {
    if (common)
      trts <- rownames(x$ranking.matrix.common)
    else if (random)
      trts <- rownames(x$ranking.matrix.random)
    #
    legendabbr(trts, treats(trts, nchar.trts), TRUE)
  }
  
  invisible()
}
