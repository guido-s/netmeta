#' Print method for objects of class netimpact
#' 
#' @description
#' Print method for objects of class \code{netimpact}.
#' 
#' @param x An object of class \code{netimpact}.
#' @param common A logical indicating whether results for the common
#'   effects model should be printed.
#' @param random A logical indicating whether results for the random
#'   effects model should be printed.
#' @param digits Minimal number of significant digits.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
#' @param nchar.studlab A numeric defining the minimum number of
#'   characters used to create unique study labels.
#' @param details.methods A logical specifying whether details on statistical
#'   methods should be printed.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param legend.studlab A logical indicating whether a legend should
#'   be printed for abbreviated study labels.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netimpact}}
#' 
#' @keywords print
#' 
#' @examples
#' data(Franchini2012)
#'
#' # Only consider first two studies (to reduce runtime of example)
#' #
#' studies <- unique(Franchini2012$Study)
#' p1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'   n = list(n1, n2, n3),
#'   mean = list(y1, y2, y3), sd = list(sd1, sd2, sd3),
#'   data = subset(Franchini2012, Study %in% studies[1:2]),
#'   studlab = Study)
#' 
#' net1 <- netmeta(p1)
#' ni <- netimpact(net1, verbose = TRUE)
#' ni
#' 
#' @method print netimpact
#' @export


print.netimpact <- function(x,
                            common = x$x$common,
                            random = x$x$random,
                            #
                            digits = gs("digits.prop"),
                            #
                            nchar.trts = x$nchar.trts,
                            nchar.studlab = x$nchar.studlab,
                            #
                            details.methods = gs("details"),
                            legend = gs("legend"),
                            legend.studlab = TRUE,
                            warn.deprecated = gs("warn.deprecated"),
                            ##
                            ...) {
  
  ##
  ##
  ## (1) Check for netimpact object and upgrade object
  ##
  ##
  chkclass(x, "netimpact")
  x <- updateversion(x)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chknumeric(digits, min = 0, length = 1)
  chklogical(details.methods)
  chklogical(legend)
  chklogical(legend.studlab)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                      warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  
  
  ##
  ##
  ## (3) Generate (abbreviated) column names
  ##
  ##
  sep.trts <- x$x$sep.trts
  ##
  cn <- colnames(x$impact.common)
  mat <- matrix(unlist(strsplit(cn, split = sep.trts)),
                ncol = 2, byrow = TRUE)
  treat1.long <- mat[, 1]
  treat2.long <- mat[, 2]
  ##
  trts <- x$x$trts
  ##
  treat1 <- as.character(factor(treat1.long, levels = trts,
                                labels = treats(trts, nchar.trts)))
  treat2 <- as.character(factor(treat2.long, levels = trts,
                                labels = treats(trts, nchar.trts)))
  
  
  ##
  ##
  ## (4) Print results
  ##
  ##
  if (common) {
    cat("Common effects model: \n\n")
    impact.common <- formatN(x$impact.common, digits = digits)
    #
    studlab <- rownames(impact.common)
    studlab.abbr <- treats(studlab, nchar.studlab)
    rownames(impact.common) <- studlab.abbr
    #
    colnames(impact.common) <- paste(treat1, treat2, sep = sep.trts)
    ##
    prmatrix(impact.common, quote = FALSE, right = TRUE)
    if (random)
      cat("\n")
  }
  ##
  ## Print results for random effects model
  ##
  if (random) {
    cat("Random effects model: \n\n")
    impact.random <- formatN(x$impact.random, digits = digits)
    #
    studlab <- rownames(impact.random)
    studlab.abbr <- treats(studlab, nchar.studlab)
    rownames(impact.random) <- studlab.abbr
    #
    colnames(impact.random) <- paste(treat1, treat2, sep = sep.trts)
    ##
    prmatrix(impact.random, quote = FALSE, right = TRUE)
  }
  #
  # Print details of network meta-analysis methods
  #
  if (details.methods & (common | random))
    cat(textmeth(x, random))
  ##
  ## Add legend with abbreviated treatment labels
  ##
  if (legend & (common | random))
    legendabbr(trts, treats(trts, nchar.trts), legend)
  #
  # Add legend with abbreviated study labels
  #
  if (legend.studlab & (common | random))
    legendabbr(studlab, studlab.abbr, TRUE, "Study label",
               if (legend) "\n" else "\nLegend:\n")
  
  
  invisible(NULL)
}
