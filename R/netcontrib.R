#' Contribution matrix in network meta-analysis
#' 
#' @description
#' This function generates the contribution of direct comparisons to
#' every network treatment comparison as a different row
#' 
#' @aliases netcontrib print.netcontrib
#' 
#' @param x An object of class \code{netmeta} or \code{netcontrib}.
#' @param method A character string indicating which method is to
#'   calculate the contribution matrix. Either \code{"stream"} or
#'   \code{"randomwalk"}, can be abbreviated.
#' @param hatmatrix.F1000 A logical indicating whether hat matrix
#'   given in F1000 article should be used for \code{method =
#'   "stream"}.
#' @param comb.fixed A logical indicating whether a contribution
#'   matrix should be printed for the fixed effects (common effects)
#'   network meta-analysis.
#' @param comb.random A logical indicating whether a contribution
#'   matrix should be printed for the random effects network
#'   meta-analysis.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param \dots Additional arguments (ignored at the moment).
#' 
#' @details
#' In network meta-analysis, it is important to assess the influence
#' of the limitations or other characteristics of individual studies
#' on the estimates obtained from the network. The contribution
#' matrix, shows how much each direct treatment effect contributes to
#' each treatment effect estimate from network meta-analysis, is
#' crucial in this context.
#' 
#' We use ideas from graph theory to derive the proportion that is
#' contributed by each direct treatment effect.  We start with the
#' ‘projection’ matrix in a two-step network meta-analysis model,
#' called the H matrix, which is analogous to the hat matrix in a
#' linear regression model.
#'
#' Two methods are implemented to estimate network contributions. If
#' argument \code{method = "stream"}, H entries are translated to
#' proportion contributions based on the observation that the rows of
#' H can be interpreted as flow networks, where a stream is defined as
#' the composition of a path and its associated flow (Papakonstantinou
#' et al., 2018). If argument \code{method = "randomwalk"}, a random
#' walk algorithm is used (Davies et al., 2021).
#' 
#' By default, treatment names are not abbreviated in
#' printouts. However, in order to get more concise printouts,
#' argument \code{nchar.trts} can be used to define the minimum number
#' of characters for abbreviated treatment names (see
#' \code{\link{abbreviate}}, argument \code{minlength}). R function
#' \code{\link{treats}} is utilised internally to create abbreviated
#' treatment names.
#' 
#' @return
#' An object of class \code{netcontrib} with corresponding
#' \code{print} function. The object is a list containing the
#' following components:
#' \item{fixed}{Numeric matrix of percentage contributions of direct
#'   comparisons for each network comparison for the fixed effects
#'   model.}
#' \item{random}{Numeric matrix of percentage contributions of direct
#'   comparisons for each network comparison for the random effects
#'   model.}
#' \item{x}{As defined above.}
#' with the contribution matrices for fixed and random NMA. Each
#' matrix has the percentage contributions of each direct comparison
#' as columns for each network comparison, direct or indirect as rows.
#' 
#' @author Theodoros Papakonstantinou \email{dev@tpapak.com}, Annabel
#'   Davies \email{annabel.davies@manchester.ac.uk}
#' 
#' @seealso \code{\link{netmeta}}
#' 
#' @references
#' Davies AL, Papakonstantinou T, Nikolakopoulou A, Rücker G, Galla T
#' (2021):
#' Network meta-analysis and random walks.
#' Available from: http://arxiv.org/abs/2107.02886
#' 
#' Papakonstantinou, T., Nikolakopoulou, A., Rücker, G., Chaimani, A.,
#' Schwarzer, G., Egger, M., Salanti, G. (2018):
#' Estimating the contribution of studies in network meta-analysis:
#' paths, flows and streams.
#' \emph{F1000Research}
#' 
#' @keywords contribution
#' 
#' @examples
#' # Use the Woods dataset
#' #
#' data("Woods2010")
#' p1 <- pairwise(treatment, event = r, n = N,
#'                studlab = author, data = Woods2010, sm = "OR")
#' 
#' net1 <- netmeta(p1)
#' cm <- netcontrib(net1)
#' cm
#'
#' netcontrib(net1, method = "r")
#' 
#' @rdname netcontrib
#' @export netcontrib


netcontrib <- function(x,
                       method = "stream",
                       hatmatrix.F1000 = TRUE,
                       comb.fixed = x$comb.fixed,
                       comb.random = x$comb.random,
                       nchar.trts = x$nchar.trts) {
  
  meta:::is.installed.package("igraph")
  ##
  method <- meta:::setchar(method, c("randomwalk", "stream"))
  meta:::chklogical(hatmatrix.F1000)
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  meta:::chknumeric(nchar.trts, min = 1, length = 1)
  
  
  res <- list(fixed =
                contribution.matrix(x, method, "fixed", hatmatrix.F1000),
              random =
                contribution.matrix(x, method, "random", hatmatrix.F1000),
              method = method,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              nchar.trts = nchar.trts,
              x = x)
  ##
  class(res) <- "netcontrib"
  ##
  res
}





#' @rdname netcontrib
#' 
#' @method print netcontrib
#' 
#' @export
#' @export print.netcontrib


print.netcontrib <- function(x,
                             comb.fixed = x$comb.fixed,
                             comb.random = x$comb.random,
                             digits = 4,
                             nchar.trts = x$nchar.trts,
                             legend = TRUE,
                             ...) {
  
  meta:::chkclass(x, "netcontrib")
  ##
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  meta:::chknumeric(nchar.trts, length = 1)
  meta:::chklogical(legend)
  
  
  matitle(x$x)
  ##
  cat(paste0("Contribution matrix (",
             if (is.null(x$method) | x$method == "stream")
               "Papakonstantinou et al., 2018, F1000Research"
             else
               "Davies et al., 2021",
             ")\n\n"))
  ##
  trts <- x$x$trts
  trts.abbr <- treats(trts, nchar.trts)
  ##
  if (comb.fixed) {
    rownames(x$fixed) <- comps(x$fixed, trts, x$x$sep.trts, nchar.trts)
    colnames(x$fixed) <- comps(x$fixed, trts, x$x$sep.trts, nchar.trts,
                               row = FALSE)
    ##
    cat("Fixed effects model:\n\n")
    prmatrix(round(x$fixed, digits))
    if (comb.random)
      cat("\n")
  }
  if (comb.random) {
    rownames(x$random) <- comps(x$random, trts, x$x$sep.trts, nchar.trts)
    colnames(x$random) <- comps(x$random, trts, x$x$sep.trts, nchar.trts,
                                row = FALSE)
    ##
    cat("Random effects model:\n\n")
    prmatrix(round(x$random, digits))
  }
  ##
  ## Add legend
  ##
  if (legend & (comb.fixed | comb.random)) {
    if (any(trts != trts.abbr)) {
      tmat <- data.frame(trts, trts.abbr)
      names(tmat) <- c("Abbreviation", "Treatment name")
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat("\nLegend:\n")
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(trts.abbr)))
    }
  }
  ##
  invisible(NULL)
}
