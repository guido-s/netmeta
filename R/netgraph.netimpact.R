#' Network graph for objects of class netimpact
#' 
#' @description
#' This function generates a graph of the evidence network.
#' 
#' @param x An object of class \code{netimpact}.
#' @param col.ignore A character string indicating color for
#'   comparisons removed from network, either \code{"transparent"} or
#'   any color defined in \code{\link[grDevices]{colours}}.
#' @param number.of.studies A logical indicating whether number of
#'   studies should be added to network graph.
#' @param main Main title.
#' @param sub Subtitle.
#' @param multiarm A logical indicating whether multi-arm studies
#'   should be marked in plot.
#' @param col.multiarm Either a function from R library colorspace or
#'   grDevice to define colors for multi-arm studies or a character
#'   vector with colors to highlight multi-arm studies.
#' @param alpha.transparency The alpha transparency of colors used to
#'   highlight multi-arm studies (0 means transparent and 1 means
#'   opaque).
#' @param col.ignore.multiarm A character string indicating color to
#'   mark multi-arm studies removed from network, either
#'   \code{"transparent"} or any color defined in
#'   \code{\link[grDevices]{colours}}.
#' @param col A single color (or vector of colors) for lines
#'   connecting treatments (edges) if argument \code{plastic =
#'   FALSE}. Length of the vector must be equal to the number of
#'   edges.
#' @param \dots Additional arguments passed on to
#'   \code{\link{netgraph.netmeta}}.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de},
#'   Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netimpact}}, \code{\link{netgraph.netmeta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' data(parkinson)
#' 
#' p1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'                n = list(n1, n2, n3),
#'                mean = list(y1, y2, y3),
#'                sd = list(sd1, sd2, sd3),
#'                data = parkinson, studlab = Study)
#' 
#' net1 <- netmeta(p1)
#' ni <- netimpact(net1, verbose = TRUE)
#' netgraph(ni, plastic = FALSE)
#' 
#' @method netgraph netimpact
#' @export
#' @export netgraph.netimpact


netgraph.netimpact <- function(x,
                               col.ignore = "red",
                               number.of.studies = TRUE,
                               main, sub,
                               multiarm = any(x$x$narms > 2),
                               col.multiarm = NULL,
                               alpha.transparency = 0.5,
                               col.ignore.multiarm =  "transparent",
                               col = "slateblue",
                               ...) {
  
  
  meta:::chkclass(x, "netimpact")
  
  
  col.ignore <- meta:::setchar(col.ignore,
                               c("transparent", colours()),
                               text = paste0("should be any color ",
                                             "defined in colours()"))
  ##
  col.ignore.multiarm <- meta:::setchar(col.ignore.multiarm,
                                        c("transparent", colours()),
                                        text = paste0("should be any color ",
                                                      "defined in colours()"))
  
  
  studlab <- x$x$studlab
  treat1  <- x$x$treat1
  treat2  <- x$x$treat2
  TE      <- x$x$TE
  seTE    <- x$x$seTE
  ##
  sep.trts <- x$x$sep.trts
  comparison <- paste(treat1, sep = sep.trts, treat2)
  
  
  comparisons <- names(x$x$prop.direct.fixed)
  studies <- x$x$studies
  narms <- x$x$narms
  ##  
  impact <- matrix(NA, ncol = x$x$k, nrow = length(comparisons))
  ##
  rownames(impact) <- comparisons
  colnames(impact) <- studies
  ##
  ## Run network meta-analyses "excluding" one study at a time
  ##
  if (multiarm) {
    mc <- multicols(studies, narms,
                    missing(col.multiarm),
                    col.multiarm, alpha.transparency)
    col.polygon <- mc$cols
    multiarm.studies <- mc$multiarm.studies
  }
  else
    col.polygon <- col.ignore.multiarm
  ##
  res <- list()
  ##
  for (i in studies) {
    ##
    seTE.i <- seTE
    seTE.i[studlab == i] <- x$seTE.ignore
    ##
    ignore.i <- x$ignored.comparisons[[i]]
    col.ignore.i <- rep(col.ignore, length(ignore.i))
    ##
    mat <- matrix(unlist(strsplit(ignore.i, split = sep.trts)),
                  ncol = 2, byrow = TRUE)
    treat1.i <- mat[, 1]
    treat2.i <- mat[, 2]
    ##
    col.polygon.i <- col.polygon
    ##
    if (multiarm)
      col.polygon.i[multiarm.studies == i] <- col.ignore.multiarm
    ##
    net.i <- x$nets[[i]]
    ##
    for (j in seq_along(treat1.i)) {
      ##
      net.i$A.matrix[treat1.i[j], treat2.i[j]] <-
        net.i$A.matrix[treat1.i[j], treat2.i[j]] - 1
      net.i$A.matrix[treat2.i[j], treat1.i[j]] <-
        net.i$A.matrix[treat2.i[j], treat1.i[j]] - 1
      ##
      if (net.i$A.matrix[treat1.i[j], treat2.i[j]] == 0)
        col.ignore.i[j] <- "transparent"
    }
    ##
    n.i <- netgraph(net.i,
                    highlight = ignore.i, col.highlight = col.ignore.i,
                    multiarm = multiarm, col.multiarm = col.polygon.i,
                    alpha.transparency = alpha.transparency,
                    number.of.studies = number.of.studies,
                    col = col, ...)
    ##
    res[[i]] <- list(nodes = n.i$nodes, edges = n.i$edges)
    ##
    if (!missing(main)) {
      if (!(is.logical(main) && length(main) == 1 && !main))
        title(main = main)
    }
    else
      title(main = paste0("Study removed: ", i))
    ##
    if (!missing(sub)) {
      if (!(is.logical(sub) && length(sub) == 1 && !sub))
        title(sub = sub)
    }
    else
      title(sub = paste0("Comparison",
                         if (length(ignore.i) > 1) "s",
                         ": ",
                         paste(paste("'", ignore.i, "'", sep = ""),
                               collapse = ", ")))
  }
  
  
  invisible(res)
}
