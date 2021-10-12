#' Scatter plot or biplot showing partially order of treatment ranks
#' 
#' @description
#' This function generates a scatter plot or biplot of P-scores with
#' an overlay describing partial order of treatment ranks.
#' 
#' @param x An object of class \code{netmeta} (mandatory).
#' @param plottype A character string indicating whether a scatter
#'   plot or biplot should be produced, either \code{"scatter"} or
#'   \code{"biplot"}. Can be abbreviated.
#' @param pooled A character string indicating whether scatter plot
#'   should be drawn for fixed (\code{"fixed"}) or random effects
#'   model (\code{"random"}). Can be abbreviated.
#' @param dim A character string indicating whether a 2- or
#'   3-dimensional plot should be produced, either \code{"2d"} or
#'   \code{"3d"}. Can be abbreviated.
#' @param sel.x A numeric specifying number of outcome to use for the
#'   x-axis in a scatterplot (argument is not considered for a
#'   biplot).
#' @param sel.y A numeric specifying number of outcome to use for the
#'   y-axis in a scatterplot (argument is not considered for a
#'   biplot).
#' @param sel.z A numeric specifying number of outcome to use for the
#'   z-axis in a scatterplot (argument is not considered for a
#'   biplot).
#' @param cex The magnification to be used for treatment labels and
#'   points.
#' @param col Colour(s) of treatment labels and points.
#' @param cex.text The magnification to be used for treatment labels.
#' @param col.text Colour(s) of treatment labels.
#' @param adj.x Value(s) in [0, 1] to specify adjustment of treatment
#'   labels on x-axis (only considered in 2-D plots); see
#'   \code{\link[graphics]{text}}.
#' @param adj.y Value(s) in [0, 1] to specify adjustment of treatment
#'   labels on y-axis (only considered in 2-D plots); see
#'   \code{\link[graphics]{text}}.
#' @param offset.x Offset(s) of treatment labels on x-axis (only
#'   considered in 2-D plots).
#' @param offset.y Offset(s) of treatment labels on y-axis (only
#'   considered in 2-D plots).
#' @param pch Plot symbol(s) for points; no points printed if equal to
#'   \code{NULL}.
#' @param cex.points Magnification(s) to be used for points.
#' @param col.points Colour(s) of points.
#' @param col.lines Line colour.
#' @param lty.lines Line type.
#' @param lwd.lines Line width.
#' @param arrows A logical indicating whether arrows should be printed
#'   (only considered in 2-D plots).
#' @param length Length of arrows; see \code{\link[graphics]{arrows}}.
#' @param grid A logical indicating whether grid lines should be added
#'   to plot.
#' @param col.grid Colour of grid lines.
#' @param lty.grid Line type of grid lines.
#' @param lwd.grid Line width of grid lines.
#' @param \dots Additional graphical arguments.
#' 
#' @details
#' By default (arguments \code{plottype = "scatter"} and \code{dim =
#' "2d"}), a scatter plot is created showing P-scores (see
#' \code{\link{netrank}}) for the first two outcomes considered in the
#' generation of a partially ordered set of treatment ranks (using
#' \code{\link{netposet}}). In addition to the P-scores, the partially
#' order of treatment ranks is shown as lines connecting treatments
#' which is analogous to a Hasse diagram. If argument \code{dim =
#' "3d"}), a 3-D scatter plot is generated showing P-scores for the
#' first three outcomes.
#' 
#' To overcome the restriction of two or three dimension, a biplot
#' (Gabriel, 1971) can be generated using argument \code{plottype =
#' "biplot"}. This is essentially a scatter plot using the first two
#' (\code{dim = "2d"}) or three (\code{dim = "3d"}) components in a
#' principal components analysis (using
#' \code{\link[stats]{prcomp}}). Note, if only two / three outcomes
#' are considered in a \code{netposet} object, a 2-D / 3-D scatter
#' plot is generated instead of a biplot as a principal component
#' analysis is superfluous in such a situation.
#' 
#' Arguments \code{sel.x} and \code{sel.y} can be used to select
#' different outcomes to show on x- and y-axis in a 2-D scatter plot;
#' argument \code{sel.z} can be used accordingly in a 3-D scatter
#' plot.  These arguments are ignored for a biplot.
#' 
#' Note, in order to generate 3-D plots (argument \code{dim = "3d"}),
#' R package \bold{rgl} is necessary. Note, under macOS the X.Org X
#' Window System must be available (see
#' \url{https://www.xquartz.org}).
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netposet}}, \code{\link{hasse}},
#'   \code{\link{netrank}}, \code{\link{netmeta}}
#'   
#' @references
#' Carlsen L, Bruggemann R (2014):
#' Partial order methodology: a valuable tool in chemometrics.
#' \emph{Journal of Chemometrics},
#' \bold{28}, 226--34
#' 
#' Gabriel KR (1971):
#' The biplot graphic display of matrices with application to
#' principal component analysis.
#' \emph{Biometrika},
#' \bold{58}, 453--67
#' 
#' @keywords hplot
#' 
#' @examples
#' \dontrun{
#' data(Linde2015)
#' 
#' # Define order of treatments
#' #
#' trts <- c("TCA", "SSRI", "SNRI", "NRI",
#'           "Low-dose SARI", "NaSSa", "rMAO-A", "Hypericum",
#'           "Placebo")
#' #
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
#' net1 <- netmeta(p1, fixed = FALSE,
#'                 seq = trts, ref = "Placebo", small.values = "bad")
#' 
#' # (2) Early remission
#' #
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'                event = list(remi1, remi2, remi3),
#'                n = list(n1, n2, n3),
#'                studlab = id, data = Linde2015, sm = "OR")
#' #
#' net2 <- netmeta(p2, fixed = FALSE,
#'                 seq = trts, ref = "Placebo", small.values = "bad")
#' 
#' # Partial order of treatment rankings
#' #
#' po2 <- netposet(netrank(net1), netrank(net2), outcomes = outcomes)
#' 
#' # Scatter plot
#' #
#' plot(po2)
#' 
#' # Same scatter plot as only two outcomes considered in netposet()
#' #
#' plot(po2, "biplot")
#'
#' 
#' # Consider three outcomes
#' #
#' # Outcome labels
#' #
#' outcomes <- c("Early response", "Early remission", "Lost to follow-up")
#' 
#' # (3) Loss to follow-up
#' #
#' p3 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'                event = list(loss1, loss2, loss3),
#'                n = list(n1, n2, n3),
#'                studlab = id, data = Linde2015, sm = "OR")
#' #
#' net3 <- netmeta(p3, fixed = FALSE,
#'                 seq = trts, ref = "Placebo", small.values = "good")
#' 
#' # Partial order of treatment rankings (with three outcomes) 
#' #
#' po3 <- netposet(netrank(net1), netrank(net2), netrank(net3),
#'                 outcomes = outcomes)
#' 
#' # Hasse diagram
#' #
#' hasse(po3)
#' 
#' # Scatter plot
#' #
#' plot(po3)
#' 
#' # Biplot (reverse limits of y-axis as biplot is upside down)
#' #
#' plot(po3, "bi", xlim = c(-1, 1.7), ylim = c(2.5, -2.5))
#' }
#' 
#' @method plot netposet
#' @export


plot.netposet <- function(x,
                          plottype = "scatter",
                          pooled = ifelse(x$random, "random", "fixed"),
                          dim = "2d",
                          sel.x = 1, sel.y = 2, sel.z = 3,
                          cex = 1, col = "black",
                          cex.text = cex, col.text = col,
                          adj.x = 0, adj.y = 1,
                          offset.x = 0.005, offset.y = -0.005,
                          pch = NULL, cex.points = cex, col.points = col,
                          col.lines = "black", lty.lines = 1, lwd.lines = 1,
                          arrows = FALSE,
                          length = 0.05,
                          grid = TRUE,
                          col.grid = "gray", lty.grid = 2, lwd.grid = 1,
                          ...) {
  
  
  chkclass(x, "netposet")
  x <- updateversion(x)
  ##
  pooled <- setchar(pooled, c("fixed", "random"))
  
  
  if (pooled == "fixed") {
    p.matrix <- x$P.fixed
    M0 <- x$M0.fixed
  }
  else {
    p.matrix <- x$P.random
    M0 <- x$M0.random
  }
  ##
  outcomes <- colnames(p.matrix)
  treatments <- rownames(p.matrix)
  ##
  n.outcomes   <- length(outcomes)
  n.treatments <- length(treatments)
  
  
  dim <- setchar(dim, c("2d", "3d"))
  is_2d <- dim == "2d"
  is_3d <- !is_2d
  ##
  if (is_3d & n.outcomes == 2) {
    warning("Scatter plot in 2-D generated as only two outcomes considered ",
            "in netposet().")
    is_2d <- TRUE
    is_3d <- FALSE
  }
  ##
  plottype <- setchar(plottype, c("scatter", "biplot"))
  is_biplot <- plottype == "biplot"
  ##
  if (is_biplot) {
    if (is_2d & n.outcomes == 2) {
      is_biplot <- FALSE
      warning("Scatter plot instead of biplot generated as only ",
              "two outcomes considered in netposet().")
    }
    if (is_3d & n.outcomes == 3) {
      is_biplot <- FALSE
      warning("Scatter plot instead of biplot generated as only ",
              "three outcomes considered in netposet().")
    }
  }
  ##
  if (is_biplot) {
    outcomes <- paste("Principal Component", 1:3)
    sel.x <- 2
    sel.y <- 1
    ##
    if (n.outcomes > 2 & is_3d) {
      sel.x <- 2
      sel.y <- 3
      sel.z <- 1
    }
    p.matrix <- prcomp(p.matrix, scale = TRUE)$x
  }
  else {
    sel.x <- as.numeric(setchar(sel.x, seq_len(n.outcomes)))
    sel.y <- as.numeric(setchar(sel.y, seq_len(n.outcomes)))
    ##
    if (n.outcomes > 2)
      sel.z <- as.numeric(setchar(sel.z, seq_len(n.outcomes)))
  }
  ##
  if (is_3d & !is.installed.package("rgl", stop = FALSE)) {
    warning(paste0("2-D plot generated as package 'rgl' is missing.",
                   "\n  ",
                   "Please install package 'rgl' in order to ",
                   "produce 3-D plots\n  ",
                   "(R command: 'install.packages(\"rgl\")').",
                   if (length(grep("darwin", R.Version()$os)) == 1)
                     paste0("\n  Note, macOS users have to install ",
                            "XQuartz, see https://www.xquartz.org/.")
                   ))
    dim <- "2d"
    is_2d <- TRUE
    is_3d <- FALSE
  }
  ##
  use_pch <- !is.null(pch)
  ##
  chknumeric(cex.text, min = 0, zero = 0)
  chknumeric(adj.x)
  chknumeric(adj.y)
  chknumeric(offset.x)
  chknumeric(offset.y)
  chknumeric(lty.lines, min = 0, zero = 0, length = 1)
  chknumeric(lwd.lines, min = 0, zero = 0, length = 1)
  chklogical(arrows)
  chknumeric(length, min = 0, zero = 0, length = 1)
  if (use_pch)
    chknumeric(cex.points, min = 0, zero = 0)
  chklogical(grid)
  chknumeric(lty.grid, min = 0, zero = 0, length = 1)
  chknumeric(lwd.grid, min = 0, zero = 0, length = 1)
  
  
  if (!(length(cex.text) %in% c(1, n.treatments)))
    stop("Argument 'cex.text' must be a single numeric or vector of length ",
         n.treatments, ".")
  ##
  if (!(length(col.text) %in% c(1, n.treatments)))
    stop("Argument 'col.text' must be a character string or vector of length ",
         n.treatments, ".")
  ##
  if (!(length(adj.x) %in% c(1, n.treatments)))
    stop("Argument 'adj.x' must be a single numeric or vector of length ",
         n.treatments, ".")
  ##
  if (!(length(adj.y) %in% c(1, n.treatments)))
    stop("Argument 'adj.y' must be a single numeric or vector of length ",
         n.treatments, ".")
  ##
  if (!(length(offset.x) %in% c(1, n.treatments)))
    stop("Argument 'offset.x' must be a single numeric or vector of length ",
         n.treatments, ".")
  ##
  if (!(length(offset.y) %in% c(1, n.treatments)))
    stop("Argument 'offset.y' must be a single numeric or vector of length ",
         n.treatments, ".")
  ##
  if (use_pch) {
    if (!(length(pch) %in% c(1, n.treatments)))
      stop("Argument 'pch' must be a single numeric or vector of length ",
           n.treatments, ".")
    ##
    if (!(length(cex.points) %in% c(1, n.treatments)))
      stop("Argument 'cex.points' must be a single numeric or vector of length ",
           n.treatments, ".")
    ##
    if (!(length(col.points) %in% c(1, n.treatments)))
      stop("Argument 'col.points' must be a character string or vector of length ",
           n.treatments, ".")
  }
  
  
  if (length(cex.text) == 1)
    cex.text <- rep(cex.text, n.treatments)
  ##
  if (length(col.text) == 1)
    col.text <- rep(col.text, n.treatments)
  ##
  if (length(adj.x) == 1)
    adj.x <- rep(adj.x, n.treatments)
  ##
  if (length(adj.y) == 1)
    adj.y <- rep(adj.y, n.treatments)
  ##
  if (length(offset.x) == 1)
    offset.x <- rep(offset.x, n.treatments)
  ##
  if (length(offset.y) == 1)
    offset.y <- rep(offset.y, n.treatments)
  ##
  if (use_pch) {
    if (length(pch) == 1)
      pch <- rep(pch, n.treatments)
    ##
    if (length(cex.points) == 1)
      cex.points <- rep(cex.points, n.treatments)
    ##
    if (length(col.points) == 1)
      col.points <- rep(col.points, n.treatments)
  }
  
  
  seq.treats <- seq_len(n.treatments)
  ##
  if (is_2d) {
    ##
    ## 2-D plot
    ##
    xvals <- p.matrix[, sel.x]
    yvals <- p.matrix[, sel.y]
    ##
    if (is_biplot) {
      plot(xvals, yvals,
           type = "n",
           xlab = outcomes[sel.x], ylab = outcomes[sel.y],
           ...)
    }
    else {
      ##
      plot(xvals, yvals,
           type = "n",
           xlab = outcomes[sel.x], ylab = outcomes[sel.y],
           xlim = c(0, 1), ylim = c(0, 1),
           ...)
    }
    
    if (grid & !is_biplot) {
      for (i in seq.treats) {
        lines(x = c(xvals[i], xvals[i]),
              y = c(0, yvals[i]),
              col = col.grid, lty = lty.grid, lwd = lwd.grid)
        ##
        lines(x = c(0, xvals[i]),
              y = c(yvals[i], yvals[i]),
              col = col.grid, lty = lty.grid, lwd = lwd.grid)
      }
    }
    ##
    if (use_pch)
      for (i in seq.treats)
        points(xvals[i], yvals[i],
               pch = pch[i], col = col.points[i], cex = cex.points[i])
    ##
    for (i in seq.treats) {
      for (j in seq.treats) {
        if (M0[i, j] == 1) {
          if (arrows)
            arrows(xvals[i], yvals[i],
                   xvals[j], yvals[j],
                   col = col.lines, lty = lty.lines, lwd = lwd.lines,
                   length = length)
          else
            lines(c(xvals[i], xvals[j]),
                  c(yvals[i], yvals[j]),
                   col = col.lines, lty = lty.lines, lwd = lwd.lines)
        }
      }
    }
    ##
    for (i in seq.treats)
      text(xvals[i] + offset.x[i],
           yvals[i] + offset.y[i],
           treatments[i],
           adj = c(adj.x[i], adj.y[i]),
           col = col.text[i], cex = cex.text[i])
  }
  else {
    ##
    ## 3-D plot
    ##
    xvals <- p.matrix[, sel.x]
    yvals <- p.matrix[, sel.y]
    zvals <- p.matrix[, sel.z]
    ##
    rgl::plot3d(xvals, yvals, zvals,
                xlab = outcomes[sel.x],
                ylab = outcomes[sel.y],
                zlab = outcomes[sel.z])
    ##
    if (use_pch)
      for (i in seq.treats)
        rgl::points3d(xvals[i], yvals[i], zvals[i],
                      pch = pch[i], col = col.points[i], cex = cex.points[i])
    ##
    for (i in seq.treats)
      for (j in seq.treats)
        if (M0[i, j] == 1)
          rgl::lines3d(x = c(xvals[i], xvals[j]),
                       y = c(yvals[i], yvals[j]),
                       z = c(zvals[i], zvals[j]),
                       lwd = 2)
    ##
    for (i in seq.treats)
      rgl::text3d(xvals[i], yvals[i], zvals[i],
                  treatments[i],
                  col = col.text[i], cex = cex.text[i])
  }
  
  invisible(NULL)
}
