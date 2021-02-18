#' Network graph
#' 
#' @description
#' This function generates a graph of the evidence network.
#' 
#' @param x An object of class \code{netmeta} (mandatory).
#' @param seq A character or numerical vector specifying the sequence
#'   of treatments arrangement (anticlockwise if \code{start.layout =
#'   "circle"}).
#' @param labels An optional vector with treatment labels.
#' @param cex The magnification to be used for treatment labels.
#' @param col A single color (or vector of colors) for lines
#'   connecting treatments (edges) if argument \code{plastic =
#'   FALSE}. Length of the vector must be equal to the number of edges
#'   (see list element 'comparisons' in \code{\link{netmeta}}).
#' @param adj One, two, or three values in [0, 1] (or a vector /
#'   matrix with length / number of rows equal to the number of
#'   treatments) specifying the x (and optionally y and z) adjustment
#'   for treatment labels.
#' @param offset Distance between edges (i.e. treatments) in graph and
#'   treatment labels for 2-D plots (value of 0.0175 corresponds to a
#'   difference of 1.75\% of the range on x- and y-axis).
#' @param srt.labels The character string \code{"orthogonal"} (can be
#'   abbreviated), a single numeric or numerical vector with value(s)
#'   between -180 and 180 specifying the angle to rotate treatment
#'   labels (see Details).
#' @param scale Additional space added outside of edges
#'   (i.e. treatments).  Increase this value for larger treatment
#'   labels (value of 1.10 corresponds to an additional space of 10\%
#'   around the network graph).
#' @param plastic A logical indicating whether the appearance of the
#'   comparisons should be in '3D look' (not to be confused with
#'   argument \code{dim}).
#' @param thickness Either a character variable to determine the
#'   method to plot line widths (see Details) or a matrix of the same
#'   dimension and row and column names as argument \code{A.matrix}
#'   with information on line width.
#' @param lwd A numeric for scaling the line width of comparisons.
#' @param lwd.max Maximum line width in network graph. The connection
#'   with the largest value according to argument \code{thickness}
#'   will be set to this value.
#' @param lwd.min Minimum line width in network graph. All connections
#'   with line widths below this values will be set to \code{lwd.min}.
#' @param dim A character string indicating whether a 2- or
#'   3-dimensional plot should be produced, either \code{"2d"} or
#'   \code{"3d"}.
#' @param highlight A character vector identifying comparisons that
#'   should be marked in the network graph, e.g. \code{highlight =
#'   "treat1:treat2"}.
#' @param col.highlight Color(s) to highlight the comparisons given by
#'   \code{highlight}.
#' @param scale.highlight Scaling factor(s) for the line width(s) to
#'   highlight the comparisons given by \code{highlight}.
#' @param multiarm A logical indicating whether multi-arm studies
#'   should be marked in plot.
#' @param col.multiarm Either a function from R package colorspace or
#'   grDevice to define colors for multi-arm studies or a character
#'   vector with colors to highlight multi-arm studies.
#' @param alpha.transparency The alpha transparency of colors used to
#'   highlight multi-arm studies (0 means transparent and 1 means
#'   opaque).
#' @param points A logical indicating whether points should be printed
#'   at nodes (i.e. treatments) of the network graph.
#' @param col.points,cex.points,pch.points,bg.points Corresponding
#'   color, size, type, and background color for points. Can be a
#'   vector with length equal to the number of treatments.
#' @param number.of.studies A logical indicating whether number of
#'   studies should be added to network graph.
#' @param cex.number.of.studies The magnification to be used for
#'   number of studies.
#' @param col.number.of.studies Color for number of studies.
#' @param bg.number.of.studies Color for shadow around number of
#'   studies.
#' @param pos.number.of.studies A single value (or vector of values)
#'   in [0, 1] specifying the position of the number of studies on the
#'   lines connecting treatments (edges). Length of the vector must be
#'   equal to the number of edges.
#' @param start.layout A character string indicating which starting
#'   layout is used if \code{iterate = TRUE}. If "circle" (default),
#'   the iteration starts with a circular ordering of the vertices; if
#'   "eigen", eigenvectors of the Laplacian matrix are used,
#'   calculated via generic function \code{\link{eigen}} (spectral
#'   decomposition); if "prcomp", eigenvectors of the Laplacian matrix
#'   are calculated via generic function \code{\link{prcomp}}
#'   (principal component analysis); if "random", a random layout is
#'   used, drawn from a bivariate normal.
#' @param eig1 A numeric indicating which eigenvector is used as x
#'   coordinate if \code{start = "eigen"} or \code{"prcomp"} and
#'   \code{iterate = TRUE}.  Default is 2, the eigenvector to the
#'   second-smallest eigenvalue of the Laplacian matrix.
#' @param eig2 A numeric indicating which eigenvector is used as
#'   y-coordinate if \code{start = "eigen"} or \code{"prcomp"} and
#'   \code{iterate = TRUE}.  Default is 3, the eigenvector to the
#'   third-smallest eigenvalue of the Laplacian matrix.
#' @param eig3 A numeric indicating which eigenvector is used as
#'   z-coordinate if \code{start = "eigen"} or \code{"prcomp"} and
#'   \code{iterate = TRUE}.  Default is 4, the eigenvector to the
#'   fourth-smallest eigenvalue of the Laplacian matrix.
#' @param iterate A logical indicating whether the stress majorization
#'   algorithm is carried out for optimization of the layout.
#' @param tol A numeric for the tolerance for convergence if
#'   \code{iterate = TRUE}.
#' @param maxit An integer defining the maximum number of iteration
#'   steps if \code{iterate = TRUE}.
#' @param allfigures A logical indicating whether all iteration steps
#'   are shown if \code{iterate = TRUE}. May slow down calculations if
#'   set to \code{TRUE} (especially if \code{plastic = TRUE}).
#' @param A.matrix Adjacency matrix (\emph{n}x\emph{n}) characterizing
#'   the structure of the network graph. Row and column names must be
#'   the same set of values as provided by argument \code{seq}.
#' @param N.matrix Neighborhood matrix (\emph{n}x\emph{n}) replacing
#'   A.matrix if neighborhood is to be specified differently from node
#'   adjacency in the network graph, for example content-based. Row
#'   and column names must be the same set of values as provided by
#'   argument \code{seq}.
#' @param D.matrix Distance matrix (\emph{n}x\emph{n}) replacing
#'   A.matrix and N.matrix if distances should be provided
#'   directly. Row and column names must be the same set of values as
#'   provided by argument \code{seq}.
#' @param xpos Vector (\emph{n}) of x coordinates.
#' @param ypos Vector (\emph{n}) of y coordinates.
#' @param zpos Vector (\emph{n}) of z coordinates.
#' @param figure A logical indicating whether network graph should be
#'   shown.
#' @param \dots Additional graphical arguments.
#' 
#' @details
#' The network is laid out in the plane, where the nodes in the graph
#' layout correspond to the treatments and edges display the observed
#' treatment comparisons. For the default setting, nodes are placed on
#' a circle.  Other starting layouts are "eigen", "prcomp", and
#' "random" (Rücker & Schwarzer 2015). If \code{iterate = TRUE}, the
#' layout is further optimized using the stress majorization
#' algorithm. This algorithm specifies an 'ideal' distance (e.g., the
#' graph distance) between two nodes in the plane. In the optimal
#' layout, these distances are best approximated in the sense of least
#' squares.  Starting from an initial layout, the optimum is
#' approximated in an iterative process called stress majorization
#' (Kamada and Kawai 1989, Michailidis and de Leeuw 2001, Hu
#' 2012). The starting layout can be chosen as a circle or coming from
#' eigenvectors of the Laplacian matrix (corresponding to Hall's
#' algorithm, Hall 1970), calculated in different ways, or
#' random. Moreover, it can be chosen whether the iteration steps are
#' shown (argument \code{allfigures = TRUE}).
#' 
#' An optimized circular presentation which typically has a reduced
#' (sometimes minimal) number of crossings can be achieved by using
#' argument \code{seq = "optimal"} in combination with argument
#' \code{start.layout}. Note, is is not possible of prespecify the
#' best value for argument \code{start.layout} for any situation as
#' the result depends on the network structure.
#' 
#' Argument \code{thickness} providing the line width of the nodes
#' (comparisons) can be a matrix of the same dimension as argument
#' \code{A.matrix} or any of the following character variables:
#' \itemize{
#' \item Same line width (argument \code{lwd}) for all comparisons
#'   (\code{thickness = "equal"})
#' \item Proportional to number of studies comparing two treatments
#'   (\code{thickness = "number.of.studies"})
#' \item Proportional to inverse standard error of fixed effects model
#'   comparing two treatments (\code{thickness = "se.fixed"})
#' \item Proportional to inverse standard error of random effects
#'   model comparing two treatments (\code{thickness = "se.random"})
#' \item Weight from fixed effects model comparing two treatments
#'   (\code{thickness = "w.fixed"})
#' \item Weight from random effects model comparing two treatments
#'   (\code{thickness = "w.random"})
#' }
#'
#' Only evidence from direct treatment comparisons is considered to
#' determine the line width if argument \code{thickness} is equal to
#' any but the first method. By default, \code{thickness = "se.fixed"}
#' is used if \code{start.layout = "circle"}, \code{iterate = FALSE},
#' and \code{plastic = TRUE}. Otherwise, the same line width is used.
#'
#' Argument \code{srt.labels} can be used to specific the rotation (in
#' degrees) of the treatment labels. If \code{srt.labels} is equal to
#' \code{"orthogonal"}, treatment labels are orthogonal to the
#' circle. If \code{srt.labels} is a single numeric, all labels are
#' rotated by this degree. If \code{srt.labels} is a numeric vector,
#' it must be of the same length as the number of treatments and
#' labels are rotated counter-clockwise starting on the right
#' side. Finally, if \code{srt.labels} is a named numeric vector, it
#' must be of the same length as the number of treatments and the
#' names must be equal to the treatment names (and treatment labels
#' are rotated according to the specified values).
#' 
#' Further, a couple of graphical parameters can be specified, such as
#' color and appearance of the edges (treatments) and the nodes
#' (comparisons), whether special comparisons should be highlighted
#' and whether multi-arm studies should be indicated as colored
#' polygons. By default, if R package colorspace is available the
#' \code{\link[colorspace]{sequential_hcl}} function is used to
#' highlight multi-arm studies; otherwise the \code{\link{rainbow}} is
#' used.
#' 
#' In order to generate 3-D plots (argument \code{dim = "3d"}), R
#' package \bold{rgl} is necessary. Note, under macOS the X.Org X
#' Window System must be available (see
#' \url{https://www.xquartz.org}).
#'
#' @return A data frame containing the following columns:
#' \item{labels}{Treatment labels.}
#' \item{seq}{Sequence of treatment labels.}
#' \item{xpos}{Position of treatment / edge on x-axis.}
#' \item{ypos}{Position of treatment / edge on y-axis.}
#' \item{zpos}{Position of treatment / edge on z-axis (for 3-D
#'   plots).}
#' \item{xpos.labels}{Position of treatment labels on x-axis (for 2-D
#'   plots).}
#' \item{ypos.labels}{Position of treatment labels on y-axis (for 2-D
#'   plots).}
#' \item{adj.x}{Adjustment for treatment label on x-axis.}
#' \item{adj.y}{Adjustment for treatment label on y-axis.}
#' \item{adj.z}{Adjustment for treatment label on z-axis (for 3-D
#'   plots).}
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Ulrike
#'   Krahn \email{ulrike.krahn@@bayer.com}, Jochem König
#'   \email{koenigjo@@uni-mainz.de}, Guido Schwarzer
#'   \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}
#' 
#' @references
#'
#' Hall KM (1970):
#' An r-dimensional quadratic placement algorithm.
#' \emph{Management Science},
#' \bold{17}, 219--29
#' 
#' Hu Y (2012):
#' \emph{Combinatorial Scientific Computing}, Chapter Algorithms for
#' Visualizing Large Networks, pages 525--49.
#' Chapman and Hall / CRC,  Computational Science.
#' 
#' Kamada T, Kawai S (1989):
#' An algorithm for drawing general undirected graphs.
#' \emph{Information Processing Letters},
#' \bold{31}, 7--15
#' 
#' Krahn U, Binder H, König J (2013):
#' A graphical tool for locating inconsistency in network meta-analyses.
#' \emph{BMC Medical Research Methodology},
#' \bold{13}, 35
#' 
#' Michailidis G, de Leeuw J (2001):
#' Data visualization through graph drawing.
#' \emph{Computational Statistics},
#' \bold{16}, 435--50
#' 
#' Rücker G, Schwarzer G (2016):
#' Automated drawing of network plots in network meta-analysis.
#' \emph{Research Synthesis Methods},
#' \bold{7}, 94--107
#' 
#' @keywords hplot
#'
#' @examples
#' data(Senn2013)
#' 
#' # Generation of an object of class 'netmeta' with reference
#' # treatment 'plac'
#' #
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'                 data = Senn2013, sm = "MD", reference = "plac")
#' 
#' # Network graph with default settings
#' #
#' netgraph(net1)
#' 
#' \dontrun{
#' # Network graph with specified order of the treatments and one
#' # highlighted comparison
#' #
#' trts <- c("plac", "benf", "migl", "acar", "sulf",
#'           "metf", "rosi", "piog", "sita", "vild")
#' netgraph(net1, highlight = "rosi:plac", seq = trts)
#' 
#' # Same network graph using argument 'seq' in netmeta function
#' #
#' net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'                 data = Senn2013, sm = "MD", reference = "plac",
#'                 seq = trts)
#' netgraph(net2, highlight = "rosi:plac")
#' 
#' # Network graph optimized, starting from a circle, with multi-arm
#' # study colored
#' #
#' netgraph(net1, start = "circle", iterate = TRUE, col.multiarm = "purple")
#'
#' # Network graph optimized, starting from a circle, with multi-arm
#' # study colored and all intermediate iteration steps visible
#' #
#' netgraph(net1, start = "circle", iterate = TRUE, col.multiarm = "purple",
#'          allfigures = TRUE)
#' 
#' # Network graph optimized, starting from Laplacian eigenvectors,
#' # with multi-arm study colored
#' #
#' netgraph(net1, start = "eigen", col.multiarm = "purple")
#' 
#' # Network graph optimized, starting from different Laplacian
#' # eigenvectors, with multi-arm study colored
#' #
#' netgraph(net1, start = "prcomp", col.multiarm = "purple")
#' 
#' # Network graph optimized, starting from random initial layout,
#' # with multi-arm study colored
#' #
#' netgraph(net1, start = "random", col.multiarm = "purple")
#' 
#' # Network graph without plastic look and one highlighted comparison
#' #
#' netgraph(net1, plastic = FALSE, highlight = "rosi:plac")
#' 
#' # Network graph without plastic look and comparisons with same
#' # thickness
#' #
#' netgraph(net1, plastic = FALSE, thickness = FALSE)
#' 
#' # Network graph with changed labels and specified order of the
#' # treatments
#' #
#' netgraph(net1, seq = c(1, 3, 5, 2, 9, 4, 7, 6, 8, 10),
#'          labels = LETTERS[1:10])
#' 
#' # Rotate treatment labels (orthogonal to circle)
#' #
#' netgraph(net1, srt.labels = "o")
#' 
#' # Network graph in 3-D (opens a new device, where you may rotate and
#' # zoom the plot using the mouse / the mouse wheel).
#' # The rgl package must be installed for 3-D plots.
#' #
#' netgraph(net1, dim = "3d")
#' }
#' 
#' @method netgraph netmeta
#' @export
#' @export netgraph.netmeta


netgraph.netmeta <- function(x, seq = x$seq,
                             labels = x$trts,
                             cex = 1, adj = NULL, srt.labels = 0,
                             offset = if (!is.null(adj) && all(unique(adj) == 0.5)) 0 else 0.0175,
                             scale = 1.10,
                             col = "slateblue", plastic, thickness,
                             lwd = 5, lwd.min = lwd / 2.5, lwd.max = lwd * 4,
                             dim = "2d",
                             ##
                             highlight = NULL, col.highlight = "red2",
                             scale.highlight = 1,
                             ##
                             multiarm = any(x$narms > 2),
                             col.multiarm = NULL,
                             alpha.transparency = 0.5,
                             ##
                             points = FALSE, col.points = "red",
                             cex.points = 1, pch.points = 20,
                             bg.points = "gray",
                             ##
                             number.of.studies = FALSE,
                             cex.number.of.studies = cex,
                             col.number.of.studies = "white",
                             bg.number.of.studies = "black",
                             pos.number.of.studies = 0.5,
                             ##
                             start.layout = ifelse(dim == "2d", "circle", "eigen"),
                             eig1 = 2, eig2 = 3, eig3 = 4,
                             iterate,
                             tol = 0.0001, maxit = 500, allfigures = FALSE,
                             A.matrix = x$A.matrix,
                             N.matrix = sign(A.matrix),
                             D.matrix = netdistance(N.matrix),
                             ##
                             xpos = NULL, ypos = NULL, zpos = NULL,
                             figure = TRUE,
                             ...) {
  
  
  meta:::chkclass(x, "netmeta")
  ##
  x <- upgradenetmeta(x)
  ##
  n.edges <- sum(x$A.matrix[upper.tri(x$A.matrix)] > 0)
  n.trts <- length(x$trts)


  setchar <- meta:::setchar
  chknumeric <- meta:::chknumeric
  chklogical <- meta:::chklogical
  ##
  if (length(srt.labels) == 1 && is.character(srt.labels)) {
    srt.labels <- setchar(srt.labels, "orthogonal",
                          "should be equal to 'orthogonal' or numeric (vector)")
    if (!missing(iterate) && iterate == TRUE) {
      warning("Orthogonal labels not supported if argument 'iterate = TRUE'")
      srt.labels <- 0
    }
    ##
    srtfunc <- function(ntrt) {
      s <- 180 * (2 * (1:ntrt) / ntrt - 1)
      for (i in 1:ntrt) {
        if (i < ntrt / 4)
          s[i] <- 360 * i / ntrt
        if (i > 3 * ntrt / 4)
          s[i] <- 360 * (i / ntrt - 1)
      }
      s
    }
    srt.labels <- srtfunc(x$n)
  }
  chknumeric(srt.labels, min = -180, max = 180)
  chknumeric(lwd, min = 0, zero = TRUE, length = 1)
  chknumeric(lwd.min, min = 0, zero = TRUE, length = 1)
  chknumeric(lwd.max, min = 0, zero = TRUE, length = 1)
  chknumeric(scale.highlight, min = 0, zero = TRUE)
  ##
  if (lwd.min > lwd.max)
    stop("Argument 'lwd.min' must be smaller than 'lwd.max'.")
  ##
  dim <- setchar(dim, c("2d", "3d"))
  is_2d <- dim == "2d"
  is_3d <- !is_2d
  ##
  if (is_3d & !meta:::is.installed.package("rgl", stop = FALSE)) {
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
  missing.start.layout <- missing(start.layout)
  start.layout <-
    setchar(start.layout, c("eigen", "prcomp", "circle", "random"))
  ##
  mf <- match.call()
  ##
  if (!missing(seq) & is.null(seq))
    stop("Argument 'seq' must be not NULL.")
  ##
  if (!missing(labels)) {
    ##
    labels <- eval(mf[[match("labels", names(mf))]],
                   x, enclos = sys.frame(sys.parent()))
    ##
    if (is.null(labels))
      stop("Argument 'labels' must be not NULL.")
  }
  ##
  ## Colors of edges
  ##
  if (is.matrix(col)) {
    if ((dim(col)[1] != dim(A.matrix)[1]) |
        (dim(col)[2] != dim(A.matrix)[2]))
      stop("Dimension of argument 'A.matrix' and 'col' are different.")
    if (is.null(dimnames(col)))
      stop("Matrix 'col' must have row and column names identical ",
           "to argument 'A.matrix'.")
    else {
      if (any(rownames(col) != rownames(A.matrix)))
        stop("Row names of matrix 'col' must be identical to ",
             "argument 'A.matrix'.")
      if (any(colnames(col) != colnames(A.matrix)))
        stop("Column names of matrix 'col' must be identical to ",
             "argument 'A.matrix'.")
    }
    ##
    col <- col[lower.tri(col)]
    col <- col[!is.na(col)]
  }
  ##
  n.col <- length(col)
  ##
  if (n.col == 1)
    col <- rep(col, n.edges)
  else if (n.col != n.edges)
    stop("Length of argument 'col' (",
         n.col, ") is different from the number of direct pairwise comparisons (",
         n.edges, ")")
  ##
  n.pos <- length(pos.number.of.studies)
  if (n.pos == 1)
    pos.number.of.studies <- rep(pos.number.of.studies, n.edges)
  else if (n.pos != n.edges)
    stop("Length of argument 'pos.number.of.studies' (",
         n.pos, ") is different from the number of direct pairwise comparisons (",
         n.edges, ")")
  ##
  if (!missing(cex.points))
    cex.points <- eval(mf[[match("cex.points", names(mf))]],
                       x, enclos = sys.frame(sys.parent()))
  ##
  if (length(cex.points) == 1)
    cex.points <- rep(cex.points, n.trts)
  else if (length(cex.points) != n.trts)
    stop("Length of argument 'cex.points' must be equal to the number of treatments.")
  ##
  if (length(col.points) == 1)
    col.points <- rep(col.points, n.trts)
  else if (length(col.points) != n.trts)
    stop("Length of argument 'col.points' must be equal to the number of treatments.")
  ##
  if (length(pch.points) == 1)
    pch.points <- rep(pch.points, n.trts)
  else if (length(pch.points) != n.trts)
    stop("Length of argument 'pch.points' must be equal to number of treatments.")
  ##
  if (length(bg.points) == 1)
    bg.points <- rep(bg.points, n.trts)
  else if (length(bg.points) != n.trts)
    stop("Length of argument 'bg.points' must be equal to the number of treatments.")
  ##
  chklogical(figure)
  ##
  if (is.null(seq) | (length(seq) == 1 & x$d == 1))
    seq1 <- 1:length(labels)
  else if (length(seq) == 1 & x$d > 1) {
    seq <- setchar(seq, "optimal", "should be equal to 'optimal' or a permutation of treatments")
    ##
    if (missing.start.layout)
      start.layout <- "eigen"
    ##
    seq1 <- optcircle(x, start.layout = start.layout)$seq
    ##
    start.layout <- "circle"
  }
  else if (!(start.layout == "circle" & (missing(iterate) || iterate == FALSE))) {
    seq1 <- 1:length(labels)
    ##
    if (!missing(seq) & !is.null(seq) & (is.null(xpos) & is.null(ypos)))
      warning("Argument 'seq' only considered if start.layout=\"circle\" and iterate=FALSE.")
  }
  else
    seq1 <- charmatch(setseq(seq, x$trts), x$trts)
  
  
  if (missing(iterate))
    iterate <- ifelse(start.layout == "circle", FALSE, TRUE)
  else if (length(seq) == 1 && seq == "optimal") {
    warning("Argument 'iterate' ignored as argument 'seq' is equal to \"optimal\".")
    iterate <- FALSE
  }
  ##
  if (!missing(allfigures) && length(seq) == 1 && seq == "optimal") {
    warning("Argument 'allfigures' ignored as argument 'seq' is equal to \"optimal\".")
    allfigures <- FALSE
  }
  
  
  addargs <- names(list(...))
  ##
  if ("highlight.split" %in% addargs)
    warning("Argument 'highlight.split' has been removed from ",
            "R function netgraph.\n  This argument has been replaced by ",
            "argument 'sep.trts' in R function netmeta.")
  ##
  highlight.split <- x$sep.trts
  ##
  if (is.null(highlight.split))
    highlight.split <- ":"
  ##
  n.high <- 1
  if (!is.null(highlight)) {
    n.high <- length(highlight)
    ##
    if (!missing(col.highlight))
      if (length(col.highlight) != 1 && length(col.highlight) != n.high)
        stop("Argument 'col.highlight' must be a single value or ",
             "of same length as argument 'highlight'.", call. = FALSE)
    ##
    if (!missing(scale.highlight))
      if (length(scale.highlight) != 1 && length(scale.highlight) != n.high)
        stop("Argument 'scale.highlight' must be a single value or ",
             "of same length as argument 'highlight'.", call. = FALSE)
  }
  ##
  if (length(col.highlight) == 1 & n.high > 1)
    col.highlight <- rep(col.highlight, n.high)
  
  
  if (missing(plastic))
    if (start.layout == "circle" & iterate == FALSE & is_2d)
      plastic <- TRUE
    else
      plastic <- FALSE


  if (missing(thickness)) {
    if (start.layout == "circle" & iterate == FALSE & plastic == TRUE) {
      thick <- "se.fixed"
      thickness <- "se.fixed"
    }
    else {
      thick <- "equal"
      thickness <- "equal"
    }
  }
  else {
    if (!is.matrix(thickness)) {
      if (length(thickness) == 1 & is.character(thickness))
        thick <- setchar(thickness,
                         c("equal", "number.of.studies",
                           "se.fixed", "se.random", "w.fixed", "w.random"))
      ##
      else if (length(thickness) == 1 & is.logical(thickness)) {
        if (thickness)
          thick <- "se.fixed"
        else
          thick <- "equal"
      }
    }
    else {
      if ((dim(thickness)[1] != dim(A.matrix)[1]) |
          (dim(thickness)[2] != dim(A.matrix)[2]))
        stop("Dimension of argument 'A.matrix' and 'thickness' are different.")
      if (is.null(dimnames(thickness)))
        stop("Matrix 'thickness' must have row and column names identical to ",
             "argument 'A.matrix'.")
      else {
        if (any(rownames(thickness) != rownames(A.matrix)))
          stop("Row names of matrix 'thickness' must be identical to ",
               "argument 'A.matrix'.")
        if (any(colnames(thickness) != colnames(A.matrix)))
          stop("Column names of matrix 'thickness' must be identical to ",
               "argument 'A.matrix'.")
      }
      ##
      W.matrix <- thickness
      thick <- "matrix"
    }
  }


  if (allfigures & is_3d) {
    warning("Argument 'allfigures' set to FALSE for 3-D network plot.")
    allfigures <- FALSE
  }


  col.matrix <- matrix("", nrow = n.trts, ncol = n.trts)
  dimnames(col.matrix) <- dimnames(A.matrix)
  col.matrix <- t(col.matrix)
  col.matrix[lower.tri(col.matrix) & t(A.matrix) > 0] <- col
  tcm <- col.matrix
  col.matrix <- t(col.matrix)
  col.matrix[lower.tri(col.matrix)] <- tcm[lower.tri(tcm)]
  ##
  pos.matrix <- matrix(NA, nrow = n.trts, ncol = n.trts)
  dimnames(pos.matrix) <- dimnames(A.matrix)
  pos.matrix <- t(pos.matrix)
  pos.matrix[lower.tri(pos.matrix) & t(A.matrix) > 0] <- pos.number.of.studies
  tam <- pos.matrix
  pos.matrix <- t(pos.matrix)
  pos.matrix[lower.tri(pos.matrix)] <- tam[lower.tri(tam)]
  ##
  A.matrix <- A.matrix[seq1, seq1]
  N.matrix <- N.matrix[seq1, seq1]
  D.matrix <- D.matrix[seq1, seq1]
  ##
  col.matrix <- col.matrix[seq1, seq1]
  pos.matrix <- pos.matrix[seq1, seq1]
  ##
  if (thick == "matrix")
    W.matrix <- W.matrix[seq1, seq1]
  ##
  trts <- x$trts[seq1]
  labels <- labels[seq1]
  ##
  col.points <- col.points[seq1]
  cex.points <- cex.points[seq1]
  pch.points <- pch.points[seq1]
  bg.points  <- bg.points[seq1]
  
  A.sign <- sign(A.matrix)
  
  if ((is_2d & (is.null(xpos) & is.null(ypos))) |
      (is_3d & (is.null(xpos) & is.null(ypos) & is.null(zpos)))) {
    stressdata <- stress(x,
                         A.matrix = A.matrix,
                         N.matrix = N.matrix,
                         D.matrix = D.matrix,
                         ##
                         dim = dim,
                         start.layout = start.layout,
                         iterate = iterate,
                         eig1 = eig1, eig2 = eig2, eig3 = eig3,
                         tol = tol,
                         maxit = maxit,
                         ##
                         allfigures = allfigures,
                         ##
                         seq = seq,
                         ##
                         labels = labels,
                         cex = cex,
                         col = col,
                         adj = adj,
                         srt.labels = srt.labels,
                         offset = offset,
                         scale = scale,
                         ##
                         plastic = plastic,
                         thickness = thickness,
                         lwd = lwd,
                         lwd.min = lwd.min,
                         lwd.max = lwd.max,
                         ##
                         highlight = highlight,
                         col.highlight = col.highlight,
                         scale.highlight = scale.highlight,
                         ## multiarm
                         col.multiarm = col.multiarm,
                         alpha.transparency = alpha.transparency,
                         ##
                         points = points, col.points = col.points,
                         cex.points = cex.points, pch.points = pch.points,
                         bg.points = bg.points,
                         ##
                         number.of.studies = number.of.studies,
                         cex.number.of.studies = cex.number.of.studies,
                         col.number.of.studies = col.number.of.studies,
                         bg.number.of.studies = bg.number.of.studies,
                         pos.number.of.studies = pos.number.of.studies,
                         ##
                         ...)
    ##
    xpos <- stressdata$x
    ypos <- stressdata$y
    if (is_3d)
      zpos <- stressdata$z
  }


  if (allfigures)
    return(invisible(NULL))


  n <- dim(A.matrix)[1]
  d <- scale * max(abs(c(min(c(xpos, ypos), na.rm = TRUE),
                         max(c(xpos, ypos), na.rm = TRUE))))


  ##
  ##
  ## Generate datasets for plotting
  ##
  ##
  ##
  ## Dataset for nodes
  ##
  dat.nodes <- data.frame(trts, labels, seq, srt = NA,
                          xpos, ypos, zpos = NA,
                          xpos.labels = NA, ypos.labels = NA,
                          cex = cex.points,
                          col = col.points,
                          pch = pch.points,
                          bg = bg.points,
                          stringsAsFactors = FALSE)
  if (is_2d)
    dat.nodes$zpos <- NULL
  else {
    dat.nodes$zpos <- zpos
    dat.nodes$zpos.labels <- NA
  }
  ##
  if (is.null(adj)) {
    dat.nodes$adj.x <- NA
    dat.nodes$adj.y <- NA
    if (!is_2d)
      dat.nodes$adj.z <- NA
    ##
    dat.nodes$adj.x[dat.nodes$xpos >= 0] <- 0
    dat.nodes$adj.x[dat.nodes$xpos <  0] <- 1
    ##
    dat.nodes$adj.y[dat.nodes$ypos >  0] <- 0
    dat.nodes$adj.y[dat.nodes$ypos <= 0] <- 1
    ##
    if (!is_2d) {
      dat.nodes$adj.z[dat.nodes$zpos >  0] <- 0
      dat.nodes$adj.z[dat.nodes$zpos <= 0] <- 1
    }
  }
  else {
    dat.nodes$adj.x <- NA
    dat.nodes$adj.y <- NA
    ##
    if (length(adj) == 1) {
      dat.nodes$adj.x <- adj
      dat.nodes$adj.y <- adj
      if (!is_2d)
        dat.nodes$adj.z <- adj
    }
    else if (length(adj) == 2) {
      dat.nodes$adj.x <- adj[1]
      dat.nodes$adj.y <- adj[2]
      if (!is_2d)
        dat.nodes$adj.z <- 0.5
    }
    else if (length(adj) == 3 & !is_2d) {
      dat.nodes$adj.x <- adj[1]
      dat.nodes$adj.y <- adj[2]
      dat.nodes$adj.z <- adj[3]
    }
    else if (is.vector(adj)) {
      if (length(adj) != length(labels))
        stop("Length of vector 'adj' must be equal to number of treatments.")
      ##
      names(adj) <- x$trts
      dat.nodes$adj.x <- adj[seq1]
      dat.nodes$adj.y <- adj[seq1]
      ##
      if (!is_2d)
        dat.nodes$adj.z <- adj[seq1]
    }
    else if (is.matrix(adj)) {
      if (nrow(adj) != length(labels))
        stop("Number of rows of matrix 'adj' must be equal to number of treatments.")
      rownames(adj) <- x$trts
      dat.nodes$adj.x <- adj[seq1, 1]
      dat.nodes$adj.y <- adj[seq1, 2]
      ##
      if (!is_2d & ncol(adj) >= 3)
        dat.nodes$adj.z <- adj[seq1, 3]
    }
  }
  ##
  if (is_2d) {
    offset <- offset * 2 * d
    ##
    if (length(offset) == 1) {
      offset.x <- offset
      offset.y <- offset
    }
    else if (length(offset) == 2) {
      offset.x <- offset[1]
      offset.y <- offset[2]
    }
    else if (is.vector(offset)) {
      if (length(offset) != length(labels))
        stop("Length of vector 'offset' must be equal to number of treatments.")
      ##
      names(offset) <- x$trts
      offset.x <- offset[seq1]
      offset.y <- offset[seq1]
    }
    else if (is.matrix(adj)) {
      if (nrow(offset) != length(labels))
        stop("Number of rows of matrix 'offset' must be equal to number of treatments.")
      ##
      rownames(offset) <- x$trts
      offset.x <- offset[seq1, 1]
      offset.y <- offset[seq1, 2]
    }
    ##
    dat.nodes$xpos.labels <- dat.nodes$xpos - offset.x +
      2 * (dat.nodes$adj.x == 0) * offset.x
    dat.nodes$ypos.labels <- dat.nodes$ypos - offset.y +
      2 * (dat.nodes$adj.y == 0) * offset.y
  }
  else {
    dat.nodes$xpos.labels <- dat.nodes$xpos
    dat.nodes$ypos.labels <- dat.nodes$ypos
    dat.nodes$zpos.labels <- dat.nodes$zpos
  }
  ##
  if (length(srt.labels) == 1)
    dat.nodes$srt <- srt.labels
  else {
    if (length(srt.labels) != length(labels))
      stop("Length of vector 'srt.labels' must be equal to",
           "number of treatments.",
           eval. = FALSE)
    if (is.null(names(srt.labels)))
      dat.nodes$srt <- srt.labels
    else {
      ## Check names of named vector 'srt.labels'
      names.srt.labels <- names(srt.labels)
      names.srt.labels <- setseq(names.srt.labels, dat.nodes$trts,
                                 paste0("Names of vector provided in ",
                                        "argument 'srt.labels'"))
      dat.nodes$srt <- srt.labels[dat.nodes$trts]
    }
  }
  ##
  ## Dataset for edges
  ##
  dat.edges <- data.frame(treat1 = rep("", n.edges),
                          treat2 = "",
                          n.stud = NA,
                          xpos = NA, ypos = NA,
                          adj = NA, pos.number.of.studies,
                          col = "",
                          stringsAsFactors = FALSE)
  ##
  comp.i <- 1
  ##
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (A.sign[i, j] > 0) {
        ##
        dat.edges$treat1[comp.i] <- rownames(A.matrix)[i]
        dat.edges$treat2[comp.i] <- colnames(A.matrix)[j]
        dat.edges$n.stud[comp.i] <- A.matrix[i, j]
        dat.edges$adj[comp.i] <- lambda <- pos.matrix[i, j]
        dat.edges$xpos[comp.i] <- lambda * xpos[i] + (1 - lambda) * xpos[j]
        dat.edges$ypos[comp.i] <- lambda * ypos[i] + (1 - lambda) * ypos[j]
        dat.edges$col[comp.i] <- col.matrix[i, j]
        ##
        comp.i <- comp.i + 1
      }
    }
  }
  
  
  ##
  ## Define coloured regions for multi-arm studies
  ##
  if (multiarm) {
    mc <- multicols(x$studies, x$narms, missing(col.multiarm),
                    col.multiarm, alpha.transparency)
    col.polygon <- mc$cols
    multiarm.studies <- mc$multiarm.studies
    n.multi <- length(multiarm.studies)
  }
  ##
  ## Define line width
  ##
  if (thick == "number.of.studies") {
    W.matrix <- lwd.max * A.matrix / max(A.matrix)
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "equal") {
    W.matrix <- lwd * A.sign
  }
  else if (thick == "se.fixed") {
    IV.matrix <- x$seTE.direct.fixed[seq1, seq1]
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max * min(IV.matrix, na.rm = TRUE) / IV.matrix
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "se.random") {
    IV.matrix <- x$seTE.direct.random[seq1, seq1]
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max * min(IV.matrix, na.rm = TRUE) / IV.matrix
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "w.fixed") {
    IV.matrix <- 1 / x$seTE.direct.fixed[seq1, seq1]^2
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max * IV.matrix / max(IV.matrix, na.rm = TRUE)
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "w.random") {
    IV.matrix <- 1 / x$seTE.direct.random[seq1, seq1]^2
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max * IV.matrix / max(IV.matrix, na.rm = TRUE)
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "matrix") {
    W.matrix[is.infinite(W.matrix)] <- NA
    if (min(W.matrix[W.matrix != 0], na.rm = TRUE) == max(W.matrix[W.matrix != 0], na.rm = TRUE))
      W.matrix <- lwd * W.matrix
    else
      W.matrix <- lwd.max * W.matrix / max(W.matrix, na.rm = TRUE)
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }





  ##
  ##
  ## Plot graph
  ##
  ##
  if (figure) {
    range <- c(-d, d)
    ##
    if (is_2d) {
      oldpar <- par(xpd = TRUE, pty = "s")
      on.exit(par(oldpar))
      ##
      plot(xpos, ypos,
           xlim = range, ylim = range,
           type = "n", axes = FALSE, bty = "n",
           xlab = "", ylab = "",
           ...)
      ##
      ## Add coloured regions for multi-arm studies
      ##
      if (multiarm) {
        ##
        if (n.multi > 0) {
          multiarm.treat <- vector("list", n.multi)
          if (length(col.polygon) == 1)
            col.polygon <- rep(col.polygon, n.multi)
          for (i in 1:n.multi) {
            treat1 <- x$treat1[x$studlab %in% multiarm.studies[i]]
            treat2 <- x$treat2[x$studlab %in% multiarm.studies[i]]
            multiarm.treat[[i]] <- sort(unique(c(treat2, treat1)))
            ##
            dat.multi <- dat.nodes[dat.nodes$trts %in% multiarm.treat[[i]], ]
            if (nrow(dat.multi) == 0)
              dat.multi <- dat.nodes[dat.nodes$trts %in% multiarm.treat[[i]], ]
            ##
            ## Clockwise ordering of polygon coordinates
            ##
            polysort <- function(x, y) {
              xnorm <- (x - mean(x)) / sd(x) # Normalise coordinate x
              ynorm <- (y - mean(y)) / sd(y) # Normalise coordinate y
              r <- sqrt(xnorm^2 + ynorm^2)   # Calculate polar coordinates
              cosphi <- xnorm / r
              sinphi <- ynorm / r
              s <- as.numeric(sinphi > 0) # Define angles to lie in [0, 2 * pi]
              phi <- acos(cosphi)
              alpha <- s * phi + (1 - s) * (2 * pi - phi)
              ##
              res <- order(alpha)
              res
            }
            ##
            dat.multi <- dat.multi[polysort(dat.multi$xpos, dat.multi$ypos), ]
            ##
            polygon(dat.multi$xpos, dat.multi$ypos,
                    col = col.polygon[i], border = NA)
          }
        }
      }
      ##
      ## Define lines
      ##
      if (plastic) {
        n.plastic <- 30
        lwd.multiply <- rep(NA, n.plastic)
        cols <- rep("", n.plastic)
        cols.highlight <- matrix("", nrow = n.high, ncol = n.plastic)
        scales.highlight <- matrix(scale.highlight,
                                   nrow = n.high, ncol = n.plastic)
        ##
        j <- 0
        for (i in n.plastic:1) {
          j <- j + 1
          lwd.multiply[j] <- sin(pi * i / 2 / n.plastic)
          cols[j] <- paste("gray", round(100 * (1 - i / n.plastic)), sep = "")
          cols.highlight[, j] <- paste("gray", round(100 * (1 - i / n.plastic)),
                                       sep = "")
        }
        ##
        for (h in seq_len(n.high)) {
          col.high.h <- col.highlight[h]
          if (col.high.h != "transparent") {
            if (nchar(col.high.h) > 1 &
                substring(col.high.h, nchar(col.high.h)) %in% 1:4)
              col.high.h <- substring(col.high.h, 1, nchar(col.high.h) - 1)
            ##
            cols.highlight[h, 1:12] <- rep(paste(col.high.h, 4:1, sep = ""),
                                           rep(3, 4))
            cols.highlight[h, 13:15] <- rep(col.high.h, 3)
          }
          else {
            cols.highlight[h, ] <- col.high.h
          }
        }
      }
      else {
        lwd.multiply <- 1
        cols <- col
        cols.highlight <- matrix(col.highlight, nrow = n.high, ncol = 1)
        scales.highlight <- matrix(scale.highlight, nrow = n.high, ncol = 1)
      }
      ##
      ## Add highlighted comparisons
      ##
      A.sign.add.lines <- A.sign
      ##
      if (!is.null(highlight)) {
        high.i <- 0
        for (high in highlight) {
          high.i <- high.i + 1
          highs <- unlist(compsplit(high, split = highlight.split))
          if (length(highs) != 2)
            stop("Wrong format for argument 'highlight' ",
                 "(see helpfile of plotgraph command).")
          ##
          if (sum(dat.nodes$trts %in% highs) != 2)
            stop(paste0("Argument 'highlight' must contain two of ",
                        "the following values ",
                        "(separated by \":\"):\n  ",
                        paste(paste("'", dat.nodes$trts, "'", sep = ""),
                              collapse = " - "), sep = ""))
          ##
          dat.high <- dat.nodes[dat.nodes$trts %in% highs, ]
          ##
          if (is_2d)
            for (n.plines in 1:length(lwd.multiply)) {
              lines(dat.high$xpos, dat.high$ypos,
                    lwd = W.matrix[trts == highs[1], trts == highs[2]] *
                      lwd.multiply[n.plines] *
                      scales.highlight[high.i, n.plines],
                    col = cols.highlight[high.i, n.plines])
            }
          ##
          A.sign.add.lines[trts == highs[1], trts == highs[2]] <- 0
          A.sign.add.lines[trts == highs[2], trts == highs[1]] <- 0
        }
      }
      ##
      ## Add lines
      ##
      comp.i <- 1
      ##
      for (n.plines in 1:length(lwd.multiply)) {
        for (i in 1:(n - 1)) {
          for (j in (i + 1):n) {
            ##
            if (plastic)
              col.ij <- cols[n.plines]
            else
              col.ij <- col.matrix[i, j]
            ##
            if (A.sign.add.lines[i, j] > 0) {
              lines(c(xpos[i], xpos[j]), c(ypos[i], ypos[j]),
                    lwd = W.matrix[i, j] * lwd.multiply[n.plines],
                    col = col.ij)
              ##
              comp.i <- comp.i + 1
            }
          }
        }
      }
      ##
      ## Add points for labels
      ##
      if (points)
        points(xpos, ypos,
               pch = pch.points, cex = cex.points, col = col.points,
               bg = bg.points)
      ##
      ## Print treatment labels
      ##
      if (!is.null(labels))
        for (i in 1:n)
          text(dat.nodes$xpos.labels[i], dat.nodes$ypos.labels[i],
               labels = dat.nodes$labels[i],
               cex = cex,
               adj = c(dat.nodes$adj.x[i], dat.nodes$adj.y[i]),
               srt = dat.nodes$srt[i])
      ##
      ## Print number of treatments
      ##
      if (number.of.studies) {
        comp.i <- 1
        ##
        for (i in 1:(n - 1)) {
          for (j in (i + 1):n) {
            if (A.sign[i, j] > 0) {
              ##
              shadowtext(dat.edges$xpos[comp.i],
                         dat.edges$ypos[comp.i],
                         labels = dat.edges$n.stud[comp.i],
                         cex = cex.number.of.studies,
                         col = col.number.of.studies,
                         bg = bg.number.of.studies)
              ##
              comp.i <- comp.i + 1
            }
          }
        }
      }
    }
    else {
      rgl::plot3d(xpos, ypos, zpos,
                  size = 10, col = col.points, cex = cex.points,
                  bg = bg.points,
                  axes = FALSE, box = FALSE,
                  xlab = "", ylab = "", zlab = "")
      ##
      ## Add points for labels
      ##
      if (points)
        rgl::points3d(xpos, ypos, zpos,
                      pch = pch.points, cex = cex.points, col = col.points,
                      bg = bg.points)
      ##
      ## Print treatment labels
      ##
      if (!is.null(labels))
        for (i in 1:n)
          rgl::text3d(dat.nodes$xpos.labels[i], dat.nodes$ypos.labels[i],
                      dat.nodes$zpos.labels[i],
                      texts = dat.nodes$labels[i],
                      cex = cex,
                      adj = c(dat.nodes$adj.x[i], dat.nodes$adj.y[i]))
      ##
      ## Add highlighted comparisons
      ##
      if (!is.null(highlight)) {
        for (high in highlight) {
          highs <- unlist(compsplit(high, split = highlight.split))
          if (length(highs) != 2)
            stop("Wrong format for argument 'highlight' (see helpfile of plotgraph command).")
          ##
          if (sum(dat.nodes$trts %in% highs) != 2)
            stop(paste("Argument 'highlight' must contain two of the following values (separated by \":\"):\n  ",
                       paste(paste("'", dat.nodes$trts, "'", sep = ""),
                             collapse = " - "), sep = ""))
          ##
          dat.high <- dat.nodes[dat.nodes$trts %in% highs, ]
          ##
          rgl::lines3d(dat.high$xpos * (1 + 1e-4), dat.high$ypos * (1 + 1e-4),
                       dat.high$zpos * (1 + 1e-4),
                       lwd = W.matrix[trts == highs[1], trts == highs[2]],
                       col = col.highlight)
        }
      }
      ##
      ## Add coloured regions for multi-arm studies
      ##
      if (multiarm) {
        ##
        morethan3 <- FALSE
        ##
        if (n.multi > 0) {
          multiarm.treat <- vector("list", n.multi)
          if (length(col.polygon) == 1)
            col.polygon <- rep(col.polygon, n.multi)
          for (i in 1:n.multi) {
            treat1 <- x$treat1[x$studlab %in% multiarm.studies[i]]
            treat2 <- x$treat2[x$studlab %in% multiarm.studies[i]]
            multiarm.treat[[i]] <- sort(unique(c(treat2, treat1)))
            ##
            dat.multi <- dat.nodes[dat.nodes$trts %in% multiarm.treat[[i]], ]
            if (nrow(dat.multi) == 0)
              dat.multi <- dat.nodes[dat.nodes$trts %in% multiarm.treat[[i]], ]
            if (nrow(dat.multi) == 3)
              rgl::triangles3d(dat.multi$xpos, dat.multi$ypos, dat.multi$zpos,
                               col = col.polygon[i])
            else
              morethan3 <- TRUE
          }
        }
        if (morethan3)
          warning("Multi-arm studies with more than three treatments not shown in 3-D plot.")
      }
      ##
      ## Draw lines
      ##
      for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
          if (A.sign[i, j] > 0) {
            rgl::lines3d(c(xpos[i], xpos[j]), c(ypos[i], ypos[j]), c(zpos[i], zpos[j]),
                         lwd = W.matrix[i, j],
                         col = col)
          }
        }
      }
    }
  }
  
  
  dat.nodes$xpos[is.zero(dat.nodes$xpos)] <- 0
  dat.nodes$ypos[is.zero(dat.nodes$ypos)] <- 0
  ##
  if (!is_2d) {
    dat.nodes$zpos[is.zero(dat.nodes$zpos)] <- 0
    ##
    dat.nodes$xpos.labels <- NULL
    dat.nodes$ypos.labels <- NULL
  }
  else {
    dat.nodes$xpos.labels[is.zero(dat.nodes$xpos.labels)] <- 0
    dat.nodes$ypos.labels[is.zero(dat.nodes$ypos.labels)] <- 0
  }
  ##
  dat.nodes$zpos.labels <- NULL


  dat.edges$xpos[is.zero(dat.edges$xpos)] <- 0
  dat.edges$ypos[is.zero(dat.edges$ypos)] <- 0
  
  
  invisible(list(nodes = dat.nodes, edges = dat.edges))
}
