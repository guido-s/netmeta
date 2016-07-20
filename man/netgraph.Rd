\name{netgraph}

\alias{netgraph}

\title{Network graph}

\description{
  This function generates a graph of the evidence network.
}

\usage{
netgraph(x, seq = x$seq,
         labels = rownames(x$TE.fixed),
         cex = 1, col = "slateblue", offset = 0.0175, scale = 1.10,
         plastic, thickness, lwd = 5, lwd.min = lwd/2.5, lwd.max = lwd*4,
         dim = "2d",
         highlight = NULL, col.highlight = "red2",
         lwd.highlight = lwd, highlight.split = ":",
         multiarm = any(x$narms > 2), col.multiarm = NULL,
         alpha.transparency = 0.5,
         points = FALSE, col.points = "red", cex.points = 1, pch.points = 20,
         number.of.studies = FALSE,
         cex.number.of.studies = cex,
         col.number.of.studies = "white",
         bg.number.of.studies = "black",
         start.layout = ifelse(dim  == "2d", "circle", "eigen"),
         eig1 = 2, eig2 = 3, eig3 = 4,
         iterate, tol = 0.0001, maxit = 500, allfigures = FALSE,
         A.matrix = x$A.matrix, N.matrix = sign(A.matrix),
         xpos = NULL, ypos = NULL, zpos = NULL,
         ...)
}

\arguments{
  \item{x}{An object of class \code{netmeta} (mandatory).}
  \item{seq}{A character or numerical vector specifying the sequence
    of treatments arrangement (anticlockwise if
    \code{start.layout = "circle"}).}
  \item{labels}{An optional vector with treatment labels.}
  \item{cex}{The magnification to be used for treatment labels.}
  \item{col}{Color of lines connecting treatments if argument
    \code{plastic = FALSE}.}
  \item{offset}{Distance between edges (i.e. treatments) in graph and
    treatment labels (value of 0.0175 corresponds to a difference of
    1.75\% of the range on x- and y-axis).}
  \item{scale}{Additional space added outside of edges
    (i.e. treatments). Increase this value for larger treatment labels
    (value of 1.10 corresponds to an additional space of 10\% around
    the network graph).}
  \item{plastic}{A logical indicating whether the appearance of the
    comparisons should be in '3D look' (not to be confused with
    argument \code{dim}).}
  \item{thickness}{Either a character variable to determine the method
    to plot line widths (see Details) or a matrix of the same
    dimension and row and column names as argument \code{A.matrix}
    with information on line width.}
  \item{lwd}{A numeric for scaling the line width of comparisons.}
  \item{lwd.max}{Maximum line width in network graph. The connection
    with the largest value according to argument \code{thickness} will
    be set to this value.}
  \item{lwd.min}{Minimum line width in network graph. All connections
    with line widths below this values will be set to \code{lwd.min}.}
  \item{dim}{A character string indicating whether a 2- or
    3-dimensional plot should be produced, either \code{"2d"} or
    \code{"3d"}.}
  \item{highlight}{A character vector identifying comparisons that
    should be marked in the network graph,
    e.g. \code{highlight = "treat1:treat2"}.}
  \item{col.highlight}{Color for highlighting the comparisons given by
    \code{highlight}.}
  \item{lwd.highlight}{A numeric for the line width for highlighting
    the comparisons given by \code{highlight}.}
  \item{highlight.split}{A character defining splitting criterion
    (only necessary if colon is used in treatment labels).}
  \item{multiarm}{A logical indicating whether multi-arm studies
    should marked in plot.}
  \item{col.multiarm}{Either a function from R library colorspace or
    grDevice to define colors for multi-arm studies or a character
    vector with colors to highlight multi-arm studies.}
  \item{alpha.transparency}{The alpha transparency of colors used to
    highlight multi-arm studies (0 means transparent and 1 means
    opaque).}
  \item{points}{A logical indicating whether points should be printed
    at nodes (i.e. treatments) of the network graph.}
  \item{col.points, cex.points, pch.points}{Corresponding color, size,
    type for points.}
  \item{number.of.studies}{A logical indicating whether number of
    studies should be added to network graph.}
  \item{cex.number.of.studies}{The magnification to be used for number
    of studies.}
  \item{col.number.of.studies}{Color for number of studies.}
  \item{bg.number.of.studies}{Color for shadow around number of
    studies.}
  \item{start.layout}{A character string indicating which starting
    layout is used if \code{iterate = TRUE}. If "circle" (default), the
    iteration starts with a circular ordering of the vertices; if
    "eigen", eigenvectors of the Laplacian matrix are used, calculated
    via generic function \code{\link{eigen}} (spectral decomposition);
    if "prcomp", eigenvectors of the Laplacian matrix are calculated
    via generic function \code{\link{prcomp}} (principal component
    analysis); if "random", a random layout is used, drawn from a
    bivariate normal.}
  \item{eig1}{A numeric indicating which eigenvector is used as x
    coordinate if \code{start = "eigen"} or \code{"prcomp"} and
    \code{iterate = TRUE}. Default is 2, the eigenvector to the
    second-smallest eigenvalue of the Laplacian matrix.}
  \item{eig2}{A numeric indicating which eigenvector is used as
    y-coordinate if \code{start = "eigen"} or \code{"prcomp"} and
    \code{iterate = TRUE}. Default is 3, the eigenvector to the
    third-smallest eigenvalue of the Laplacian matrix.}
  \item{eig3}{A numeric indicating which eigenvector is used as
    z-coordinate if \code{start = "eigen"} or \code{"prcomp"} and
    \code{iterate = TRUE}. Default is 4, the eigenvector to the
    fourth-smallest eigenvalue of the Laplacian matrix.}
  \item{iterate}{A logical indicating whether the stress majorization
    algorithm is carried out for optimization of the layout.}
  \item{tol}{A numeric for the tolerance for convergence if
    \code{iterate = TRUE}.}
  \item{maxit}{An integer defining the maximum number of iteration
    steps if \code{iterate = TRUE}.}
  \item{allfigures}{A logical indicating whether all iteration steps
    are shown if \code{iterate = TRUE}. May slow down calculations if
    set to \code{TRUE} (especially if \code{plastic = TRUE}).}
  \item{A.matrix}{Adjacency matrix (\emph{n}x\emph{n}) characterizing
    the structure of the network graph. Row and column names must be
    the same set of values as provided by argument \code{seq}.}
  \item{N.matrix}{Neighborhood matrix (\emph{n}x\emph{n}) replacing
    A.matrix if neighborhood is to be specified differently from node
    adjacency in the network graph, for example content-based. Row and
    column names must be the same set of values as provided by
    argument \code{seq}.}
  \item{xpos}{Vector (\emph{n}) of x coordinates.}
  \item{ypos}{Vector (\emph{n}) of y coordinates.}
  \item{zpos}{Vector (\emph{n}) of z coordinates.}
  \item{\dots}{Additional graphical arguments.}
}

\details{
  The network is laid out in the plane, where the nodes in the graph
  layout correspond to the treatments and edges display the observed
  treatment comparisons. For the default setting, nodes are placed on
  a circle.  Other starting layouts are "eigen", "prcomp", and
  "random" (Rücker & Schwarzer 2015). If \code{iterate = TRUE}, the
  layout is further optimized using the stress majorization
  algorithm. This algorithm specifies an 'ideal' distance (e.g., the
  graph distance) between two nodes in the plane. In the optimal
  layout, these distances are best approximated in the sense of least
  squares. Starting from an initial layout, the optimum is
  approximated in an iterative process called stress majorization
  (Kamada and Kawai 1989, Michailidis and de Leeuw 2001, Hu 2012). The
  starting layout can be chosen as a circle or coming from
  eigenvectors of the Laplacian matrix (corresponding to Hall's
  algorithm, Hall 1970), calculated in different ways, or
  random. Moreover, it can be chosen whether the iteration steps are
  shown (argument \code{allfigures = TRUE}).

  Argument \code{thickness} providing the line width of the nodes
  (comparisons) can be a matrix of the same dimension as argument
  \code{A.matrix} or any of the following character variables:
   \itemize{
   \item Same line width (argument \code{lwd}) for all comparisons
     (\code{thickness = "equal"})
   \item Proportional to number of studies comparing two treatments
     (\code{thickness = "number.of.studies"})
   \item Proportional to inverse standard error of fixed effect model
     comparing two treatments (\code{thickness = "se.fixed"})
   \item Proportional to inverse standard error of random effects
     model comparing two treatments (\code{thickness = "se.random"})
   \item Weight from fixed effect model comparing two treatments
     (\code{thickness = "w.fixed"})
   \item Weight from random effects model comparing two treatments
     (\code{thickness = "w.random"})
  }
  Only evidence from direct treatment comparisons is considered to
  determine the line width if argument \code{thickness} is equal to
  any but the first method. By default, \code{thickness = "se.fixed"} is
  used if \code{start.layout = "circle"}, \code{iterate = FALSE}, and
  \code{plastic = TRUE}. Otherwise, the same line width is used.
    
  Further, a couple of graphical parameters can be specified, such as
  color and appearance of the edges (treatments) and the nodes
  (comparisons), whether special comparisons should be highlighted and
  whether multi-arm studies should be indicated as colored
  polygons. By default, if R package colorspace is available the
  \code{\link[colorspace]{sequential_hcl}} function is used to
  highlight multi-arm studies; otherwise the \code{\link{rainbow}} is
  used.
}

\references{
  Hall, K.M. (1970).
  An r-dimensional quadratic placement algorithm.
  \emph{Management Science}, \bold{17}, 219--229.

  Hu, Y. (2012).
  \emph{Combinatorial Scientific Computing}, Chapter Algorithms for
  Visualizing Large Networks, pages 525--549.
  Chapman and Hall/CRC Computational Science.

  Kamada, T. and Kawai, S. (1989).
  An algorithm for drawing general undirected graphs.
  \emph{Information Processing Letters}, \bold{31}(1), 7--15.

  Krahn U, Binder H, König J (2013),
  A graphical tool for locating inconsistency in network
  meta-analyses.
  \emph{BMC Medical Research Methodology}, \bold{13}, 35.
  
  Michailidis, G. and de Leeuw, J. (2001).
  Data visualization through graph drawing.
  \emph{Computational Statistics}, \bold{16}(3), 435--450.

  Rücker G & Schwarzer G (2015),
  Automated drawing of network plots in network meta-analysis.
  \emph{Research Synthesis Methods},
  First published online 2015 Jun 9.
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Ulrike Krahn \email{ulrike.krahn@bayer.com}, Jochem König \email{koenigjo@uni-mainz.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{netmeta}}}

\examples{
data(Senn2013)

# Generation of an object of class 'netmeta' with reference treatment 'plac'
#
net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
        data = Senn2013, sm = "MD", reference = "plac")

# Network graph with default settings
#
netgraph(net1)

# Network graph with specified order of the treatments and one
# highlighted comparison
#
trts <- c("plac", "benf", "migl", "acar", "sulf",
          "metf", "rosi", "piog", "sita", "vild")
netgraph(net1, highlight = "acar:metf", seq = trts)

# Same network graph using argument 'seq' in netmeta function
#
net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data = Senn2013, sm = "MD", reference = "plac",
                seq = trts)
netgraph(net2, highlight = "acar:metf")

# Network graph optimized, starting from a circle, with multi-arm
# study colored
#
netgraph(net1, start = "circle", iterate = TRUE, col.multiarm = "purple")

# Network graph optimized, starting from a circle, with multi-arm
# study colored and all intermediate iteration steps visible
#
\dontrun{netgraph(net1, start = "circle", iterate = TRUE, col.multiarm = "purple",
         allfigures = TRUE)}

# Network graph optimized, starting from Laplacian eigenvectors, with
# multi-arm study colored
#
netgraph(net1, start = "eigen", col.multiarm = "purple")

# Network graph optimized, starting from different Laplacian
# eigenvectors, with multi-arm study colored
#
netgraph(net1, start = "prcomp", col.multiarm = "purple")

# Network graph optimized, starting from random initial layout, with
# multi-arm study colored
#
netgraph(net1, start = "random", col.multiarm = "purple")

# Network graph without plastic look and one highlighted comparison
#
netgraph(net1, plastic = FALSE, highlight = "acar:metf")

# Network graph without plastic look and comparisons with same
# thickness
netgraph(net1, plastic = FALSE, thickness = FALSE)

# Network graph with changed labels and specified order of the
# treatments
#
netgraph(net1, seq = c(1, 3, 5, 2, 9, 4, 7, 6, 8, 10),
         labels = LETTERS[1:10])

\dontrun{
# Network graph in 3-D (opens a new device, where you may rotate and
# zoom the plot using the mouse / the mouse wheel).
#
netgraph(net1, dim = "3d")
}
}

\keyword{Network meta-analysis}
\keyword{Network graph}
\keyword{hplot}
