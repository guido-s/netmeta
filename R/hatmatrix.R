#' Derive hat matrix from network meta-analysis
#' 
#' @description
#' Auxiliary function to derive hat matrix from network meta-analysis
#' 
#' @param x A \code{\link{netmeta}} object.
#' @param method A character string indicating which method is used to
#'   derive the hat matrix. Either \code{"Ruecker"}, \code{"Krahn"} or
#'   \code{"Davies"} (can be abbreviated, see Details).
#' @param type A character string indicating which specific hat matrix
#'   should be derived (can be abbreviated, see Details).
#' @param common A logical indicating whether a hat matrix should be
#'   printed for the common effects network meta-analysis.
#' @param random A logical indicating whether a hat matrix should be
#'   printed for the random effects network meta-analysis.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param nchar.studlab A numeric defining the minimum number of
#'   characters used to create unique study labels.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param legend.studlab A logical indicating whether a legend should
#'   be printed for abbreviated study labels.
#' @param \dots Additional arguments (ignored).
#' 
#' @details
#' This auxiliary function can be used to derive various hat matrices
#' from a network meta-analysis object.
#' 
#' \subsection{Hat matrix by Rücker (2012)}{
#' 
#' This hat matrix is estimated if \code{method = "Ruecker"}.
#' 
#' Let \emph{n} be the number of different treatments (nodes,
#' vertices) in a network and let \emph{m} be the number of existing
#' comparisons (edges) between the treatments. If there are only
#' two-arm studies, \emph{m} is equal to the number of studies,
#' \emph{k}. Let seTE.adj.common and seTE.adj.random be the vectors of
#' adjusted standard errors under the common and random effects model
#' (see \code{\link{netmeta}}). Let \strong{W} be the \emph{m} x
#' \emph{m} diagonal matrix that contains the inverse variance 1 /
#' seTE.adj.common\eqn{^2} or 1 / seTE.adj.random\eqn{^2}.
#'
#' The given comparisons define the network structure. Therefrom an
#' \emph{m} x \emph{n} design matrix X (edge-vertex incidence matrix) is
#' formed; for more precise information, see Rücker (2012). Moreover,
#' the \emph{n} x \emph{n} Laplacian matrix \strong{L} and its
#' Moore-Penrose pseudoinverse \strong{L}\eqn{^+} are calculated (both
#' matrices play an important role in graph theory and electrical
#' network theory). Using these matrices, the variances based on both
#' direct and indirect comparisons can be estimated. The hat matrix
#' \strong{H} is estimated by \strong{H =
#' XL}\eqn{^+}\strong{X}\eqn{^T}\strong{W = X(X}\eqn{^T}\strong{W
#' X)}\eqn{^+}\strong{X}\eqn{^T}\strong{W}.
#' }
#' 
#' \subsection{Hat matrices by Krahn et al. (2013)}{
#' 
#' One of the following hat matrices is estimated if \code{method
#' = "Krahn"}.
#'
#' Use of \code{type = "design"} (default) results in a hat matrix of
#' dimension \emph{n(n-1)/2} x \emph{d}, where \emph{d} is the sum of
#' the number of independent comparisons from each design.
#'
#' Use of \code{type = "studies"} results in a hat matrix of dimension
#' \emph{n(n-1)/2} x \emph{l}, where \emph{l} is the number of
#' independent pairwise comparisons, i.e., a three-arm study
#' contributes two pairwise comparisons.
#' }
#' 
#' \subsection{Hat matrices by Davies et al. (2021)}{
#' 
#' One of three hat matrices is estimated if \code{method = "Davies"}.
#'
#' Here, we focus on the hat matrix of the aggregate (two-step)
#' version of the graph theoretical NMA model. In the first step, a
#' pairwise meta-analysis is performed across each edge using the
#' adjusted weights (these account for correlations due to multi-arm
#' trials). From this we obtain direct treatment effect estimates (and
#' corresponding aggregate weights) associated with each edge. In step
#' two, we combine these direct estimates in a network meta analysis
#' to obtain the network estimates.  This is done using weighted least
#' squares regression. The hat matrix associated with this second step
#' is called the \emph{aggregate hat matrix}.
#' 
#' All three versions of the aggregate hat matrix contain the same
#' information: the second two can be derived directly from the
#' first. They differ in their dimensionality.
#' 
#' Each row of the hat matrix that represents a treatment comparison
#' (\emph{ij}) describes the flow of evidence through each edge for
#' that comparison. This defines a directed acyclic 'flow graph' from
#' node \emph{i} to node \emph{j}.
#'
#' (1) Use of \code{type = "short"} (default) results in a hat matrix
#' of dimension \emph{e} x \emph{e}, where \emph{e} is the number of
#' (unique) edges (direct comparisons) in the network. This is the
#' aggregate hat matrix described in Davies et al. (2021). Each row
#' and column represents a pair of treatments for which there is at
#' least one direct comparison.
#'
#' (2) Use of \code{type = "long"} results in a hat matrix of
#' dimension \emph{n(n-1)/2} x \emph{e}. There is a row for every
#' possible pair of treatments in the network - regardless of whether
#' there is direct evidence for this comparison. Each column
#' represents a pair of treatments for which there is at least one
#' direct comparison. The extra rows can be calculated from the short
#' hat matrix using consistency equations.
#' 
#'
#' (3) Use of \code{type = "full"} results in a hat matrix of
#' dimension \emph{n(n-1)/2} x \emph{n(n-1)/2}. In comparison to the
#' long hat matrix, columns of zeroes are added for comparisons that
#' do not have any direct evidence. Therefore, there is a row and
#' column for every pair of treatments in the network. This hat matrix
#' is used to calculate the transition matrices for the random walk in
#' \code{\link{netcontrib}}.
#' }
#'
#' @return
#' A list with two hat matrices: \code{common} (common effects model)
#' and \code{random} (random effects model).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netcontrib}},
#'   \code{\link{netheat}}
#' 
#' @references
#' Davies AL, Papakonstantinou T, Nikolakopoulou A, Rücker G, Galla T
#' (2021):
#' Network meta-analysis and random walks.
#' Available from: http://arxiv.org/abs/2107.02886
#' 
#' Krahn U, Binder H, König J (2013):
#' A graphical tool for locating inconsistency in network meta-analyses.
#' \emph{BMC Medical Research Methodology},
#' \bold{13}, 35
#'
#' Rücker G (2012):
#' Network meta-analysis, electrical networks and graph theory.
#' \emph{Research Synthesis Methods},
#' \bold{3}, 312--24
#' 
#' @examples
#' data(Dong2013)
#' # Only consider first ten studies for concise output
#' first10 <- subset(Dong2013, id <= 10)
#' p1 <- pairwise(treatment, death, randomized, studlab = id,
#'   data = first10, sm = "OR")
#' net1 <- netmeta(p1, common = FALSE)
#' 
#' hatmatrix(net1)
#' hatmatrix(net1, method = "k")
#' hatmatrix(net1, method = "k", type = "studies")
#' hatmatrix(net1, method = "d")
#' hatmatrix(net1, method = "d", type = "long")
#' hatmatrix(net1, method = "d", type = "full")
#' 
#' @export hatmatrix


hatmatrix <- function(x, method = "Ruecker", type,
                      common = x$common,
                      random = x$random,
                      nchar.trts = x$nchar.trts,
                      nchar.studlab = x$nchar.studlab) {
  
  chkclass(x, "netmeta")
  x <- updateversion(x)
  ##
  method <- setchar(method, c("Ruecker", "Krahn", "Davies"))
  ##
  if (x$n == 2 & method == "Krahn") {
    warning("Hat matrix by Krahn et al. (2013) not available for ",
            "network meta-analysis with only two treatments.",
            call. = FALSE)
    return(NULL)
  }
  ##
  if (missing(type)) {
    if (method == "Ruecker")
      type <- ""
    else if (method == "Krahn")
      type <- "design"
    else if (method == "Davies")
      type <- "short"
  }
  else {
    if (method == "Ruecker") {
      warning("Argument 'type' ignored for argument method = \"Ruecker\".",
              call. = FALSE)
      type <- ""
    }
    else if (method == "Krahn")
      type <- setchar(type, c("design", "studies"))
    else if (method == "Davies")
      type <- setchar(type, c("short", "long", "full"))
  }
  ##
  common <- replaceNULL(x$common, x$comb.fixed)
  common <- replaceNULL(x$common, x$fixed)
  chklogical(common)
  ##
  random <- replaceNULL(x$random, x$comb.random)
  chklogical(random)
  ##
  nchar.trts <- replaceNULL(nchar.trts, 666)
  chknumeric(nchar.trts, min = 1, length = 1)
  ##
  nchar.studlab <- replaceNULL(nchar.studlab, 666)
  chknumeric(nchar.studlab, min = 1, length = 1)
  
  
  res <- list()
  ##
  if (method == "Ruecker") {
    res$common <- x$H.matrix.common
    res$random <- x$H.matrix.random
  }
  else if (method == "Krahn") {
    if (type == "design") {
      res$common <- nma.krahn(x)$H
      res$random <- nma.krahn(x, tau.preset = x$tau)$H
    }
    else if (type == "studies") {
      res$common <- nma.krahn(x)$H.studies
      res$random <- nma.krahn(x, tau.preset = x$tau)$H.studies
    }
  }
  else if (method == "Davies") {
    res$common <- hatmatrix.aggr(x, "common", type)
    res$random <- hatmatrix.aggr(x, "random", type)
  }
  ##
  res$method <- method
  res$type <- type
  ##
  res$x <- x
  res$x$common <- common
  res$x$random <- random
  res$x$nchar.trts <- nchar.trts
  res$x$nchar.studlab <- nchar.studlab
  ##
  res$version <- packageDescription("netmeta")$Version
  ##
  class(res) <- "hatmatrix"
  ##
  res
}





#' @rdname hatmatrix
#' @method print hatmatrix
#' @export


print.hatmatrix <- function(x,
                            common = x$x$common,
                            random = x$x$random,
                            nchar.trts = x$x$nchar.trts,
                            nchar.studlab = x$x$nchar.studlab,
                            digits = gs("digits"),
                            legend = TRUE,
                            legend.studlab = TRUE,
                            ...) {
  
  chkclass(x, "hatmatrix")
  x <- updateversion(x)
  ##
  chklogical(common)
  common.logical <- common
  chklogical(random)
  random.logical <- random
  ##
  chknumeric(nchar.trts, length = 1)
  nchar.studlab <- replaceNULL(nchar.studlab, 666)
  chknumeric(nchar.studlab, length = 1)
  chknumeric(digits, length = 1)
  chklogical(legend)
  chklogical(legend.studlab)
  ##
  matitle(x$x)
  ##
  cat(paste0("Hat matrix (",
             if (x$method == "Ruecker")
               "R\u00FCcker, 2012, Research Synthesis Method"
             else if (x$method == "Krahn")
               "Krahn et al., 2013"
             else if (x$method == "Davies")
               "Davies et al., 2021",
             ")\n\n"))
  ##
  trts <- x$x$trts
  sep.trts <- x$x$sep.trts
  anystudy.r <- anystudy.c <- FALSE
  anycomp.r <- anycomp.c <- FALSE
  ##
  common <- x$common
  random <- x$random
  ##
  legend <- legend & (common.logical | random.logical)
  legend.studlab <- legend.studlab & (common.logical | random.logical)
  ##
  if (common.logical) {
    compnames <- any(grepl(sep.trts, rownames(common), fixed = TRUE))
    anystudy.r <- anystudy.r | !compnames
    anycomp.r  <- anycomp.r  | compnames
    ##
    if (compnames)
      rownames(common) <- comps(common, trts, sep.trts, nchar.trts)
    else
      rownames(common) <- treats(common, nchar.studlab)
    ##
    compnames <- any(grepl(sep.trts, colnames(common), fixed = TRUE))
    anystudy.c <- anystudy.c | !compnames
    anycomp.c  <- anycomp.c  | compnames
    ##
    if (compnames)
      colnames(common) <- comps(common, trts, sep.trts, nchar.trts,
                                row = FALSE)
    else
      colnames(common) <- treats(common, nchar.studlab, row = FALSE)
    ##
    cat("Common effects model:\n\n")
    prmatrix(round(common, digits))
    if (random.logical)
      cat("\n")
  }
  if (random.logical) {
    compnames <- any(grepl(sep.trts, rownames(random), fixed = TRUE))
    anystudy.r <- anystudy.r | !compnames
    anycomp.r  <- anycomp.r  | compnames
    ##
    if (compnames)
      rownames(random) <- comps(random, trts, sep.trts, nchar.trts)
    else
      rownames(random) <- treats(random, nchar.studlab)
    ##
    compnames <- any(grepl(sep.trts, colnames(random), fixed = TRUE))
    anystudy.c <- anystudy.c | !compnames
    anycomp.c  <- anycomp.c  | compnames
    ##
    if (compnames)
      colnames(random) <- comps(random, trts, sep.trts, nchar.trts,
                                row = FALSE)
    else
      colnames(random) <- treats(random, nchar.studlab, row = FALSE)
    ##
    cat("Random effects model:\n\n")
    prmatrix(round(random, digits))
  }
  ##
  ## Add legend with abbreviated treatment labels
  ##
  legend <- legendabbr(trts, treats(trts, nchar.trts),
                       legend && (anycomp.r | anycomp.c))
  ##
  ## Add legend with abbreviated study labels
  ##
  if (legend.studlab && (anystudy.r | anystudy.c)) {
    if (anystudy.r) {
      if (common.logical) {
        studlab <- rownames(x$common)
        studlab.abbr <- rownames(common)
      }
      else {
        studlab <- rownames(x$random)
        studlab.abbr <- rownames(random)
      }
    }
    ##
    if (anystudy.c) {
      if (common.logical) {
        studlab <- colnames(x$common)
        studlab.abbr <- colnames(common)
      }
      else {
        studlab <- colnames(x$random)
        studlab.abbr <- colnames(random)
      }
    }
    ##
    if (legend)
      legendabbr(studlab, studlab.abbr, TRUE, "Study label", "\n")
    else
      legendabbr(studlab, studlab.abbr, TRUE, "Study label")
  }
  ##
  invisible(NULL)
}





hatmatrix.aggr <- function(x, model, type) {
  
  model <- setchar(model, c("common", "random"))
  type <- setchar(type, c("full", "long", "short"))
  
  ## Create aggregate B matrix
  if (!is.null(x$B.matrix.aggr))
    B <- x$B.matrix.aggr
  else
    B <- createB(x$treat1.pos, x$treat2.pos, x$n, aggr = TRUE)
  ## Number of edges (direct comparisons) = no. of rows of B
  e <- nrow(B)
  ##
  ## Aggregate weights
  ##
  ## Create a matrix of the aggregate weights e x e
  ##
  W <- matrix(0, nrow = e, ncol = e)
  for (i in 1:e){
    ## idxA element of row i of B that corresponds to baseline treatment
    ## idxB element of row i of B that corresponds to the other treatment
    for (j in seq_len(x$n))
      if (B[i, j] == 1)
        idxA <- j
      else if (B[i, j] == -1)
        idxB <- j
    ## Find all rows of x$B.matrix (dimension m x n) that have the 
    ## same elements as row i of B and calculate the inverse variance 
    ## of the weighted mean
    WAB <- 0.0
    for (k in seq_len(x$m)) {
      if (x$B.matrix[k, idxA] == 1 & x$B.matrix[k, idxB] == -1) {
        if (model == "common")
          WAB <- WAB + 1.0 / (x$seTE.adj.common[k])^2
        else if (model == "random")
          WAB <- WAB + 1.0 / (x$seTE.adj.random[k])^2
      }
    }
    ##
    W[i, i] = WAB
  }
  ##
  ## Laplacian of the aggregate model
  ##
  L <- t(B) %*% W %*% B
  ## Pseudo-Inverse of L 
  L.plus <- invmat(L)
  L.plus[is.zero(L.plus)] <- 0
  ##
  ## Aggregate Hat matrix
  ##
  ## type = "short"
  ## e x e hat matrix
  H.short <- B %*% L.plus %*% t(B) %*% W
  ##
  ## type = "long"
  ## n (n - 1) / 2 x e hat matrix
  ## Get edge incidence matrix for a fully connected network with n
  ## nodes
  B.full <- createB(ncol = x$n)
  H.long <- B.full %*% L.plus %*% t(B) %*% W
  ## type = "full"
  ## n (n - 1) / 2 x n (n - 1) / 2 hat matrix
  ## Insert columns of zeroes for missing edges
  c <- (x$n) * (x$n - 1.0) / 2.0
  H.full <- matrix(0, nrow = c, ncol = c)
  ##
  k <- 0 # counts H.full rows and columns
  l <- 1 # counts B rows and H.long columns
  for (i in 1:(x$n - 1)) {
    for (j in (i + 1):x$n){
      k <- k + 1
      ## is there a direct comparison of ij?
      if (l <= e){
        ## once l is > e (no. of edges) we have reached the end of B
        if (B[l, i] == 1 & B[l, j] == -1){
          for (m in 1:c)
            H.full[m, k] <- H.long[m, l]
          l <- l+1
        }
      }
    }
  }
  ##
  if (type == "short")
    H <- H.short
  else if (type == "long")
    H <- H.long
  else if (type == "full")
    H <- H.full
  else
    return(NULL)
  ##
  ## Row and col names
  ##
  if (type == "short")
    rownames(H) <- colnames(H) <- x$comparisons
  else {
    allcomps <- ""
    k <- 0
    for (i in 1:(x$n - 1)) {
      for (j in (i + 1):x$n) {
        k <- k + 1
        allcomps[k] <- paste(x$trts[i], x$trts[j], sep = x$sep.trts)
      }
    }
    ##
    if (type == "long") {
      rownames(H) <- allcomps
      colnames(H) <- x$comparisons
    }
    else if (type == "full")
      rownames(H) <- colnames(H) <- allcomps
  }
  
  H
}





hatmatrix.F1000 <- function(x, model) {
  ##
  ## H matrix
  ##
  if (model == "common")
    krahn <- nma.krahn(x)
  else if (model == "random")
    krahn <- nma.krahn(x, tau.preset = x$tau)
  ##
  X.full <- krahn$X.full
  direct <- krahn$direct
  X <- krahn$X.full[direct$comparison, , drop = FALSE]
  Vd <- diag(direct$seTE^2,
             nrow = length(direct$seTE),
             ncol = length(direct$seTE))
  H <-
    X.full %*% solve(t(X) %*% solve(Vd) %*% X) %*% t(X) %*% solve(Vd)
  ##
  colnames(H) <- rownames(X)
  ##
  H
}
