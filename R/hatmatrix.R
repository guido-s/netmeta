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
#' @param comb.fixed A logical indicating whether a hat matrix should
#'   be printed for the fixed effects (common effects) network
#'   meta-analysis.
#' @param comb.random A logical indicating whether a hat matrix should
#'   be printed for the random effects network meta-analysis.
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
#' \emph{k}. Let seTE.adj.fixed and seTE.adj.random be the vectors of
#' adjusted standard errors under the fixed effect and random effects
#' model (see \code{\link{netmeta}}). Let \strong{W} be the
#' \emph{m} x \emph{m} diagonal matrix that contains the inverse
#' variance 1 / seTE.adj.fixed\eqn{^2} or 1 / seTE.adj.random\eqn{^2}.
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
#' One of the following hat matrices is estimated if \code{method
#' = "Davies"}.
#'
#' Use of \code{type = "short"} (default) results in a hat matrix of
#' dimension \emph{e} x \emph{e}, where \emph{e} is the number of
#' (unique) edges (direct comparisons) in the network.
#'
#' Use of \code{type = "long"} results in a hat matrix of dimension
#' \emph{n(n-1)/2} x \emph{e}. This hat matrix describes the flow of
#' evidence through each direct comparison for every possible pair of
#' treatments (regardless of whether there is direct evidence for this
#' pair).
#'
#' Use of \code{type = "full"} results in a hat matrix of dimension
#' \emph{n(n-1)/2} x \emph{n(n-1)/2}. In comparison to the long hat
#' matrix, columns of zeroes are added for comparisons that do not
#' have any direct evidence. This hat matrix is used to calculate the
#' transition matrices for the random walk in
#' \code{\link{netcontrib}}.
#' }
#'
#' @return
#' A list with two hat matrices: \code{fixed} (fixed effect model) and
#' \code{random} (random effects model).
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
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
#'                data = first10, sm = "OR")
#' net1 <- netmeta(p1, comb.fixed = FALSE)
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
                      comb.fixed = x$comb.fixed,
                      comb.random = x$comb.random) {
  
  
  meta:::chkclass(x, "netmeta")
  ##
  setchar <- meta:::setchar
  ##
  method <- setchar(method, c("Ruecker", "Krahn", "Davies"))
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
  

  res <- list()
  ##
  if (method == "Ruecker") {
    res$fixed <- x$H.matrix.fixed
    res$random <- x$H.matrix.random
  }
  else if (method == "Krahn") {
    if (type == "design") {
      res$fixed <- nma.krahn(x)$H
      res$random <- nma.krahn(x, tau.preset = x$tau)$H
    }
    else if (type == "studies") {
      res$fixed <- nma.krahn(x)$H.studies
      res$random <- nma.krahn(x, tau.preset = x$tau)$H.studies
    }
  }
  else if (method == "Davies") {
    res$fixed <- hatmatrix.aggr(x, "fixed", type)
    res$random <- hatmatrix.aggr(x, "random", type)
    ##
    ## Row and col names
    ##
    if (type == "short") {
      rownames(res$fixed) <- colnames(res$fixed) <- x$comparisons
      rownames(res$random) <- colnames(res$random) <- x$comparisons
    }
    else {
      allcomps <- ""
      k <- 0
      for (i in 1:(x$n - 1)) {
        for (j in (i + 1):x$n) {
          k <- k + 1
          allcomps[k] <- paste(x$trts[i], x$trts[j], sep = x$sep.trts)
        }
      }
    }
    if (type == "long") {
      rownames(res$fixed) <- rownames(res$random) <- allcomps
      colnames(res$fixed) <- x$comparisons
      colnames(res$random) <- x$comparisons
    }
    else if (type == "full") {
      rownames(res$fixed) <- colnames(res$fixed) <- allcomps
      rownames(res$random) <- colnames(res$random) <- allcomps
    }
  }
  ##
  res$method <- method
  res$type <- type
  res$comb.fixed <- comb.fixed
  res$comb.random <- comb.random
  res$x <- x
  ##
  class(res) <- "hatmatrix"
  ##
  res
}





#' @rdname hatmatrix
#' @method print hatmatrix
#' @export
#' @export print.hatmatrix


print.hatmatrix <- function(x,
                            comb.fixed = x$comb.fixed,
                            comb.random = x$comb.random,
                            nchar.trts = x$x$nchar.trts,
                            nchar.studlab = 666,
                            digits = gs("digits"),
                            legend = TRUE,
                            legend.studlab = TRUE,
                            ...) {
  
  meta:::chkclass(x, "hatmatrix")
  ##
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  meta:::chknumeric(nchar.trts, length = 1)
  meta:::chknumeric(nchar.studlab, length = 1)
  meta:::chknumeric(digits, length = 1)
  meta:::chklogical(legend)
  meta:::chklogical(legend.studlab)
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
  trts.abbr <- treats(trts, nchar.trts)
  sep.trts <- x$x$sep.trts
  anystudy.r <- anystudy.c <- FALSE
  anycomp.r <- anycomp.c <- FALSE
  ##
  fixed <- x$fixed
  random <- x$random
  ##
  legend <- legend & (comb.fixed | comb.random)
  legend.studlab <- legend.studlab & (comb.fixed | comb.random)
  ##
  if (comb.fixed) {
    compnames <- any(grepl(sep.trts, rownames(fixed), fixed = TRUE))
    anystudy.r <- anystudy.r | !compnames
    anycomp.r  <- anycomp.r  | compnames
    ##
    if (compnames)
      rownames(fixed) <- comps(fixed, trts, sep.trts, nchar.trts)
    else
      rownames(fixed) <- treats(fixed, nchar.studlab)
    ##
    compnames <- any(grepl(sep.trts, colnames(fixed), fixed = TRUE))
    anystudy.c <- anystudy.c | !compnames
    anycomp.c  <- anycomp.c  | compnames
    ##
    if (compnames)
      colnames(fixed) <- comps(fixed, trts, sep.trts, nchar.trts,
                               row = FALSE)
    else
      colnames(fixed) <- treats(fixed, nchar.studlab, row = FALSE)
    ##
    cat("Fixed effects model:\n\n")
    prmatrix(round(fixed, digits))
    if (comb.random)
      cat("\n")
  }
  if (comb.random) {
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
  ## Add legend
  ##
  if (legend & (anycomp.r | anycomp.c)) {
    if (any(trts != trts.abbr)) {
      tmat <- data.frame(trts.abbr, trts)
      names(tmat) <- c("Abbreviation", "Treatment name")
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat("\nLegend:\n")
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(trts.abbr)))
      ##
      if (legend.studlab & (anystudy.r | anystudy.c))
        cat("\n")
    }
    else
      legend <- FALSE
  }
  ##
  if (legend.studlab & (anystudy.r | anystudy.c)) {
    if (anystudy.r) {
      if (comb.fixed) {
        studlab <- rownames(x$fixed)
        studlab.abbr <- rownames(fixed)
      }
      else {
        studlab <- rownames(x$random)
        studlab.abbr <- rownames(random)
      }
    }
    ##
    if (anystudy.c) {
      if (comb.fixed) {
        studlab <- colnames(x$fixed)
        studlab.abbr <- colnames(fixed)
      }
      else {
        studlab <- colnames(x$random)
        studlab.abbr <- colnames(random)
      }
    }
    ##
    if (any(studlab != studlab.abbr)) {
      tmat <- data.frame(studlab.abbr, studlab)
      names(tmat) <- c("Abbreviation", "Study label")
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      if (!(legend & (anycomp.r | anycomp.c)))
        cat("\nLegend:\n")
      prmatrix(unique(tmat), quote = FALSE, right = TRUE,
               rowlab = rep("", nrow(unique(tmat))))
    }
  }
  ##
  invisible(NULL)
}





hatmatrix.aggr <- function(x, model, type) {
  
  model <- meta:::setchar(model, c("fixed", "random"))
  type <- meta:::setchar(type, c("full", "long", "short"))
  
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
          if (model == "fixed")
            WAB <- WAB + 1.0 / (x$seTE.adj.fixed[k])^2
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
    return(H.short)
  else if (type == "long")
    return(H.long)
  else if (type == "full")
    return(H.full)
  else
    return(NULL)
}





hatmatrix.F1000 <- function(x, model) {
  ##
  ## H matrix
  ##
  if (model == "fixed")
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
