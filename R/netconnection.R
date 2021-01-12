#' Get information on network connectivity (number of subnetworks,
#' distance matrix)
#' 
#' @description
#' To determine the network structure and to test whether a given
#' network is fully connected. Network information is provided as a
#' triple of vectors \code{treat1}, \code{treat2}, and \code{studlab}
#' where each row corresponds to an existing pairwise treatment
#' comparison (\code{treat1}, \code{treat2}) in a study
#' (\code{studlab}). The function calculates the number of subnetworks
#' (connectivity components; value of 1 corresponds to a fully
#' connected network) and the distance matrix (in block-diagonal form
#' in the case of subnetworks). If some treatments are combinations of
#' other treatments or have common components, an analysis based on
#' the additive network meta-analysis model might be possible, see
#' \link{discomb} function.
#' 
#' @aliases netconnection print.netconnection
#' 
#' @param treat1 Label / number for first treatment.
#' @param treat2 Label / number for second treatment.
#' @param studlab An optional - but important! - vector with study
#'   labels (see Details).
#' @param data An optional data frame containing the study
#'   information.
#' @param subset An optional vector specifying a subset of studies to
#'   be used.
#' @param title Title of meta-analysis / systematic review.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param warn A logical indicating whether warnings should be
#'   printed.
#' @param x An object of class \code{netconnection}.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param details A logical indicating whether to print the distance
#'   matrix.
#' @param ... Additional arguments (ignored at the moment)
#' 
#' @return
#' An object of class \code{netconnection} with corresponding
#' \code{print} function. The object is a list containing the
#' following components:
#' \item{treat1, treat2, studlab, title, warn, nchar.trts}{As defined
#'   above.}
#' \item{k}{Total number of studies.}
#' \item{m}{Total number of pairwise comparisons.}
#' \item{n}{Total number of treatments.}
#' \item{n.subnets}{Number of subnetworks; equal to 1 for a fully
#'   connected network.}
#' \item{D.matrix}{Distance matrix.}
#' \item{A.matrix}{Adjacency matrix.}
#' \item{L.matrix}{Laplace matrix.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netdistance}},
#'   \code{\link{discomb}}
#' 
#' @examples
#' data(Senn2013)
#' 
#' nc1 <- netconnection(treat1, treat2, studlab, data = Senn2013)
#' nc1
#' 
#' # Extract number of (sub)networks
#' #
#' nc1$n.subnets
#' 
#' # Extract distance matrix
#' #
#' nc1$D.matrix
#' 
#' # Conduct network meta-analysis (results not shown)
#' #
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab, data = Senn2013)
#' 
#' # Artificial example with two subnetworks
#' #
#' t1 <- c("G", "B", "B", "D", "A", "F")
#' t2 <- c("B", "C", "E", "E", "H", "A")
#' #
#' nc2 <- netconnection(t1, t2)
#' print(nc2, details = TRUE)
#' 
#' # Number of subnetworks
#' #
#' nc2$n.subnets
#' 
#' # Extract distance matrix
#' #
#' nc2$D.matrix
#' 
#' # Conduct network meta-analysis (results in an error message due to
#' # unconnected network)
#' try(net2 <- netmeta(1:6, 1:6, t1, t2, 1:6))
#' 
#' # Conduct network meta-analysis on first subnetwork
#' #
#' net2.1 <- netmeta(1:6, 1:6, t1, t2, 1:6,
#'                   subset = (t1 %in% c("A", "F", "H") & t2 %in% c("A", "F", "H")))
#' 
#' # Conduct network meta-analysis on second subnetwork
#' #
#' net2.2 <- netmeta(1:6, 1:6, t1, t2, 1:6,
#'                   subset = !(t1 %in% c("A", "F", "H") & t2 %in% c("A", "F", "H")))
#' 
#' summary(net2.1)
#' summary(net2.2)
#' 
#' @rdname netconnection
#' @export netconnection


netconnection <- function(treat1, treat2, studlab,
                          data = NULL, subset = NULL,
                          nchar.trts = 666,
                          title = "", warn = FALSE) {
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  if (missing(treat1))
    stop("Argument 'treat1' is mandatory.")
  if (missing(treat2))
    stop("Argument 'treat2' is mandatory.")
  ##
  meta:::chknumeric(nchar.trts, min = 1, length = 1)
  ##
  meta:::chklogical(warn)
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  if (is.null(data))
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch treat1, treat2, studlab from data:
  ##
  treat1 <- eval(mf[[match("treat1", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  k.All <- length(treat1)
  ##
  treat2 <- eval(mf[[match("treat2", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  if (length(studlab) != 0)
    studlab <- as.character(studlab)
  else {
    if (warn)
      warning("No information given for argument 'studlab'. ",
              "Assuming that comparisons are from independent studies.")
    studlab <- as.character(seq_along(treat1))
  }
  ##
  ## Catch subset from data:
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  missing.subset <- is.null(subset)
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  fun <- "netconnection"
  ##
  meta:::chklength(treat2, k.All, fun,
                   text = paste0("Arguments 'treat1' and 'treat2' ",
                                 "must have the same length."))
  meta:::chklength(studlab, k.All, fun,
                   text = paste0("Arguments 'treat1' and 'studlab' ",
                                 "must have the same length."))
  ##
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  
  
  ##
  ##
  ## (4) Use subset for analysis
  ##
  ##
  if (!missing.subset) {
    if ((is.logical(subset) & (sum(subset) > k.All)) ||
        (length(subset) > k.All))
      stop("Length of subset is larger than number of studies.")
    ##
    treat1 <- treat1[subset]
    treat2 <- treat2[subset]
    studlab <- studlab[subset]
  }
  
  
  ##
  ##
  ## (5) Additional checks
  ##
  ##
  if (any(treat1 == treat2))
    stop("Treatments must be different (arguments 'treat1' and 'treat2').")
  ##
  ## Check for correct number of comparisons
  ##
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)
      abs(x - round(x)) < tol
  tabnarms <- table(studlab)
  sel.narms <- !is.wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  ##
  if (sum(sel.narms) == 1)
    stop("Study '", names(tabnarms)[sel.narms],
         "' has a wrong number of comparisons.",
         "\n  Please provide data for all treatment comparisons ",
         "(two-arm: 1; three-arm: 3; four-arm: 6, ...).")
  if (sum(sel.narms) > 1)
    stop("The following studies have a wrong number of comparisons: ",
         paste(paste0("'", names(tabnarms)[sel.narms], "'"),
               collapse = ", "),
         "\n  Please provide data for all treatment comparisons ",
         "(two-arm: 1; three-arm: 3; four-arm: 6, ...).")
  
  
  ##
  ##
  ## (6) Determine (sub)network(s)
  ##
  ##
  treats <- as.factor(c(as.character(treat1), as.character(treat2)))
  trts <- levels(treats)
  ##
  n <- length(trts)   # Number of treatments
  m <- length(treat1) # Number of comparisons
  ##
  ## Edge-vertex incidence matrix
  ##
  treat1.pos <- treats[1:m]
  treat2.pos <- treats[(m + 1):(2 * m)]
  B <- createB(treat1.pos, treat2.pos, ncol = n)
  ##
  rownames(B) <- studlab
  colnames(B) <- trts
  ##
  L.mult <- t(B) %*% B             # Laplacian matrix with multiplicity
  A <- diag(diag(L.mult)) - L.mult # Adjacency matrix
  D <- netdistance(A)              # Distance matrix
  L <- diag(rowSums(A)) - A        # Laplacian matrix without multiplicity
  ##
  n.subsets <- as.integer(table(round(eigen(L)$values, 10) == 0)[2])
  ##
  if (n.subsets > 1) {
    ##
    ## Block diagonal matrix in case of sub-networks
    ##
    maxdist <- dim(D)[1]
    ##
    D2 <- D
    D2[is.infinite(D2)] <- maxdist
    order.D <- hclust(dist(D2))$order
    ##
    D <- D[order.D, order.D]
    A <- A[order.D, order.D]
    L <- L[order.D, order.D]
  }
  
  
  designs <- designs(treat1, treat2, studlab)
  
  
  res <- list(treat1 = treat1,
              treat2 = treat2,
              studlab = studlab,
              ##
              k = length(unique(studlab)),
              m = m,
              n = n,
              n.subnets = n.subsets,
              d = length(unique(designs$design)),
              ##
              D.matrix = D,
              A.matrix = A,
              L.matrix = L,
              ##
              designs = unique(sort(designs$design)),
              ##
              nchar.trts = nchar.trts,
              ##
              title = title,
              ##
              warn = warn,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )

  class(res) <- "netconnection"

  res
}





#' @rdname netconnection
#' @method print netconnection
#' @export
#' @export print.netconnection


print.netconnection <- function(x,
                                digits = max(4, .Options$digits - 3),
                                nchar.trts = x$nchar.trts,
                                details = FALSE,
                                ...) {
  
  meta:::chkclass(x, "netconnection")
  ##
  if (is.null(nchar.trts))
    nchar.trts <- 666
  
  
  meta:::chknumeric(digits, length = 1)
  meta:::chknumeric(nchar.trts, min = 1, length = 1)
  meta:::chklogical(details)
  
  
  matitle(x)
  ##
  cat(paste("Number of studies: k = ", x$k, "\n", sep = ""))
  cat(paste("Number of treatments: n = ", x$n, "\n", sep = ""))
  cat(paste("Number of pairwise comparisons: m = ", x$m, "\n", sep = ""))
  if (!is.null(x$d))
    cat(paste("Number of designs: d = ", x$d, "\n", sep = ""))
  ##
  cat("Number of subnetworks: ", x$n.subnets, "\n", sep = "")

  if (details) {
    cat("\nDistance matrix:\n")
    
    D <- round(x$D.matrix, digits = digits)
    D[is.infinite(D)] <- "."
    ##
    if (x$n.subnets == 1)
      diag(D) <- "."
    ##
    rownames(D) <- treats(rownames(D), nchar.trts)
    colnames(D) <- treats(colnames(D), nchar.trts)
    ##
    prmatrix(D, quote = FALSE, right = TRUE)
    
    
    if (any(rownames(x$D.matrix) != rownames(D))) {
      abbr <- rownames(D)
      full <- rownames(x$D.matrix)
      ##
      tmat <- data.frame(abbr, full)
      names(tmat) <- c("Abbreviation", "Treatment name")
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat("\nLegend:\n")
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(abbr)))
    }
  }
  
  
  invisible(NULL)
}
