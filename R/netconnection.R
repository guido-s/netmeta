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
#' @aliases netconnection netconnection.default netconnection.pairwise
#'   print.netconnection
#' 
#' @param data A data frame, e.g., created with
#'   \code{\link{pairwise}}.
#' @param treat1 Label / number for first treatment (ignored if
#'   \code{data} was created with \code{\link{pairwise}}).
#' @param treat2 Label / number for second treatment (ignored if
#'   \code{data} was created with \code{\link{pairwise}}).
#' @param studlab Study labels (ignored if \code{data} was created
#'   with \code{\link{pairwise}}).
#' @param subset An optional vector specifying a subset of studies to
#'   be used.
#' @param title Title of meta-analysis / systematic review.
#' @param sep.trts A character used in comparison names as separator
#'   between treatment labels.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param warn A logical indicating whether warnings should be
#'   printed.
#' @param x An object of class \code{netconnection}.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param details A logical indicating whether to print the distance
#'   matrix.
#' @param details.disconnected A logical indicating whether to print
#'   more details for disconnected networks.
#' @param \dots Additional arguments (ignored at the moment)
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
#' @author Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}, Guido
#'   Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
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
#' \dontrun{
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
#' net2.1 <- netmeta(1:6, 1:6, t1, t2, 1:6, subset = nc2$subnet == 1)
#' 
#' # Conduct network meta-analysis on second subnetwork
#' #
#' net2.2 <- netmeta(1:6, 1:6, t1, t2, 1:6, subset = nc2$subnet == 2)
#' 
#' net2.1
#' net2.2
#' }
#' 
#' @rdname netconnection
#' @method netconnection default
#' @export


netconnection.default <- function(data = NULL, treat1, treat2, studlab = NULL,
                                  subset = NULL,
                                  sep.trts = ":",
                                  nchar.trts = 666,
                                  title = "", details.disconnected = FALSE,
                                  warn = FALSE, ...) {
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  ## Catch treat1
  ##
  if (!nulldata & !is.data.frame(data)) {
    if (!missing(treat1) & !missing(treat2) & !missing(studlab))
      stop("Argument 'data' must be a data frame.")
    else if (!missing(treat1) & !missing(treat2)) {
      treat1 <- catch("data", mc, sfsp, sfsp)
      treat2 <- catch("treat1", mc, sfsp, sfsp)
      studlab <- catch("treat2", mc, sfsp, sfsp)
    }
    else if (!missing(treat1)) {
      treat1 <- catch("data", mc, sfsp, sfsp)
      treat2 <- catch("treat1", mc, sfsp, sfsp)
    }
  }
  else {
    ##
    if (nulldata)
      data <- sfsp
    ##
    treat1 <- catch("treat1", mc, data, sfsp)
    treat2 <- catch("treat2", mc, data, sfsp)
    studlab <- catch("studlab", mc, data, sfsp)
    subset <- catch("subset", mc, data, sfsp)
  }
  ##
  if (length(studlab) != 0)
    studlab <- as.character(studlab)
  else {
    if (warn)
      warning("No information given for argument 'studlab'. ",
              "Assuming that comparisons are from independent studies.")
    studlab <- as.character(seq_along(treat1))
  }
  ##
  chknumeric(nchar.trts, min = 1, length = 1)
  ##
  chklogical(details.disconnected)
  chklogical(warn)
  
  
  ##
  ##
  ## (2) Check length of essential variables
  ##
  ##
  
  fun <- "netconnection"
  ##
  k.All <- length(treat1)
  ##
  missing.subset <- is.null(subset)
  ##
  chklength(treat2, k.All, fun,
                   text = paste0("Arguments 'treat1' and 'treat2' ",
                                 "must have the same length."))
  chklength(studlab, k.All, fun,
                   text = paste0("Arguments 'treat1' and 'studlab' ",
                                 "must have the same length."))
  ##
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  
  
  ##
  ##
  ## (3) Use subset for analysis
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
  ## (4) Additional checks
  ##
  ##
  
  if (any(treat1 == treat2))
    stop("Treatments must be different (arguments 'treat1' and 'treat2').")
  ##
  ## Check for correct number of comparisons
  ##
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
  labels <- sort(unique(c(treat1, treat2)))
  ##
  if (compmatch(labels, sep.trts)) {
    if (!missing(sep.trts))
      warning("Separator '", sep.trts,
              "' used in at least one treatment label. ",
              "Try to use predefined separators: ",
              "':', '-', '_', '/', '+', '.', '|', '*'.",
              call. = FALSE)
    ##
    if (!compmatch(labels, ":"))
      sep.trts <- ":"
    else if (!compmatch(labels, "-"))
      sep.trts <- "-"
    else if (!compmatch(labels, "_"))
      sep.trts <- "_"
    else if (!compmatch(labels, "/"))
      sep.trts <- "/"
    else if (!compmatch(labels, "+"))
      sep.trts <- "+"
    else if (!compmatch(labels, "."))
      sep.trts <- "-"
    else if (!compmatch(labels, "|"))
      sep.trts <- "|"
    else if (!compmatch(labels, "*"))
      sep.trts <- "*"
    else
      stop("All predefined separators (':', '-', '_', '/', '+', ",
           "'.', '|', '*') are used in at least one treatment label.",
           "\n   Please specify a different character that should be ",
           "used as separator (argument 'sep.trts').",
           call. = FALSE)
  }
  
  
  ##
  ##
  ## (5) Determine (sub)network(s)
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
  ## Block diagonal matrix in case of sub-networks
  ##
  maxdist <- nrow(D)
  D2 <- D
  D2[is.infinite(D2)] <- maxdist
  o <- hclust(dist(D2))$order
  ##
  D <- D[o, o]
  A <- A[o, o]
  L <- L[o, o]
  ##
  A.loop <- A
  D.loop <- D
  id.treats <- character(0)
  id.subnets <- numeric(0)
  more.subnets <- TRUE
  subnet.i <- 0
  dat.subnet <- data.frame()
  ##
  while (more.subnets) {
    subnet.i <- subnet.i + 1
    n.i <- seq_len(nrow(D.loop))
    #
    next.subnet <- min(max(n.i) + 1, n.i[is.infinite(D.loop[, 1])])
    sel.i <- seq_len(next.subnet - 1)
    #
    id.treats <- c(id.treats, rownames(D.loop)[sel.i])
    id.subnets <- c(id.subnets, rep(subnet.i, length(sel.i)))
    ##
    A.i <- A.loop[sel.i, sel.i]
    A.i <- A.i[order(rownames(A.i)), order(colnames(A.i))]
    ##
    for (row.i in 1:(ncol(A.i) - 1)) {
      for (col.i in 2:ncol(A.i)) {
        if (col.i > row.i) {
          if (A.i[row.i, col.i] > 0) {
            trt1.i <- rownames(A.i)[row.i]
            trt2.i <- rownames(A.i)[col.i]
            #
            dat.subnet <-
              rbind(
                dat.subnet,
                data.frame(subnet = subnet.i,
                           comparison = paste(trt1.i, trt2.i, sep = sep.trts),
                           treat1 = trt1.i, treat2 = trt2.i))
          }
        }
      }
    }
    #
    A.loop <- A.loop[-sel.i, -sel.i]
    D.loop <- D.loop[-sel.i, -sel.i]
    #
    more.subnets <- nrow(D.loop) > 0
  }
  #
  # Order matrices within subnets by treatments
  # (data frame 'dat.subnet' is already sorted by subnetwork and treatments)
  #
  seq <- NULL
  for (i in unique(dat.subnet$subnet)) {
    dat.i <- dat.subnet[dat.subnet$subnet == i, ]
    seq <- c(seq, unique(c(dat.i$treat1, dat.i$treat2)))
  }
  #
  A <- A[seq, seq]
  D <- D[seq, seq]
  L <- L[seq, seq]
  #
  # Add subnetwork number to comparisons
  #
  subnet <- rep(NA, length(treat1))
  #
  for (i in seq_along(id.treats))
    subnet[treat1 == id.treats[i]] <- id.subnets[i]
  ##
  comparisons <- dat.subnet$comparison
  subnet.comparisons <- dat.subnet$subnet
  #
  designs <- designs(treat1, treat2, studlab)
  
  res <- list(treat1 = treat1,
              treat2 = treat2,
              studlab = studlab,
              design = designs$design,
              subnet = subnet,
              #
              k = length(unique(studlab)),
              m = m,
              n = n,
              n.subnets = n.subsets,
              d = length(unique(designs$design)),
              #
              seq = seq,
              #
              D.matrix = D,
              A.matrix = A,
              L.matrix = L,
              ##
              designs = unique(sort(designs$design)),
              comparisons = comparisons,
              subnet.comparisons = subnet.comparisons,
              ##
              nchar.trts = nchar.trts,
              ##
              title = title,
              ##
              details.disconnected = details.disconnected,
              ##
              warn = warn,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )

  class(res) <- "netconnection"

  res
}





#' @rdname netconnection
#' @method netconnection pairwise
#' @export


netconnection.pairwise <- function(data,
                                   treat1, treat2, studlab = NULL,
                                   subset = NULL,
                                   sep.trts = ":",
                                   nchar.trts = 666,
                                   title = "", details.disconnected = FALSE,
                                   warn = FALSE,
                                   ...) {
  
  #
  #
  # (1) Check arguments
  #
  #
  
  chkclass(data, "pairwise")
  #
  chklogical(warn)
  #
  # Arguments 'treat1', 'treat2' and 'studlab' ignored
  #
  if (warn) {
    if (!missing(treat1))
      warning("Argument 'treat1' ignored as argument 'data' is an ",
              "object created with pairwise().",
              call. = FALSE)
    #
    if (!missing(treat2))
      warning("Argument 'treat2' ignored as argument 'data' is an ",
              "object created with pairwise().",
              call. = FALSE)
    #
    if (!missing(studlab))
      warning("Argument 'studlab' ignored as argument 'data' is an ",
              "object created with pairwise().",
              call. = FALSE)
  }
  #
  treat1 <- data$treat1
  treat2 <- data$treat2
  studlab <- data$studlab
  #
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  #
  missing.subset <- missing(subset)
  if (!missing.subset) {
    sfsp <- sys.frame(sys.parent())
    mc <- match.call()
    subset <- catch("subset", mc, data, sfsp)
    #
    k.All <- length(treat1)
    #
    if ((is.logical(subset) & (sum(subset) > k.All)) ||
        (length(subset) > k.All))
      stop("Length of subset is larger than number of studies.")
    #
    treat1 <- treat1[subset]
    treat2 <- treat2[subset]
    studlab <- studlab[subset]
  }
  #
  chknumeric(nchar.trts, min = 1, length = 1)
  #
  chklogical(details.disconnected)
  
  
  #
  #
  # (2) Call netconnection.default()
  #
  #
  
  res <- netconnection(treat1 = treat1, treat2 = treat2,
                       studlab = studlab,
                       sep.trts = sep.trts, nchar.trts = nchar.trts,
                       title = title,
                       details.disconnected = details.disconnected,
                       warn = warn,
                       ...)
  #
  res
}





#' @rdname netconnection
#' @method netconnection netmeta
#' @export

netconnection.netmeta <- function(data,
                                  sep.trts = data$sep.trts,
                                  nchar.trts = data$nchar.trts,
                                  title = data$title,
                                  details.disconnected = FALSE,
                                  warn = FALSE, ...) {
  
  chkclass(data, "netmeta")
    
  res <- netconnection(treat1 = data$treat1, treat2 = data$treat2,
                       studlab = data$studlab,
                       sep.trts = sep.trts, nchar.trts = nchar.trts,
                       title = title)
  #
  res
}





#' @rdname netconnection
#' @method netconnection netcomb
#' @export

netconnection.netcomb <- function(data,
                                  sep.trts = data$sep.trts,
                                  nchar.trts,
                                  title = data$title,
                                  details.disconnected = FALSE,
                                  warn = FALSE, ...) {
  
  chkclass(data, "netcomb")
  
  if (inherits(data, "discomb")) {
    if (missing(nchar.trts))
      nchar.trts <- data$nchar.comps
    #
    return(netconnection(treat1 = data$treat1, treat2 = data$treat2,
                         studlab = data$studlab,
                         sep.trts = sep.trts, nchar.trts = nchar.trts,
                         title = title))
  }
  else {
    if (missing(nchar.trts))
      nchar.trts <- data$nchar.trts
    #
    return(netconnection(treat1 = data$treat1, treat2 = data$treat2,
                         studlab = data$studlab,
                         sep.trts = sep.trts, nchar.trts = nchar.trts,
                         title = title))
  }
}





#' @rdname netconnection
#' @method print netconnection
#' @export


print.netconnection <- function(x,
                                digits = max(4, .Options$digits - 3),
                                nchar.trts = x$nchar.trts,
                                details = FALSE,
                                details.disconnected = x$details.disconnected,
                                ...) {
  
  chkclass(x, "netconnection")
  ##
  if (is.null(nchar.trts))
    nchar.trts <- 666
  
  
  chknumeric(digits, length = 1)
  chknumeric(nchar.trts, min = 1, length = 1)
  chklogical(details)
  details.disconnected <- replaceNULL(details.disconnected, FALSE)
  chklogical(details.disconnected)
  
  
  matitle(x)
  ##
  cat("Number of studies: k = ", x$k, "\n", sep = "")
  cat("Number of pairwise comparisons: m = ", x$m, "\n", sep = "")
  cat("Number of treatments: n = ", x$n, "\n", sep = "")
  if (!is.null(x$d))
    cat("Number of designs: d = ", x$d, "\n", sep = "")
  ##
  cat("Number of subnetworks: ", x$n.subnets, "\n", sep = "")
  ##
  if (x$n.subnets > 1) {
    f <- function(x) length(unique(x))
    d <- as.data.frame(x)
    k.subset <- tapply(d$studlab, d$subnet, f)
    ##
    m <- as.matrix(
      data.frame(subnetwork = names(k.subset),
                 k = as.vector(k.subset),
                 m = as.vector(tapply(d$studlab, d$subnet, length)),
                 n = as.vector(tapply(c(d$treat1, d$treat2),
                                      c(d$subnet, d$subnet), f)))
    )
    rownames(m) <- rep("", nrow(m))
    ##
    cat("\nDetails on subnetworks: \n")
    prmatrix(m, quote = FALSE, right = TRUE)
    ##
    if (details.disconnected) {
      cat("\n")
      for (i in seq_len(x$n.subnets)) {
        d.i <- subset(d, d$subnet == i)
        cat("Subnetwork ", i, ":\n", sep = "")
        print(sort(unique(c(d.i$treat1, d.i$treat2))))
      }
    }
  }
  
  
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
    
    diff.rownames <- rownames(x$D.matrix) != rownames(D)
    if (any(diff.rownames)) {
      abbr <- rownames(D)
      full <- rownames(x$D.matrix)
      ##
      tmat <- data.frame(abbr, full)
      names(tmat) <- c("Abbreviation", "Treatment name")
      tmat <- tmat[diff.rownames, ]
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat("\nLegend:\n")
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(abbr)))
    }
  }
  
  invisible(NULL)
}





#' @rdname netconnection
#' @export


netconnection <- function(data, ...)
  UseMethod("netconnection")
