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
  meta:::chknumeric(nchar.trts, min = 1, single = TRUE)
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
      warning("No information given for argument 'studlab'. Assuming that comparisons are from independent studies.")
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
                   text = "Arguments 'treat1' and 'treat2' must have the same length.")
  meta:::chklength(studlab, k.All, fun,
                   text = "Arguments 'treat1' and 'studlab' must have the same length.")
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
    stop(paste("Study '", names(tabnarms)[sel.narms],
               "' has a wrong number of comparisons.",
               "\n  Please provide data for all treatment comparisons (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
               sep = ""))
  if (sum(sel.narms) > 1)
    stop(paste("The following studies have a wrong number of comparisons: ",
               paste(paste("'", names(tabnarms)[sel.narms], "'", sep = ""),
                     collapse = ", "),
               "\n  Please provide data for all treatment comparisons (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
               sep = ""))
  
  
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
  
  
  res <- list(treat1 = treat1,
              treat2 = treat2,
              studlab = studlab,
              ##
              k = length(unique(studlab)),
              m = m,
              n = n,
              n.subnets = n.subsets,
              ##
              D.matrix = D,
              A.matrix = A,
              L.matrix = L,
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
