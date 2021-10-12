bySummary <- function(x, index1, index2,
                      data = NULL, subset = NULL,
                      na.rm = TRUE, ties.method = "random", sep,
                      long = TRUE) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chklogical(na.rm)
  ##
  ties.method <- setchar(ties.method, c("first", "last", "random"))
  ##
  chklogical(long)
  ##
  nulldata <- is.null(data)
  ##
  if (!nulldata & inherits(data, "netmeta")) {
    data <- data$data
    if (is.null(subset) & !is.null(data$.subset))
      subset <- data$.subset
  }
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch 'x', 'index1', and 'index2' from data:
  ##
  x <- eval(mf[[match("x", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  chknull(x)
  k.All <- length(x)
  ##
  index1 <- eval(mf[[match("index1", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  chknull(index1)
  ##
  missing.index2 <- missing(index2)
  if (!missing.index2) {
    index2 <- eval(mf[[match("index2", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
    chknull(index2)
  }
  ##
  ## Catch 'subset' from data:
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  missing.subset <- is.null(subset)
  ##
  if (missing(sep))
    sep <- ":"
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  fun <- "calcSummary"
  ##
  chklength(index1, k.All, fun)
  if (!missing.index2)
    chklength(index2, k.All, fun)
  
  
  ##
  ##
  ## (4) Subset
  ##
  ##
  if (!missing.subset)
    if ((is.logical(subset) & (sum(subset) > k.All)) ||
        (length(subset) > k.All))
      stop("Length of argument 'subset' is larger than number of observations.")
  
  
  ##
  ##
  ## (5) Use subset for analysis
  ##
  ##
  if (!missing.subset) {
    x <- x[subset]
    index1 <- index1[subset]
    if (!missing.index2)
      index2 <- index2[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(x)
  ##
  if (k.all == 0)
    stop("No observations in dataset / subset.")
  
  
  ##
  ##
  ## (6) Determine summaries
  ##
  ##
  if (missing(index2))
    indices <- index1
  else
    indices <- paste(index1, index2, sep = sep)
  ##
  sum    <- by(x, indices, sum, na.rm = na.rm)
  min    <- by(x, indices, min, na.rm = na.rm)
  max    <- by(x, indices, max, na.rm = na.rm)
  mean   <- by(x, indices, mean, na.rm = na.rm)
  median <- by(x, indices, median, na.rm = na.rm)
  ##
  ## Mode
  ##
  tab.x <- table(indices, x)
  max.x <- max.col(tab.x, ties.method = ties.method)
  mode <- colnames(tab.x)[max.x]
  if (is.numeric(x))
    mode <- as.numeric(mode)
  ##
  indices.x <- rownames(sum)
  
  
  ##
  ##
  ## (7) Merge summaries with original data / subset
  ##
  ##
  tdat <- data.frame(indices = indices.x,
                     sum = as.vector(sum),
                     min = as.vector(min),
                     max = as.vector(max),
                     mean = as.vector(mean),
                     median = as.vector(median),
                     mode = mode,
                     stringsAsFactors = FALSE)
  ##
  if (long) {
    res <- data.frame(index1, index2 = index1,
                      indices, x, order = seq(along = index1),
                      stringsAsFactors = FALSE)
    ##
    if (missing(index2))
      res$index2 <- NULL
    else
      res$index2 <- index2
    ##
    res <- merge(res, tdat, by = "indices")
    ##
    res <- res[order(res$order), ]
    ##
    res$indices <- NULL
    ## res$order <- NULL
  }
  else
    res <- tdat
  
  
  res
}
