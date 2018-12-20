netmatrix <- function(x, var, levels, labels = levels,
                      func = "mode", ties.method = "random") {
  
  
  meta:::chkclass(x, "netmeta")
  ##
  setchar <- meta:::setchar
  ##
  func <- setchar(func,
                  c("mode", "min", "max", "mean", "median", "sum"))
  ties.method <- setchar(ties.method, c("first", "last", "random"))
  
  
  ##
  ## Catch var, treat1, treat2 from data:
  ##
  mf <- match.call()
  data <- x$data
  ##
  var <- eval(mf[[match("var", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  treat1 <- data[[".treat1"]]
  treat2 <- data[[".treat2"]]
  ##
  if (length(var) == 1 & length(treat1) > 1)
    var <- rep(var, length(treat1))
  ##
  subset <- data$.subset
  excl <- data$.excl
  drop <- data$.drop
  ##
  if (is.null(subset))
    subset <- rep(TRUE, nrow(data))
  if (is.null(excl))
    excl <- rep(FALSE, nrow(data))
  if (is.null(drop))
    drop <- rep(FALSE, nrow(data))
  ##
  treat1 <- treat1[subset & !excl & !drop]
  treat2 <- treat2[subset & !excl & !drop]
  var <- var[subset & !excl & !drop]
  
  
  dat <- bySummary(var, treat1, treat2,
                   ties.method = ties.method,
                   sep = x$sep.trts, long = FALSE)
  ##
  selfirst <- function(x) x[1]
  selsecond <- function(x) x[2]
  split <- strsplit(dat$indices, x$sep.trts)
  dat$idx1 <- unlist(lapply(split, selfirst))
  dat$idx2 <- unlist(lapply(split, selsecond))
  ##
  dat <- dat[, c("idx1", "idx2", func)]
  ##
  if (missing(levels))
    levels <- sort(unique(dat[[func]]))
  ##
  if (length(levels) != length(labels))
    stop("Different lengths of arguments 'levels' and 'labels'.")
  ##
  dat[[func]] <- as.character(factor(dat[[func]],
                                     levels = levels,
                                     labels = labels))
  ##
  if (is.numeric(labels))
    dat[[func]] <- as.numeric(dat[[func]])
  
  
  res <- x$A.matrix
  res[!is.na(res)] <- NA
  ##
  for (i in seq(along = dat$idx1)) {
    res[dat$idx1[i], dat$idx2[i]] <- dat[[func]][i]
    res[dat$idx2[i], dat$idx1[i]] <- dat[[func]][i]
  }
  ##
  res
}
