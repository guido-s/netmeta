#' Create a matrix with additional information for pairwise
#' comparisons
#' 
#' @description
#' Auxiliary function to create a matrix with additional information
#' for pairwise comparisons
#' 
#' @param x A \code{\link{netmeta}} object.
#' @param var Variable with additional information.
#' @param levels An optional vector of the values that \code{var}
#'   might have taken (see \code{\link{factor}}).
#' @param labels An optional vector with labels for \code{var} (see
#'   \code{\link{factor}}).
#' @param func A character string with the function name to summarize
#'   values within pairwise comparisons; see Details.
#' @param ties.method A character string describing how ties are
#'   handled if \code{func = "mode"}; see Details.
#' 
#' @details
#' For each pairwise comparison, unique values will be calculated for
#' the variable \code{var} based on the argument \code{func}: "mode"
#' (most common value), "min" (minimum value), "max", "mean",
#' "median", and "sum". In order to determine the most common value,
#' the argument \code{ties.method} can be used in the case of ties
#' with "first" meaning that the first / smallest value will be
#' selected; similar for "last" (last / largest value) and "random"
#' (random selection).
#'
#' @return
#' A matrix with the same row and column names as the adjacency matrix
#' \code{x$A.matrix}.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netgraph.netmeta}}
#' 
#' @examples
#' data(smokingcessation)
#' # Add variable with (fictious) risk of bias values
#' # with 1 = "low risk" and 2 = "high risk"
#' #
#' smokingcessation$rob <- rep(1:2, 12)
#' 
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' net1 <- netmeta(p1, common = FALSE, ref = "A")
#' 
#' # Generate network graph with information on risk of bias
#' #
#' col.rob <- netmatrix(net1, rob, ties.method = "last",
#'   levels = 1:2, labels = c("green", "yellow"))
#' #
#' netgraph(net1, plastic = FALSE, col = col.rob,
#'   cex.points = 5, bg.points = "gray", adj = 0.5)
#' 
#' netgraph(net1, plastic = FALSE, col = col.rob,
#'   cex.points = n.trts, bg.points = "blue",
#'   labels = paste0(trts, " (n=", n.trts, ")"),
#'   offset = c(0.05, 0.035, 0.05, 0.025))
#' 
#' @export netmatrix


netmatrix <- function(x, var, levels, labels = levels,
                      func = "mode", ties.method = "random") {
  
  
  chkclass(x, "netmeta")
  ##
  func <- setchar(func,
                  c("mode", "min", "max", "mean", "median", "sum"))
  ties.method <- setchar(ties.method, c("first", "last", "random"))
  
  
  ##
  ## Catch var, treat1, treat2 from data:
  ##
  data <- x$data
  ##
  var <- catch("var", match.call(), data, sys.frame(sys.parent()))
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
  attr(res, "version") <- packageDescription("netmeta")$Version
  ##
  res
}
