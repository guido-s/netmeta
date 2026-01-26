#' Summary method for objects of class rankogram
#' 
#' @description
#' Summary method for objects of class \code{rankogram} to print
#' list of studies in subnetworks.
#' 
#' @param object An object of class \code{rankogram}.
#' @param x An object of class \code{summary.rankogram}.
#' @param common A logical indicating to print ranking probabilities
#'   and SUCRAs for the common effects model.
#' @param random A logical indicating to print ranking probabilities
#'   and SUCRAs for the random effects model.
#' @param sort A logical indicating whether treatments should be
#'   sorted by decreasing SUCRAs.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param digits Minimal number of significant digits for ranking probabilities
#'   and SUCRAs, see \code{\link{print.default}}.
#' @param digits.mean Minimal number of significant digits for mean ranks, see
#'   \code{\link{print.default}}.
#' @param details.methods A logical specifying whether details on statistical
#'   methods should be printed.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param \dots Additional arguments (passed on to
#'   \code{\link{print.rankogram}}.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' @seealso \code{\link{rankogram}}
#' 
#' @examples
#' # Examples:
#' # example(rankogram.netmeta)
#'
#' @method summary rankogram 
#' @export

summary.rankogram <- function(object, ...) {
  
  chkclass(object, "rankogram")
  
  res <- object
  #
  class(res) <- c("summary.rankogram", class(res))
  #
  res
}


#' @rdname summary.rankogram
#' @method print summary.rankogram 
#' @export

print.summary.rankogram <- function(x,
                                    common = x$common,
                                    random = x$random,
                                    sort = TRUE,
                                    nchar.trts = x$nchar.trts,
                                    digits = gs("digits.prop"),
                                    digits.mean = 2,
                                    details.methods = gs("details"),
                                    legend = gs("legend"),
                                    ...) {
  
  chkclass(x, "summary.rankogram")
  #
  chklogical(common)
  chklogical(random)
  #
  if (is.character(sort)) {
    sort <- setchar(sort, c("common", "random", "fixed"))
    sort[sort == "fixed"] <- "common"
  }
  else
    chklogical(sort)
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.mean, min = 0, length = 1)
  #
  both <- (common + random) == 2
  #
  if (!both & is.character(sort)) {
    if (common & sort == "random") {
      warning("Argument 'sort=\"random\"' ignored for common effects model.",
              call. = FALSE)
      sort <- TRUE
    }
    if (random & sort == "common") {
      warning("Argument 'sort=\"common\"' ignored for random effects model.",
              call. = FALSE)
      sort <- TRUE
    }
  }
  else if (both & !is.character(sort) && sort)
    sort <- "random"
  #
  sort.logical <- if (is.character(sort)) TRUE else sort
  #
  chknumeric(nchar.trts, length = 1)
  chklogical(details.methods)
  chklogical(legend)
  #
  class(x) <- "rankogram"
  #
  print(x,
        common = common,
        random = random,
        sort = sort.logical,
        nchar.trts = nchar.trts,
        digits = digits,
        details.methods = FALSE, legend = FALSE,
        ...)
  #
  if (common | random) {
    cat("\nSurface Under the Cumulative RAnking curve\n\n")
    #
    dat.c <- dat.r <- NULL
    #
    if (common) {
      dat.c <- data.frame(mean = formatN(x$ranking.common, digits))
      names(dat.c) <- paste0("SUCRA", if (random) " (common)")
    }
    #
    if (random) {
      dat.r <- data.frame(mean = formatN(x$ranking.random, digits))
      names(dat.r) <- paste0("SUCRA", if (common) " (random)")
    }
    #
    if (both) {
      dat <- merge(dat.c, dat.r, by = "row.names")
      rownames(dat) <- dat$Row.names
      dat$Row.names <- NULL
      #
      if (sort == "common")
        dat <- dat[rev(order(x$ranking.common)), , drop = FALSE]
      else if (sort == "random")
        dat <- dat[rev(order(x$ranking.random)), , drop = FALSE]
      else
        dat <- dat[x$x$seq, , drop = FALSE]
    }
    else if (common) {
      if (sort.logical)
        dat <- dat.c[rev(order(x$ranking.common)), , drop = FALSE]
      else
        dat <- dat.c[x$x$seq, , drop = FALSE]
    }
    else if (random) {
      if (sort.logical)
        dat <- dat.r[rev(order(x$ranking.random)), , drop = FALSE]
      else
        dat <- dat.r[x$x$seq, , drop = FALSE]
    }
    #
    rownames(dat) <- treats(rownames(dat), nchar.trts)
    #
    prmatrix(dat, quote = FALSE, right = TRUE, ...)
  }
  #
  #
  if (common | random) {
    cat("\nMean ranks\n\n")
    #
    dat.c <- dat.r <- NULL
    #
    if (common) {
      dat.c <- data.frame(mean = formatN(x$meanranks.common, digits.mean))
      names(dat.c) <- paste0("Mean", if (random) " (common)")
    }
    #
    if (random) {
      dat.r <- data.frame(mean = formatN(x$meanranks.random, digits.mean))
      names(dat.r) <- paste0("Mean", if (common) " (random)")
    }
    #
    if (both) {
      dat <- merge(dat.c, dat.r, by = "row.names")
      rownames(dat) <- dat$Row.names
      dat$Row.names <- NULL
      #
      if (sort == "common")
        dat <- dat[rev(order(x$ranking.common)), , drop = FALSE]
      else if (sort == "random")
        dat <- dat[rev(order(x$ranking.random)), , drop = FALSE]
      else
        dat <- dat[x$x$seq, , drop = FALSE]
    }
    else if (common) {
      if (sort.logical)
        dat <- dat.c[rev(order(x$ranking.common)), , drop = FALSE]
      else
        dat <- dat.c[x$x$seq, , drop = FALSE]
    }
    else if (random) {
      if (sort.logical)
        dat <- dat.r[rev(order(x$ranking.random)), , drop = FALSE]
      else
        dat <- dat.r[x$x$seq, , drop = FALSE]
    }
    #
    rownames(dat) <- treats(rownames(dat), nchar.trts)
    #
    prmatrix(dat, quote = FALSE, right = TRUE, ...)
  }
  #
  # Print details of network meta-analysis methods
  #
  if (details.methods) {
    text.details <-
      textmeth(x, random, TRUE)
    #
    cat(text.details)
  }
  #
  # Add legend with abbreviated treatment labels
  #
  if ((common | random) & legend) {
    if (random)
      trts <- rownames(x$ranking.matrix.random)
    else
      trts <- rownames(x$ranking.matrix.common)
    #
    legendabbr(trts, treats(trts, nchar.trts), TRUE)
  }
  #
  invisible(NULL)
}
