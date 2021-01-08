#' Determine the importance of individual studies in network
#' meta-analysis
#' 
#' @description
#' This function measures the importance of individual studies in
#' network meta-analysis by the reduction of the precision if the
#' study is removed / ignored from the network.
#' 
#' @param x An object of class \code{netmeta}.
#' @param seTE.ignore Assumed (large) standard error in order to
#'   mimicking the removal of individual studies from the network
#'   meta-analysis (ignored for \code{\link{netmetabin}} objects).
#' @param event.ignore Assumed event number mimicking the removal of
#'   individual studies from the network meta-analysis (considered for
#'   \code{\link{netmetabin}} objects).
#' @param verbose A logical indicating whether information on the
#'   estimation progress should be printed.
#' 
#' @return
#' An object of class \code{"netimpact"} with corresponding
#' \code{netgraph} and \code{print} function. The object is a list
#' containing the following components:
#' \item{impact.fixed}{A matrix with contributions of individual
#'   studies (columns) to comparisons (rows) under the fixed effects
#'   model.}
#' \item{impact.random}{A matrix with contributions of individual
#'   studies (columns) to comparisons (rows) under the random effects
#'   model.}
#' \item{ignored.comparisons}{List with comparisons of ignored study.}
#' \item{seTE.ignore, event.ignore, x}{As defined above.}
#' \item{nets}{List of all network meta-analyses (removing a single
#'   study).}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de},
#'   Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netmetabin}},
#'   \code{\link{netgraph.netimpact}}, \code{\link{print.netimpact}}
#' 
#' @examples
#' data(parkinson)
#' 
#' # Only consider first four studies (to reduce runtime of example)
#' #
#' p1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'                n = list(n1, n2, n3),
#'                mean = list(y1, y2, y3),
#'                sd = list(sd1, sd2, sd3),
#'                data = subset(parkinson, Study < 5),
#'                studlab = Study)
#' 
#' net1 <- netmeta(p1)
#' ni1 <- netimpact(net1, verbose = TRUE)
#' ni1
#'
#' netgraph(ni1)
#' 
#' @export netimpact


netimpact <- function(x, seTE.ignore = 1e4, event.ignore = 0.01, verbose = FALSE) {
  
  
  meta:::chkclass(x, "netmeta")
  ##
  x <- upgradenetmeta(x)
  ##
  meta:::chknumeric(seTE.ignore, min = 0, zero = TRUE, length = 1)
  meta:::chknumeric(event.ignore, min = 0, zero = TRUE, length = 1)
  
  
  studlab <- x$studlab
  treat1  <- x$treat1
  treat2  <- x$treat2
  TE      <- x$TE
  seTE    <- x$seTE
  ##
  event1 <- x$event1
  event2 <- x$event2
  n1 <- x$n1
  n2 <- x$n2
  ##
  comparison <- paste(treat1, sep = x$sep.trts, treat2)
  
  
  comparisons <- rownames(x$Cov.fixed)
  ##
  impact.fixed <- impact.random <-
    matrix(NA, ncol = length(x$studies), nrow = length(comparisons))
  ##
  rownames(impact.fixed) <- rownames(impact.random) <- comparisons
  colnames(impact.fixed) <- colnames(impact.random) <- x$studies
  ##
  ## Run network meta-analyses "excluding" one study at a time
  ##
  ignored <- nets <- list()
  ##
  for (i in x$studies) {
    ##
    ignored[[i]] <- comparison[studlab == i]
    ##
    if (verbose) {
      cat(paste0("** Removed study: ", i,
                 " **\nComparison",
                 if (length(comparison[studlab == i]) > 1) "s",
                 ": ",
                 paste(paste("'", comparison[studlab == i],
                             "'", sep = ""),
                       collapse = ", "),
                 "\n\n"))
    }
    ##
    if (!inherits(x, "netmetabin")) {
      seTE.i <- seTE
      seTE.i[studlab == i] <- seTE.ignore
      ##
      net.i <- netmeta(TE, seTE.i, treat1, treat2, studlab, tau.preset = x$tau)
    }
    else {
      event1.i <- event1
      event2.i <- event2
      n1.i <- n1
      n2.i <- n2
      event1.i[studlab == i] <- event.ignore
      event2.i[studlab == i] <- event.ignore
      n1.i[studlab == i] <- 1
      n2.i[studlab == i] <- 1
      ##
      net.i <- netmetabin(event1.i, n1.i, event2.i, n2.i, treat1, treat2, studlab,
                          method = x$method, sm = x$sm)
    }
    nets[[i]] <- net.i
    ##
    seTE.fixed <- x$seTE.fixed
    seTE.fixed.i <- net.i$seTE.fixed
    ##
    seTE.random <- x$seTE.random
    seTE.random.i <- net.i$seTE.random
    ##
    impact.fixed.i <- 1 - (lowertri(seTE.fixed) / lowertri(seTE.fixed.i))^2
    zero.fixed <- abs(impact.fixed.i) < .Machine$double.eps^0.5
    ##
    impact.fixed[, x$studies == i] <- impact.fixed.i
    impact.fixed[zero.fixed, x$studies == i] <- 0
    ##
    impact.random.i <- 1 - (lowertri(seTE.random) / lowertri(seTE.random.i))^2
    zero.random <- abs(impact.random.i) < .Machine$double.eps^0.5
    ##
    impact.random[, x$studies == i] <- impact.random.i
    impact.random[zero.random, x$studies == i] <- 0
  }
  
  
  res <- list(impact.fixed = t(impact.fixed),
              impact.random = t(impact.random),
              ignored.comparisons = ignored,
              seTE.ignore = seTE.ignore,
              x = x,
              nets = nets,
              version = packageDescription("netmeta")$Version)
  ##
  class(res) <- "netimpact"
  ##
  res
}
