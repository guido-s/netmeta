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
#'   meta-analysis.
#' @param verbose ARGUMENT / FUNKTIONALITAET NOTWENDIG ???
#' 
#' @return
#' An object of class \code{"netimpact"} with corresponding
#' \code{netgraph} function. The object is a list containing the
#' following components:
#' \item{impact}{A matrix with contributions of individual studies
#'   (columns) to comparisons (rows).}
#' \item{seTE.ignore, x}{As defined above.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de},
#'   Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netgraph.netimpact}},
#'   \code{\link{print.netimpact}}
#' 
#' @examples
#' data(parkinson)
#' 
#' p1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'                n = list(n1, n2, n3),
#'                mean = list(y1, y2, y3),
#'                sd = list(sd1, sd2, sd3),
#'                data = parkinson, studlab = Study)
#' 
#' net1 <- netmeta(p1)
#' ni1 <- netimpact(net1, verbose = TRUE)
#' ni1
#'
#' netgraph(ni1)
#' 
#' @export netimpact


netimpact <- function(x, seTE.ignore = 1e4, verbose = FALSE) {
  
  
  meta:::chkclass(x, "netmeta")
  
  
  studlab <- x$studlab
  treat1  <- x$treat1
  treat2  <- x$treat2
  TE      <- x$TE
  seTE    <- x$seTE
  ##
  comparison <- paste(treat1, sep = x$sep.trts, treat2)
  
  
  comparisons <- names(x$prop.direct.fixed)
  ##  
  impact <- matrix(NA, ncol = x$k, nrow = length(comparisons))
  ##
  rownames(impact) <- comparisons
  colnames(impact) <- x$studies
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
    seTE.i <- seTE
    seTE.i[studlab == i] <- seTE.ignore
    ##
    net.i <- netmeta(TE, seTE.i, treat1, treat2, studlab)
    nets[[i]] <- net.i
    ##
    seTE.nma <- x$seTE.fixed
    seTE.i   <- net.i$seTE.fixed
    ##
    impact.i <- 1 - (lowertri(seTE.nma) / lowertri(seTE.i))^2
    zero <- abs(impact.i) < .Machine$double.eps^0.5
    ##
    impact[, x$studies == i] <- impact.i
    impact[zero, x$studies == i] <- 0
  }
  
  
  res <- list(impact = t(impact),
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
