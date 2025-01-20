#' Determine the importance of individual studies in network
#' meta-analysis
#' 
#' @description
#' This function measures the importance of individual studies in
#' network meta-analysis by the reduction of the precision if the
#' study is removed / ignored from the network (Rücker et al., 2020).
#' 
#' @param x An object of class \code{netmeta}.
#' @param seTE.ignore Assumed (large) standard error in order to
#'   mimicking the removal of individual studies from the network
#'   meta-analysis (ignored for \code{\link{netmetabin}} objects).
#' @param event.ignore Assumed event number mimicking the removal of
#'   individual studies from the network meta-analysis (considered for
#'   \code{\link{netmetabin}} objects).
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
#' @param nchar.studlab A numeric defining the minimum number of
#'   characters used to create unique study labels.
#' @param verbose A logical indicating whether information on the
#'   estimation progress should be printed.
#' 
#' @return
#' An object of class \code{"netimpact"} with corresponding
#' \code{netgraph} and \code{print} function. The object is a list
#' containing the following components:
#' \item{impact.common}{A matrix with contributions of individual
#'   studies (columns) to comparisons (rows) under the common effects
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
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de},
#'   Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}
#' 
#' @references
#' Rücker G, Nikolakopoulou A, Papakonstantinou T, Salanti G, Riley
#' RD, Schwarzer G (2020):
#' The statistical importance of a study for a network meta-analysis
#' estimate.
#' \emph{BMC Medical Research Methodology},
#' \bold{20}, 190
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netmetabin}},
#'   \code{\link{netgraph.netimpact}}, \code{\link{print.netimpact}},
#'   \code{\link[metadat]{dat.franchini2012}}
#' 
#' @examples
#' # Only consider first two studies (to reduce runtime of example)
#' #
#' studies <- unique(dat.franchini2012$Study)
#' p1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'   n = list(n1, n2, n3),
#'   mean = list(y1, y2, y3), sd = list(sd1, sd2, sd3),
#'   data = subset(dat.franchini2012, Study %in% studies[1:2]),
#'   studlab = Study)
#' 
#' net1 <- netmeta(p1)
#' ni1 <- netimpact(net1, verbose = TRUE)
#' ni1
#'
#' netgraph(ni1)
#' 
#' @export netimpact

netimpact <- function(x,
                      seTE.ignore = 100 * max(x$seTE, na.rm = TRUE),
                      event.ignore = 0.01,
                      nchar.trts = x$nchar.trts,
                      nchar.studlab = x$nchar.studlab,
                      verbose = FALSE) {
  
  
  chkclass(x, "netmeta")
  ##
  x <- updateversion(x)
  ##
  chknumeric(seTE.ignore, min = 0, zero = TRUE, length = 1)
  chknumeric(event.ignore, min = 0, zero = TRUE, length = 1)
  #
  chknumeric(nchar.trts, min = 1, length = 1)
  chknumeric(nchar.studlab, min = 1, length = 1)
  
  
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
  
  
  comparisons <- rownames(x$Cov.common)
  ##
  impact.common <- impact.random <-
    matrix(NA, ncol = length(x$studies), nrow = length(comparisons))
  ##
  rownames(impact.common) <- rownames(impact.random) <- comparisons
  colnames(impact.common) <- colnames(impact.random) <- x$studies
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
      cat("** Removed study: ", i,
          " **\nComparison",
          if (length(comparison[studlab == i]) > 1) "s",
          ": ",
          paste(paste0("'", comparison[studlab == i], "'"), collapse = ", "),
          "\n\n",
          sep = "")
    }
    ##
    if (!inherits(x, "netmetabin")) {
      seTE.i <- seTE
      seTE.i[studlab == i] <- seTE.ignore
      ##
      net.i <- netmeta(TE, seTE.i, treat1, treat2, studlab,
                       tau.preset = x$tau,
                       tol.multiarm = x$tol.multiarm,
                       tol.multiarm.se = x$tol.multiarm.se)
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
      net.i <- netmetabin(event1.i, n1.i, event2.i, n2.i,
                          treat1, treat2, studlab,
                          method = x$method, sm = x$sm)
    }
    nets[[i]] <- net.i
    ##
    seTE.common <- x$seTE.common
    seTE.common.i <- net.i$seTE.common
    ##
    seTE.random <- x$seTE.random
    seTE.random.i <- net.i$seTE.random
    ##
    impact.common.i <- 1 - (lowertri(seTE.common) / lowertri(seTE.common.i))^2
    zero.common <- abs(impact.common.i) < .Machine$double.eps^0.5
    ##
    impact.common[, x$studies == i] <- impact.common.i
    impact.common[zero.common, x$studies == i] <- 0
    ##
    impact.random.i <- 1 - (lowertri(seTE.random) / lowertri(seTE.random.i))^2
    zero.random <- abs(impact.random.i) < .Machine$double.eps^0.5
    ##
    impact.random[, x$studies == i] <- impact.random.i
    impact.random[zero.random, x$studies == i] <- 0
  }
  
  
  res <- list(impact.common = t(impact.common),
              impact.random = t(impact.random),
              ignored.comparisons = ignored,
              seTE.ignore = seTE.ignore,
              #
              x = x,
              nets = nets,
              method.tau = x$method.tau,
              nchar.trts = nchar.trts,
              nchar.studlab = nchar.studlab,
              #
              call = match.call(),
              version = packageDescription("netmeta")$Version)
  ##
  ## Backward compatibility
  ##
  res$impact.fixed <- res$impact.common
  ##
  class(res) <- "netimpact"
  ##
  res
}
