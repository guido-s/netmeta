#' Measures for characterizing a network meta-analysis
#' 
#' @description
#' This function provides measures for quantifying the direct evidence
#' proportion, the mean path length and the minimal parallelism (the
#' latter on aggregated and study level) of mixed treatment
#' comparisons (network estimates) as well as the evidence flow per
#' design, see König et al. (2013).  These measures support the
#' critical evaluation of the network meta-analysis results by
#' rendering transparent the process of data pooling.
#' 
#' @param x An object of class \code{netmeta}.
#' @param random A logical indicating whether random effects model
#'   should be used to calculate network measures.
#' @param tau.preset An optional value for the square-root of the
#'   between-study variance \eqn{\tau^2}.
#' @param warn A logical indicating whether warnings should be
#'   printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @details
#' The direct evidence proportion gives the absolute contribution of
#' direct effect estimates combined for two-arm and multi-arm studies
#' to one network estimate.
#' 
#' Concerning indirectness, comparisons with a mean path length beyond
#' two should be interpreted with particular caution, as more than two
#' direct comparisons have to be combined serially on average.
#' 
#' Large indices of parallelism, either on study-level or on
#' aggregated level, can be considered as supporting the validity of a
#' network meta-analysis if there is only a small amount of
#' heterogeneity.
#' 
#' The network estimates for two treatments are linear combinations of
#' direct effect estimates comparing these or other treatments. The
#' linear coefficients can be seen as the generalization of weights
#' known from classical meta-analysis. These coefficients are given in
#' the projection matrix \eqn{H} of the underlying model. For
#' multi-arm studies, the coefficients depend on the choice of the
#' study-specific baseline treatment, but the absolute flow of
#' evidence can be made explicit for each design as shown in König et
#' al. (2013) and is given in \code{H.tilde}.
#' 
#' All measures are calculated based on the common effects
#' meta-analysis by default. In the case that in function
#' \code{netmeta} the argument \code{random = TRUE}, all measures
#' are calculated for a random effects model. The value of the
#' square-root of the between-study variance \eqn{\tau^2} can also be
#' prespecified by argument \code{tau.preset} in function
#' \code{netmeta}.
#'
#' @return
#' A list containing the following components:
#' \item{random, tau.preset}{As defined above.}
#' \item{proportion}{A named vector of the direct evidence proportion
#'   of each network estimate.}
#' \item{meanpath}{A named vector of the mean path length of each
#'   network estimate.}
#' \item{minpar}{A named vector of the minimal parallelism on
#'   aggregated level of each network estimate.}
#' \item{minpar.study}{A named vector of the minimal parallelism on
#'   study level of each network estimate.}
#' \item{H.tilde }{Design-based hat matrix with information on
#'   absolute evidence flow per design. The number of rows is equal to
#'   the number of possible pairwise treatment comparisons and the
#'   number of columns is equal to the number of designs.}
#' 
#' @author Ulrike Krahn \email{ulrike.krahn@@bayer.com}, Jochem König
#'   \email{koenigjo@@uni-mainz.de}
#' 
#' @seealso \link{netmeta}
#' 
#' @references
#' König J, Krahn U, Binder H (2013):
#' Visualizing the flow of evidence in network meta-analysis and
#' characterizing mixed treatment comparisons.
#' \emph{Statistics in Medicine},
#' \bold{32}, 5414--29
#' 
#' @examples
#' data(Senn2013)
#' 
#' # Conduct common effects network meta-analysis with reference
#' # treatment 'plac', i.e. placebo
#' #
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD", reference = "plac", random = FALSE)
#' 
#' # Calculate measures based on a common effects model
#' #        
#' nm1 <- netmeasures(net1)
#' 
#' # Plot of minimal parallelism versus mean path length
#' #
#' plot(nm1$meanpath, nm1$minpar, pch = "",
#'   xlab = "Mean path length", ylab = "Minimal parallelism")
#' text(nm1$meanpath, nm1$minpar, names(nm1$meanpath), cex = 0.8)
#'
#' \dontrun{
#' # Conduct random effects network meta-analysis with reference
#' # treatment 'plac', i.e. placebo
#' #
#' net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD", reference = "plac", common = FALSE)
#' 
#' # Calculate measures based on a random effects model
#' #                          
#' nm2 <- netmeasures(net2)
#' }
#' 
#' @export netmeasures


netmeasures <- function(x,
                        random = x$random | !missing(tau.preset),
                        tau.preset = x$tau.preset,
                        warn = TRUE, warn.deprecated = gs("warn.deprecated"),
                        ...) {
  
  ##
  ##
  ## (1) Check for netmeta object and upgrade object
  ##
  ##
  chkclass(x, "netmeta")
  x <- updateversion(x)
  ##
  is.bin <- inherits(x, "netmetabin")


  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  random <-
    deprecated(random, missing(random), args, "comb.random", warn.deprecated)
  chklogical(random)
  ##
  if (is.bin & random) {
    txt <-
      paste0("Argument 'random' set to FALSE for ",
             if (x$method == "MH")
               "Mantel-Haenszel method"
             else
               "non-central hypergeometric distribution",
             ".")
    if (warn)
      warning(txt, call. = FALSE)
    random <- FALSE
  }
  ##
  if (!missing(random) & !random) {
    if (!is.null(tau.preset)) {
      if (!missing(tau.preset) & tau.preset > 0)
        stop("Argument 'tau.preset' must be equal 0 if random=FALSE.")
      tau.preset <- NULL
    }
  }
  ##
  if (!is.null(tau.preset)) {
    chknumeric(tau.preset, min = 0, length = 1)
    if (!random & warn) {
      warning("Measures calculated for random effects model ",
              "(argument random=TRUE) as argument 'tau.preset' is provided.")
    }
  }
  
  
  ##
  ##
  ## (3) Calculate measures
  ##
  ##
  if (is.bin) {
    ##
    ## Direct evidence proportion
    ## (Rücker et al. 2020, BMC Med Res Meth, equation 7)
    ##
    se.direct  <- uppertri(x$seTE.direct.common)
    se.network <- uppertri(x$seTE.common)
    proportion <- se.network^2 / se.direct^2
    proportion[is.na(proportion)] <- 0
    names(proportion) <- rownames(x$Cov.common)
    ##
    meanpath <- NA
    minpar <- NA
    minpar.study <- NA
    H.tilde <- NA
  }
  else {
    x$reference.group <- ""
    ##
    if (random == FALSE & length(tau.preset) == 0) {
      nmak <- nma.krahn(x)
      if (nmak$n == 2) {
        prop <- mpath <- 1
        names(prop) <- names(mpath) <- x$designs
        ##
        res <- list(proportion = prop,
                    meanpath = mpath,
                    minpar = NULL,
                    minpar.study = NULL,
                    H.tilde = NULL,
                    random = random,
                    tau.preset = tau.preset)
        ##
        return(res)
      }
    }
    ##                                                             
    if (length(tau.preset) == 1) {
      nmak <- nma.krahn(x, tau.preset = tau.preset)
      if (nmak$n == 2) {
        prop <- mpath <- 1
        names(prop) <- names(mpath) <- x$designs
        ##
        res <- list(proportion = prop,
                    meanpath = mpath,
                    minpar = NULL,
                    minpar.study = NULL,
                    H.tilde = NULL,
                    random = random,
                    tau.preset = tau.preset)
        ##
        return(res)
      }
    }
    ##
    if (random == TRUE & length(tau.preset) == 0) {
      nmak <- nma.krahn(x, tau.preset = x$tau)
      if (nmak$n == 2) {
        prop <- mpath <- 1
        names(prop) <- names(mpath) <- x$designs
        ##
        res <- list(proportion = prop,
                    meanpath = mpath,
                    minpar = NULL,
                    minpar.study = NULL,
                    H.tilde = NULL,
                    random = random,
                    tau.preset = tau.preset)
        ##
        return(res)
      }
    }
    ##
    comparisons <- nmak$direct[, "comparison"]
    direct  <- nmak$direct
    network <- nmak$network
    ##
    H <- nmak$H
    H.studies <- nmak$H.studies
    ##
    ## Direct evidence proportion (Koenig 2013, subsection 3.4.1 and 3.4.2)
    ##
    proportion <- rep(0, nrow(H))
    names(proportion) <- rownames(H)
    proportion[comparisons] <- network[comparisons,"seTE"]^2 / direct$seTE^2
    ##
    ## Minimal parallelism (Koenig 2013, subsection 3.4.3)
    ## Mean path length (Koenig 2013, subsection 3.4.4)
    ##
    l <- lapply(split(H,
                      matrix(colnames(H), nrow = nrow(H),
                             ncol = length(colnames(H)), byrow = TRUE)),
                function(x) matrix(x, nrow = nrow(H))
                )
    ##
    H.tilde <- sapply(l,
                      function(x)
                        apply(x, 1,
                              function(x)
                                0.5 * (abs(sum(x)) + sum(abs(x)))
                              )
                      )
    ##
    rownames(H.tilde) <- rownames(H)
    ##  
    minpar  <- apply(H.tilde, 1, function(x) 1 / max(abs(x)))
    meanpath <- apply(H.tilde, 1, sum)
    ##
    ## Minimal parallelism (study-level)
    ##
    l.studies <- lapply(split(H.studies,
                              matrix(colnames(H.studies), nrow = nrow(H.studies),
                                     ncol = length(colnames(H.studies)), byrow = TRUE)),
                        function(x) matrix(x, nrow = nrow(H.studies))
                        )
    ##
    H.tilde.studies <- sapply(l.studies,
                              function(x)
                                apply(x, 1,
                                      function(x)
                                        0.5 * (abs(sum(x)) + sum(abs(x)))
                                      )
                              )
    ##
    rownames(H.tilde.studies) <- rownames(H)
    ##  
    minpar.study <- apply(H.tilde.studies, 1,
                          function(x)
                            1 / max(abs(x))
                          )
  }
  
  
  ##
  ##
  ## (3) Create netmeasures object
  ##
  ##
  res <- list(proportion = proportion,
              meanpath = meanpath,
              minpar = minpar,
              minpar.study = minpar.study,
              H.tilde = H.tilde,
              random = random,
              tau.preset = tau.preset,
              version = packageDescription("netmeta")$Version
              )
  ##
  res
}
