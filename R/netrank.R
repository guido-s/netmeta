#' Frequentist method to rank treatments in network
#' 
#' @description
#' Ranking treatments in frequentist network meta-analysis with and
#' without resampling methods.
#' 
#' @aliases netrank print.netrank
#' 
#' @param x An object of class \code{netmeta} or \code{rankogram}.
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"good"}) or
#'   harmful (\code{"bad"}) effect, can be abbreviated.
#' @param method A character string specifying whether the
#'   \code{"P-score"} or \code{"SUCRA"} ranking metric will be
#'   calculated.
#' @param nsim Number of simulations to calculate SUCRAs.
#' @param common A logical indicating whether to print P-scores or
#'   SUCRAs for the common effects model.
#' @param random A logical indicating whether to print P-scores or
#'   SUCRAs for the random effects model.
#' @param sort A logical indicating whether printout should be sorted
#'   by decreasing P-score.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments passed on to
#'   \code{\link{print.data.frame}} function (used internally).
#' 
#' @details
#' 
#' Treatments are ranked based on a network meta-analysis. Ranking is
#' performed by a ranking metric: P-score or SUCRA.
#'
#' P-scores are based solely on the point estimates and standard
#' errors of the network estimates. They measure the extent of
#' certainty that a treatment is better than another treatment,
#' averaged over all competing treatments (Rücker and Schwarzer 2015).
#'
#' The Surface Under the Cumulative RAnking curve (SUCRA) is the rank
#' of treatment \emph{i} within the range of treatments, measured on a
#' scale from 0 (worst) to 1 (best) (Salanti et al. 2011). A
#' resampling method is used to calculate SUCRAs for frequentist
#' network meta-analysis. The number of simulations is determine by
#' argument \code{nsim}.
#'
#' The interpretation of P-scores and SUCRAs is comparable.
#' 
#' The P-score of treatment \emph{i} is defined as the mean of all 1 -
#' P[\emph{j}] where P[\emph{j}] denotes the one-sided P-value of
#' accepting the alternative hypothesis that treatment \emph{i} is
#' better than one of the competing treatments \emph{j}. Thus, if
#' treatment \emph{i} is better than many other treatments, many of
#' these P-values will be small and the P-score will be large. Vice
#' versa, if treatment \emph{i} is worse than most other treatments,
#' the P-score is small.
#' 
#' The P-score of treatment \emph{i} can be interpreted as the mean
#' extent of certainty that treatment \emph{i} is better than another
#' treatment.
#'
#' @return
#' An object of class \code{netrank} with corresponding \code{print}
#' function. The object is a list containing the following components:
#' \item{ranking.common}{A named numeric vector with P-scores or SUCRAs
#'   for the common effects model.}
#' \item{Pmatrix.common}{Numeric matrix based on pairwise one-sided
#'   p-values for the common effects model.}
#' \item{ranking.random}{A named numeric vector with P-scores or
#'   SUCRAs for the random effects model.}
#' \item{Pmatrix.random}{Numeric matrix based on pairwise one-sided
#'   p-values of the random effects model.}
#' \item{small.values, method, x}{As defined above.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}, Theodoros
#'   Papakonstantinou \email{dev@@tpapak.com}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{rankogram}}
#' 
#' @references
#' Rücker G, Schwarzer G (2017):
#' Resolve conflicting rankings of outcomes in network meta-analysis:
#' Partial ordering of treatments.
#' \emph{Research Synthesis Methods},
#' \bold{8}, 526--36
#' 
#' Salanti G, Ades AE, Ioannidis JP (2011):
#' Graphical methods and numerical summaries for presenting results
#' from multiple-treatment meta-analysis: an overview and tutorial.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{64}, 163--71
#' 
#' @examples
#' data(Senn2013)
#' 
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD", random = FALSE)
#' 
#' nr1 <- netrank(net1)
#' nr1
#' print(nr1, sort = FALSE)
#'
#' \dontrun{
#' net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD")
#' 
#' nr2 <- netrank(net2)
#' nr2
#' print(nr2, sort = "common")
#' print(nr2, sort = FALSE)
#' }
#' 
#' \dontrun{
#' net3 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD")
#' 
#' nr3 <- netrank(net3, method = "SUCRA", nsim = 100)
#' nr3
#' print(nr3, sort = "common")
#' print(nr3, sort = FALSE)
#' }
#' 
#' @rdname netrank
#' @export netrank


netrank <- function(x, small.values = x$small.values, method, nsim,
                    common = x$common, random = x$random,
                    warn.deprecated = gs("warn.deprecated"),
                    ...) {
  
  ##
  ##
  ## (1) Check and upgrade object
  ##
  ##
  chkclass(x, c("netmeta", "netcomb", "rankogram"))
  x <- updateversion(x)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  if (missing(method))
    if (inherits(x, c("netmeta", "netcomb")))
      method <- "P-score"
    else
      method <- "SUCRA"
  else
    method <- setchar(method, c("P-score", "SUCRA"))
  ##
  if (is.null(small.values))
    small.values <- "good"
  else
    small.values <- setchar(small.values, c("good", "bad"))
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  
  
  if (method == "SUCRA") {
    ##
    ## SUCRAs
    ##
    if (inherits(x, "netcomb"))
      stop("netcomb object is not compatible with SUCRAs.",
           call. = FALSE)
    else if (inherits(x, "netmeta")) {
      if (missing(nsim))
        nsim <- 1000
      ##
      rnk <- rankogram(x,
                       nsim = nsim,
                       small.values = small.values,
                       common = common, random = random)
    }
    else {
      if (!missing(nsim))
        warning("Argument 'nsim' ignored for rankogram object.",
                call. = FALSE)
      rnk <- x
      x <- rnk$x
      nsim <- rnk$nsim
    }
    ##
    P.common <- NULL
    P.random <- NULL
    ##
    ranking.common <- rnk$ranking.common
    ranking.random <- rnk$ranking.random
  }
  else {
    ##
    ## P-scores
    ##
    nsim <- NULL
    ##
    TE.common <- x$TE.common
    pval.common <- x$pval.common
    ##
    TE.random <- x$TE.random
    pval.random <- x$pval.random
    
    
    ## Calculate one-sided p-values
    ##
    w.common <- (1 + sign(TE.common)) / 2
    p.common <- pval.common
    ##
    if (small.values == "good")
      P.common <-
        w.common * p.common / 2 + (1 - w.common) * (1 - p.common / 2)
    else
      P.common <-
        w.common * (1 - p.common / 2) + (1 - w.common) * p.common / 2
    ##
    w.random <- (1 + sign(TE.random)) / 2
    p.random <- pval.random
    ##
    if (small.values == "good")
      P.random <-
        w.random * p.random / 2 + (1 - w.random) * (1 - p.random / 2)
    else
      P.random <-
        w.random * (1 - p.random / 2) + (1 - w.random) * p.random / 2
    
    
    ## Row means provide P-scores
    ##
    ranking.common <- rowMeans(P.common, na.rm = TRUE)
    ##
    if (!all(is.na(TE.random)))
      ranking.random <- rowMeans(P.random, na.rm = TRUE)
    else
      ranking.random <- NA
  }
  
  
  res <- list(ranking.common = ranking.common,
              Pmatrix.common = P.common,
              ##
              ranking.random = ranking.random,
              Pmatrix.random = P.random,
              ##
              method = method,
              small.values = small.values,
              ##
              nsim = nsim,
              ##
              common = common,
              random = random,
              ##
              x = x,
              title = x$title,
              version = packageDescription("netmeta")$Version,
              ##
              Pscore = paste("'Pscore' replaced by 'ranking.common'",
                             "and 'ranking.random'."),
              Pmatrix = paste("'Pmatrix' replaced by 'Pmatrix.common'",
                              "and 'Pmatrix.random'.")
              )
  ##
  ## Backward compatibility
  ##
  res$fixed <- common
  ##
  res$ranking.fixed <- res$ranking.common
  res$Pmatrix.fixed <- res$Pmatrix.common
  ##
  ## For R package NMAoutlier:
  ##
  res$Pscore.fixed <-
    if (method == "P-score") ranking.common else NULL
  res$Pscore.random <-
    if (method == "P-score") ranking.random else NULL
  
  
  class(res) <- "netrank"
  
  res
}





#' @rdname netrank
#' @method print netrank
#' @export


print.netrank <- function(x,
                          common = x$common,
                          random = x$random,
                          sort = TRUE,
                          digits = gs("digits.prop"),
                          warn.deprecated = gs("warn.deprecated"),
                          ...) {
  
  ##
  ##
  ## (1) Check netrank object and upgrade object
  ##
  ##
  chkclass(x, "netrank")
  x <- updateversion(x)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  if (is.character(sort)) {
    sort <- setchar(sort, c("common", "random", "fixed"))
    sort[sort == "fixed"] <- "common"
  }
  else
    chklogical(sort)
  ##
  chknumeric(digits, length = 1)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  
  
  both <- (common + random) == 2
  ##
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


  if (is.null(x$method)) {
    x$method <- "P-score"
    x$ranking.common <- x$Pscore.common
    x$ranking.random <- x$Pscore.random
  }
  
  
  if (both) {
    if (is.character(sort)) {
      res.both <- data.frame(common = round(x$ranking.common, digits),
                             random = round(x$ranking.random, digits))
      res.both <- res.both[order(-res.both[, sort]), ]
    }
    else if (!sort) {
      res.both <- data.frame(common = round(x$ranking.common[x$x$seq], digits),
                             random = round(x$ranking.random[x$x$seq], digits))
    }
    ##
    colnames(res.both) <-
      paste(x$method,
            paste0("(", c(gs("text.w.common"), gs("text.w.random")), ")"))
  }
  else {
    if (sort) {
      if (common)
        res.common <-
          as.data.frame(round(x$ranking.common[order(-x$ranking.common)],
                              digits))
      if (random)
        res.random <-
          as.data.frame(round(x$ranking.random[order(-x$ranking.random)],
                              digits))
    }
    else {
      if (common)
        res.common <- as.data.frame(round(x$ranking.common[x$x$seq], digits))
      if (random)
        res.random <- as.data.frame(round(x$ranking.random[x$x$seq], digits))
    }
    ##
    if (common)
      colnames(res.common)  <- x$method
    if (random)
      colnames(res.random) <- x$method
  }
  ##
  matitle(x)
  ##
  if (both)
    prmatrix(res.both, quote = FALSE, ...)
  ##
  else if (common) {
    prmatrix(res.common, quote = FALSE, ...)
    if (random)
      cat("\n")
  }
  ##
  else if (random) {
    prmatrix(res.random, quote = FALSE, ...)
  }

  if (x$method == "SUCRA" & !is.null(x$nsim))
    cat(paste0("\n- based on ", x$nsim,
               " simulation", if (x$nsim > 1) "s", "\n"))
  
  
  invisible(NULL)
}
