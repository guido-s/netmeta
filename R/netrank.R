#' Frequentist method to rank treatments in network
#' 
#' @description
#' Ranking treatments in frequentist network meta-analysis without
#' resampling methods.
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
#' @param comb.fixed A logical indicating whether to print P-scores or
#'   SUCRAs for the fixed effects (common effects) model.
#' @param comb.random A logical indicating whether to print P-scores
#'   or SUCRAs for the random effects model.
#' @param sort A logical indicating whether printout should be sorted
#'   by decreasing P-score.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
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
#' \item{ranking.fixed}{A named numeric vector with P-scores or SUCRAs
#'   for the fixed effects model.}
#' \item{Pmatrix.fixed}{Numeric matrix based on pairwise one-sided
#'   p-values for the fixed effects model.}
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
#'                 data = Senn2013, sm = "MD",
#'                 comb.random = FALSE)
#' 
#' nr1 <- netrank(net1)
#' nr1
#' print(nr1, sort = FALSE)
#'
#' \dontrun{
#' net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'                 data = Senn2013, sm = "MD")
#' 
#' nr2 <- netrank(net2)
#' nr2
#' print(nr2, sort = "fixed")
#' print(nr2, sort = FALSE)
#' }
#' 
#' \dontrun{
#' net3 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'                 data = Senn2013, sm = "MD")
#' 
#' nr3 <- netrank(net3, method = "SUCRA", nsim = 100)
#' nr3
#' print(nr3, sort = "fixed")
#' print(nr3, sort = FALSE)
#' }
#' 
#' @rdname netrank
#' @export netrank


netrank <- function(x, small.values = x$small.values, method, nsim,
                    comb.fixed = x$comb.fixed, comb.random = x$comb.random) {
  
  ## Check for netmeta object
  ##
  meta:::chkclass(x, c("netmeta", "netcomb", "rankogram"))
  
  
  if (missing(method))
    if (inherits(x, c("netmeta", "netcomb")))
      method <- "P-score"
    else
      method <- "SUCRA"
  else
    method <- meta:::setchar(method, c("P-score", "SUCRA"))
  ##
  if (is.null(small.values))
    small.values <- "good"
  else
    small.values <- meta:::setchar(small.values, c("good", "bad"))
  ##
  if (is.null(comb.fixed))
    comb.fixed <- TRUE
  meta:::chklogical(comb.fixed)
  if (is.null(comb.random))
    comb.random <- TRUE
  meta:::chklogical(comb.random)
  
  
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
                       comb.fixed = comb.fixed, comb.random = comb.random)
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
    P.fixed <- NULL
    P.random <- NULL
    ##
    ranking.fixed <- rnk$ranking.fixed
    ranking.random <- rnk$ranking.random
  }
  else {
    ##
    ## P-scores
    ##
    nsim <- NULL
    ##
    TE.fixed <- x$TE.fixed
    pval.fixed <- x$pval.fixed
    ##
    TE.random <- x$TE.random
    pval.random <- x$pval.random
    
    
    ## Calculate one-sided p-values
    ##
    w.fixed <- (1 + sign(TE.fixed)) / 2
    p.fixed <- pval.fixed
    ##
    if (small.values == "good")
      P.fixed <- w.fixed * p.fixed / 2 + (1 - w.fixed) * (1 - p.fixed / 2)
    else
      P.fixed <- w.fixed * (1 - p.fixed / 2) + (1 - w.fixed) * p.fixed / 2
    ##
    w.random <- (1 + sign(TE.random)) / 2
    p.random <- pval.random
    ##
    if (small.values == "good")
      P.random <- w.random * p.random / 2 + (1 - w.random) * (1 - p.random / 2)
    else
      P.random <- w.random * (1 - p.random / 2) + (1 - w.random) * p.random / 2
    
    
    ## Row means provide P-scores
    ##
    ranking.fixed <- rowMeans(P.fixed, na.rm = TRUE)
    ##
    if (!all(is.na(TE.random)))
      ranking.random <- rowMeans(P.random, na.rm = TRUE)
    else
      ranking.random <- NA
  }
  
  
  res <- list(ranking.fixed = ranking.fixed,
              Pmatrix.fixed = P.fixed,
              ##
              ranking.random = ranking.random,
              Pmatrix.random = P.random,
              ##
              method = method,
              small.values = small.values,
              ##
              nsim = nsim,
              ##
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              ##
              x = x,
              title = x$title,
              version = packageDescription("netmeta")$Version,
              ##
              Pscore = paste("'Pscore' replaced by 'ranking.fixed'",
                             "and 'ranking.random'."),
              Pmatrix = paste("'Pmatrix' replaced by 'Pmatrix.fixed'",
                              "and 'Pmatrix.random'."),
              ##
              ## Backward compatibility (for R package NMAoutlier)
              ##
              Pscore.fixed = if (method == "P-score") ranking.fixed else NULL,
              Pscore.random = if (method == "P-score") ranking.random else NULL
              )

  class(res) <- "netrank"
  
  res
}





#' @rdname netrank
#' @method print netrank
#' @export
#' @export print.netrank


print.netrank <- function(x,
                          comb.fixed = x$comb.fixed,
                          comb.random = x$comb.random,
                          sort = TRUE,
                          digits = max(4, .Options$digits - 3),
                          ...) {
  
  
  meta:::chkclass(x, "netrank")
  ##
  if (is.null(comb.fixed))
    comb.fixed <- !is.null(x$ranking.fixed)
  meta:::chklogical(comb.fixed)
  if (is.null(comb.random))
    comb.random <- !is.null(x$ranking.random)
  meta:::chklogical(comb.random)
  ##
  if (is.character(sort))
    sort <- meta:::setchar(sort, c("fixed", "random"))
  else
    meta:::chklogical(sort)
  ##
  meta:::chknumeric(digits, length = 1)
  
  
  both <- (comb.fixed + comb.random) == 2
  ##
  if (!both & is.character(sort)) {
    if (comb.fixed & sort == "random") {
      warning("Argument 'sort=\"random\"' ignored for fixed effects model.",
              call. = FALSE)
      sort <- TRUE
    }
    if (comb.random & sort == "fixed") {
      warning("Argument 'sort=\"fixed\"' ignored for random effects model.",
              call. = FALSE)
      sort <- TRUE
    }
  }
  else if (both & !is.character(sort) && sort)
    sort <- "random"


  if (is.null(x$method)) {
    x$method <- "P-score"
    x$ranking.fixed <- x$Pscore.fixed
    x$ranking.random <- x$Pscore.random
  }
  
  
  if (both) {
    if (is.character(sort)) {
      res.both <- data.frame(fixed = round(x$ranking.fixed, digits),
                             random = round(x$ranking.random, digits))
      res.both <- res.both[order(-res.both[sort]), ]
    }
    else if (!sort) {
      res.both <- data.frame(fixed = round(x$ranking.fixed[x$x$seq], digits),
                             random = round(x$ranking.random[x$x$seq], digits))
    }
    ##
    colnames(res.both)  <- paste(x$method, c("(fixed)", "(random)"))
  }
  else {
    if (sort) {
      if (comb.fixed)
        res.fixed <-
          as.data.frame(round(x$ranking.fixed[order(-x$ranking.fixed)],
                              digits))
      if (comb.random)
        res.random <-
          as.data.frame(round(x$ranking.random[order(-x$ranking.random)],
                              digits))
    }
    else {
      if (comb.fixed)
        res.fixed <- as.data.frame(round(x$ranking.fixed[x$x$seq], digits))
      if (comb.random)
        res.random <- as.data.frame(round(x$ranking.random[x$x$seq], digits))
    }
    ##
    if (comb.fixed)
      colnames(res.fixed)  <- x$method
    if (comb.random)
      colnames(res.random) <- x$method
  }
  ##
  matitle(x)
  ##
  if (both)
    prmatrix(res.both, quote = FALSE, ...)
  ##
  else if (comb.fixed) {
    prmatrix(res.fixed, quote = FALSE, ...)
    if (comb.random)
      cat("\n")
  }
  ##
  else if (comb.random) {
    prmatrix(res.random, quote = FALSE, ...)
  }

  if (x$method == "SUCRA" & !is.null(x$nsim))
    cat(paste0("\n- based on ", x$nsim,
               " simulation", if (x$nsim > 1) "s", "\n"))
  
  
  invisible(NULL)
}
