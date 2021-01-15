#' Frequentist method to rank treatments in network
#' 
#' @description
#' Ranking treatments in frequentist network meta-analysis without
#' resampling methods.
#' 
#' @aliases netrank print.netrank
#' 
#' @param x An object of class \code{netmeta} (netrank function) or
#'   \code{netrank} (print function).
#' @param comb.fixed A logical indicating whether to print P-scores
#'   for the fixed effects (common effects) model.
#' @param comb.random A logical indicating whether to print P-scores
#'   for the random effects model.
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"good"}) or
#'   harmful (\code{"bad"}) effect, can be abbreviated.
#' @param sort A logical indicating whether printout should be sorted
#'   by decreasing P-score.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param \dots Additional arguments passed on to
#'   \code{\link{print.data.frame}} function (used internally).
#' 
#' @details
#' Treatments are ranked based on a network meta-analysis. Ranking is
#' performed by P-scores. P-scores are based solely on the point
#' estimates and standard errors of the network estimates. They
#' measure the extent of certainty that a treatment is better than
#' another treatment, averaged over all competing treatments (Rücker
#' and Schwarzer 2015).
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
#' treatment. This interpretation is comparable to that of the Surface
#' Under the Cumulative RAnking curve (SUCRA) which is the rank of
#' treatment \emph{i} within the range of treatments, measured on a
#' scale from 0 (worst) to 1 (best) (Salanti et al. 2011).
#'
#' @return
#' An object of class \code{netrank} with corresponding \code{print}
#' function. The object is a list containing the following components:
#' \item{Pscore.fixed}{A named numeric vector with P-scores for fixed
#'   effects model.}
#' \item{Pmatrix.fixed}{Numeric matrix based on pairwise one-sided
#'   p-values for fixed effects model.}
#' \item{Pscore.random}{A named numeric vector with P-scores for
#'   random effects model.}
#' \item{Pmatrix.random}{Numeric matrix based on pairwise one-sided
#'   p-values of random effects model.}
#' \item{small.values, x}{As defined above.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}
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
#' @rdname netrank
#' @export netrank


netrank <- function(x, small.values = x$small.values) {
  
  ## Check for netmeta object
  ##
  meta:::chkclass(x, c("netmeta", "netcomb"))
  ##
  if (is.null(small.values))
    small.values <- "good"
  small.values <- meta:::setchar(small.values, c("good", "bad"))
  
  
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
  Pscore.fixed <- rowMeans(P.fixed, na.rm = TRUE)
  ##
  if (!all(is.na(TE.random)))
    Pscore.random <- rowMeans(P.random, na.rm = TRUE)
  else
    Pscore.random <- NA
  
  
  ##
  if (!x$comb.random) {
    Pscore <- Pscore.fixed
    Pmatrix <- P.fixed
  }
  else {
    Pscore <- Pscore.random
    Pmatrix <- P.random
  }
  
  
  res <- list(Pscore.fixed = Pscore.fixed,
              Pmatrix.fixed = P.fixed,
              Pscore.random = Pscore.random,
              Pmatrix.random = P.random,
              small.values = small.values,
              x = x,
              title = x$title,
              version = packageDescription("netmeta")$Version,
              Pscore = "'Pscore' replaced by 'Pscore.fixed' and 'Pscore.random'.",
              Pmatrix = "'Pmatrix' replaced by 'Pmatrix.fixed' and 'Pmatrix.random'.")

  class(res) <- "netrank"
  
  res
}





#' @rdname netrank
#' @method print netrank
#' @export
#' @export print.netrank


print.netrank <- function(x,
                          comb.fixed = x$x$comb.fixed,
                          comb.random = x$x$comb.random,
                          sort = TRUE,
                          digits = max(4, .Options$digits - 3),
                          ...) {
  
  meta:::chklogical(comb.fixed)
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
      warning("Argument 'sort=\"random\"' ignored for P-scores of fixed effects model.",
              call. = FALSE)
      sort <- TRUE
    }
    if (comb.random & sort == "fixed") {
      warning("Argument 'sort=\"fixed\"' ignored for P-scores of random effects model.",
              call. = FALSE)
      sort <- TRUE
    }
  }
  else if (both & !is.character(sort) && sort)
    sort <- "random"
  
  
  if (both) {
    if (is.character(sort)) {
      res.both <- data.frame(fixed = round(x$Pscore.fixed, digits),
                             random = round(x$Pscore.random, digits))
      res.both <- res.both[order(-res.both[sort]), ]
    }
    else if (!sort) {
      res.both <- data.frame(fixed = round(x$Pscore.fixed[x$x$seq], digits),
                             random = round(x$Pscore.random[x$x$seq], digits))
    }
    ##
    colnames(res.both)  <- c("P-score (fixed)", "P-score (random)")
  }
  else {
    if (sort) {
      res.fixed <- as.data.frame(round(x$Pscore.fixed[order(-x$Pscore.fixed)],
                                       digits))
      res.random <- as.data.frame(round(x$Pscore.random[order(-x$Pscore.random)],
                                        digits))
    }
    else {
      res.fixed <- as.data.frame(round(x$Pscore.fixed[x$x$seq], digits))
      res.random <- as.data.frame(round(x$Pscore.random[x$x$seq], digits))
    }
    ##
    colnames(res.fixed)  <- "P-score"
    colnames(res.random) <- "P-score"
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
  
  invisible(NULL)
}
