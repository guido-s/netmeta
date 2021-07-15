#' Calculate rankogram
#'
#' @description
#' This function calculates the probabilities of each treatment being
#' at each possible rank and the SUCRAs (Surface Under the Cumulative
#' RAnking curve) in frequentist network meta-analysis.
#'
#' @param x An object of class \code{\link{netmeta}}.
#' @param nsim Number of simulations.
#' @param comb.fixed A logical indicating to compute ranking
#'   probabilities and SUCRAs for the fixed effects (common effects)
#'   model.
#' @param comb.random A logical indicating to compute ranking
#'   probabilities and SUCRAs for the random effects model.
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"good"}) or
#'   harmful (\code{"bad"}) effect, can be abbreviated.
#' @param cumulative.rankprob A logical indicating whether cumulative
#'   ranking probabilites should be printed.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param \dots Additional arguments for printing.
#'
#' @details
#' We derive a matrix showing the probability of each treatment being
#' at each possible rank. To this aim, we use resampling from a
#' multivariate normal distribution with estimated network effects as
#' means and corresponding estimated variance covariance matrix. We
#' then summarise them using the ranking metric SUCRAs (Surface Under
#' the Cumulative RAnking curve).
#'
#' @return
#' An object of class \code{rankogram} with corresponding \code{print}
#' and \code{plot} function. The object is a list containing the
#' following components:
#' \item{ranking.matrix.fixed}{Numeric matrix giving the probability
#'   of each treatment being at each possible rank for the fixed
#'   effects model.}
#' \item{ranking.fixed}{SUCRA values for the fixed effects model.}
#' \item{ranking.matrix.random}{Numeric matrix giving the probability
#'   of each treatment being at each possible rank for the random
#'   effects model.}
#' \item{ranking.random}{SUCRA values for the random effects model.}
#' \item{cumrank.matrix.fixed}{Numeric matrix giving the cumulative
#'   ranking probability of each treatment for the
#'   fixed effects model.}
#' \item{cumrank.matrix.random}{Numeric matrix giving the cumulative
#'   ranking probability of each treatment for the random effects
#'   model.}
#' \item{nsim, comb.fixed, comb.random}{As defined above},
#' \item{small.values, x}{As defined above},
#'
#' @author Theodoros Papakonstantinou \email{dev@@tpapak.com}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#'
#' @seealso \code{\link{netmeta}}, \code{\link{netrank}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP (2011):
#' Graphical methods and numerical summaries for presenting results
#' from multiple-treatment meta-analysis: an overview and tutorial.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{64}, 163--71
#'
#' @examples
#' data(Woods2010)
#' p1 <- pairwise(treatment, event = r, n = N, studlab = author,
#'                data = Woods2010, sm = "OR")
#' net1 <- netmeta(p1, small.values = "good")
#'
#' ran1 <- rankogram(net1, nsim = 100)
#' ran1
#' print(ran1, cumulative.rankprob = TRUE)
#'
#' plot(ran1)
#'
#' @rdname rankogram
#' @export rankogram

rankogram <- function(x, nsim = 1000,
                      comb.fixed = x$comb.fixed, comb.random = x$comb.random,
                      small.values = x$small.values,
                      cumulative.rankprob = FALSE,
                      nchar.trts = x$nchar.trts) {
  
  
  meta:::is.installed.package("mvtnorm")
  ##
  meta:::chkclass(x, "netmeta")
  ##
  chklogical <- meta:::chklogical
  meta:::chknumeric(nsim, min = 1, length = 1)
  chklogical(comb.fixed)
  chklogical(comb.random)
  ##
  small.values <- meta:::setchar(small.values, c("good", "bad"))
  ##
  chklogical(cumulative.rankprob)
  ##
  if (is.null(nchar.trts))
    nchar.trts <- 666
  meta:::chknumeric(nchar.trts, length = 1)
  
  
  resampling <- lapply(1:nsim,
                       function(y) matrix(NA,
                                          nrow = nrow(x$TE.random),
                                          ncol = ncol(x$TE.random),
                                          dimnames = list(rownames(x$TE.random),
                                                          colnames(x$TE.random))))
  
  rearrange <- function(y, simul) {
    resampling[[y]][lower.tri(resampling[[y]])] <-  simul[y, ]
    resampling[[y]] <- t(resampling[[y]])
    resampling[[y]][lower.tri(resampling[[y]])] <- -simul[y, ]
    resampling
  }
  
  
  if (comb.fixed) {
    simul.fixed <- mvtnorm::rmvnorm(nsim,
                                    t(x$TE.fixed)[lower.tri(x$TE.fixed)],
                                    x$Cov.fixed)
    
    resampling.fixed <- lapply(1:nsim,
                               function(y)
                                 rearrange(y, simul = simul.fixed)[[y]])
    
    if(small.values == "good") {
      rankings.fixed <- lapply(1:nsim,
                                function(y)
                                  rank(rowSums(resampling.fixed[[y]]> 0,
                                               na.rm = TRUE)))
      sortedtreats.fixed <- order(x$TE.fixed[, 1])
    }
    else if (small.values == "bad") {
      rankings.fixed <- lapply(1:nsim,
                               function(y)
                               (x$n + 1) -
                               rank(rowSums(resampling.fixed[[y]] > 0,
                                            na.rm = TRUE)))
      sortedtreats.fixed <- order(x$TE.fixed[1, ])
    }
    
    df.fixed <-
      data.frame(x = mapply(function(i){
        (rankings.fixed[[i]][sortedtreats.fixed])}, 1:nsim))
    
    ## Ranking matrix
    ##
    ranking.matrix.fixed <-
      mapply(function(x){
        mapply(function(y) {
          sum((df.fixed[x, ] == y) / nsim)
        },
        1:length(rownames(df.fixed)))
      }, 1:length(rownames(df.fixed)))
    ##
    rownames(ranking.matrix.fixed) <- rownames(df.fixed)
    colnames(ranking.matrix.fixed) <- 1:x$n
    
    ## Cumulative ranking matrix
    ##
    rank.cum.fixed <- t(mapply(function(i)
      cumsum(ranking.matrix.fixed[i,]), 1:nrow(ranking.matrix.fixed)))
    ##
    dimnames(rank.cum.fixed) <- list(rownames(ranking.matrix.fixed), 1:x$n)
    
    ## SUCRAs
    ##
    ranking.fixed <- mapply(function(i)
      sum(rank.cum.fixed[i, -length(rank.cum.fixed[i, ])]) /
      (nrow(rank.cum.fixed) - 1), 
      1:nrow(ranking.matrix.fixed))
    ##
    names(ranking.fixed) <- rownames(ranking.matrix.fixed)
    ##
    ranking.fixed <- ranking.fixed[x$trts]
  }
  
  
  if (comb.random) {
    simul.random <- mvtnorm::rmvnorm(nsim,
                                     t(x$TE.random)[lower.tri(x$TE.random)],
                                     x$Cov.random)
    
    resampling.random <- lapply(1:nsim,
                                function(y)
                                  rearrange(y, simul = simul.random)[[y]])
    
    if(small.values == "good") {
      rankings.random <- lapply(1:nsim,
                                function(y)
                                  rank(rowSums(resampling.random[[y]]> 0,
                                               na.rm = TRUE)))
      sortedtreats.random <- order(x$TE.random[, 1])
    }
    else if (small.values == "bad") {
      rankings.random <- lapply(1:nsim,
                                function(y)
                                (x$n + 1) -
                                rank(rowSums(resampling.random[[y]] > 0,
                                             na.rm = TRUE)))
      sortedtreats.random <- order(x$TE.random[1, ])
    }
    
    df.random <-
      data.frame(x = mapply(function(i){
        (rankings.random[[i]][sortedtreats.random])}, 1:nsim))
    
    ## Ranking matrix
    ##
    ranking.matrix.random <-
      mapply(function(x){
        mapply(function(y) {
          sum((df.random[x, ] == y) / nsim)
        },
        1:length(rownames(df.random)))
      }, 1:length(rownames(df.random)))
    ##
    rownames(ranking.matrix.random) <- rownames(df.random)
    colnames(ranking.matrix.random) <- 1:x$n
    
    ## Cumulative ranking matrix
    ##
    rank.cum.random <- t(mapply(function(i)
      cumsum(ranking.matrix.random[i,]), 1:nrow(ranking.matrix.random)))
    ##
    dimnames(rank.cum.random) <- list(rownames(ranking.matrix.random), 1:x$n)
    
    ## SUCRAs
    ##
    ranking.random <- mapply(function(i)
      sum(rank.cum.random[i, -length(rank.cum.random[i, ])]) /
      (nrow(rank.cum.random) - 1), 
      1:nrow(ranking.matrix.random))
    ##
    names(ranking.random) <- rownames(ranking.matrix.random)
    ##
    ranking.random <- ranking.random[x$trts]
    ##
    cumrank.matrix.random <- t(apply(ranking.matrix.random, 1, cumsum))
  }
  

  if (!comb.fixed) {
    ranking.matrix.fixed <- NULL
    ranking.fixed <- NULL
    rank.cum.fixed <- NULL
  }
  ##
  if (!comb.random) {
    ranking.matrix.random <- NULL
    ranking.random <- NULL
    rank.cum.random <- NULL
  }
  
  
  res <- list(ranking.matrix.fixed = ranking.matrix.fixed,
              ranking.fixed = ranking.fixed,
              ranking.matrix.random = ranking.matrix.random,
              ranking.random = ranking.random,
              ##
              cumrank.matrix.fixed = rank.cum.fixed,
              cumrank.matrix.random = rank.cum.random,
              ##
              nsim = nsim,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              small.values = small.values,
              cumulative.rankprob = cumulative.rankprob,
              nchar.trts = nchar.trts,
              x = x)
  
  
  class(res) <- "rankogram"
  ##
  res
}





#' @rdname rankogram
#' @method print rankogram
#' @export
#' @export print.rankogram


print.rankogram <- function(x,
                            comb.fixed = x$comb.fixed,
                            comb.random = x$comb.random,
                            cumulative.rankprob = x$cumulative.rankprob,
                            nchar.trts = x$nchar.trts,
                            digits = gs("digits.prop"),
                            legend = TRUE,
                            ...) {
  
  
  meta:::chkclass(x, "rankogram")
  ##
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  chklogical(comb.fixed)
  chklogical(comb.random)
  if (is.null(cumulative.rankprob))
    cumulative.rankprob <- FALSE
  chklogical(cumulative.rankprob)
  ##
  if (is.null(nchar.trts))
    nchar.trts <- 666
  chknumeric(nchar.trts, length = 1)
  ##
  chknumeric(digits, length = 1)
  chklogical(legend)
  
  
  if (!comb.fixed & !comb.random)
    return(invisible(NULL))


  cat(paste0(if (cumulative.rankprob)
               "Cumulative ranking probabilities" else "Rankogram",
             " (based on ", x$nsim, " simulation",
             if (x$nsim > 1) "s", ")\n\n"))
  ##
  if (comb.fixed) {
    if (cumulative.rankprob)
      fixed <- x$cumrank.matrix.fixed
    else
      fixed <- x$ranking.matrix.fixed
    rownames(fixed) <- treats(fixed, nchar.trts)
    ##
    cat("Fixed effects model: \n\n")
    prmatrix(meta:::formatN(fixed, digits), quote = FALSE, right = TRUE, ...)
    if (comb.random)
      cat("\n")
  }
  ##
  if (comb.random) {
    if (cumulative.rankprob)
      random <- x$cumrank.matrix.random
    else
      random <- x$ranking.matrix.random
    rownames(random) <-
      treats(random, nchar.trts)
    ##
    cat("Random effects model: \n\n")
    prmatrix(meta:::formatN(random, digits), quote = FALSE, right = TRUE, ...)
  }
  
  
  if (legend) {
    if (comb.fixed)
      trts <- rownames(x$ranking.matrix.fixed)
    else if (comb.random)
      trts <- rownames(x$ranking.matrix.random)
    ##
    trts.abbr <- treats(trts, nchar.trts)
    if (any(trts != trts.abbr)) {
      cat("\nLegend:\n")
      tmat <- data.frame(trts.abbr, trts)
      names(tmat) <- c("Abbreviation", "Treatment name")
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(trts.abbr)))
    }
  }
  
  
  invisible(NULL)
}
