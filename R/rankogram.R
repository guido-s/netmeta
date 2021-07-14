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
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
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
#' \item{nsim, comb.fixed, comb.random}{As defined above},
#' \item{small.values, x}{As defined above},
#'
#' @author Theodoros Papakonstantinou \email{dev@tpapak.com}
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
#'
#' plot(ran1)      
#'
#' @rdname rankogram
#' @export rankogram

rankogram <- function(x, nsim = 1000,
                      comb.fixed = x$comb.fixed, comb.random = x$comb.random,
                      small.values = x$small.values,
                      nchar.trts = x$nchar.trts) {
  
  
  meta:::is.installed.package("mvtnorm")
  ##
  meta:::chkclass(x, "netmeta")
  ##
  meta:::chknumeric(nsim, min = 1, length = 1)
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  ##
  small.values <- meta:::setchar(small.values, c("good", "bad"))
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
  }
  

  if (!comb.fixed) {
    ranking.matrix.fixed <- NULL
    ranking.fixed <- NULL
  }
  ##
  if (!comb.random) {
    ranking.matrix.random <- NULL
    ranking.random <- NULL
  }
  
  
  res <- list(ranking.matrix.fixed = ranking.matrix.fixed,
              ranking.fixed = ranking.fixed,
              ranking.matrix.random = ranking.matrix.random,
              ranking.random = ranking.random,
              nsim = nsim,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              small.values = small.values,
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
                            nchar.trts = x$nchar.trts,
                            digits = gs("digits.prop"),
                            ...) {
  
  
  meta:::chkclass(x, "rankogram")
  ##  
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  ##
  if (is.null(nchar.trts))
    nchar.trts <- 666
  meta:::chknumeric(nchar.trts, length = 1)
  ##
  meta:::chknumeric(digits, length = 1)
  

  if (!comb.fixed & !comb.random)
    return(invisible(NULL))


  cat(paste0("Rankogram (based on ", x$nsim,
             " simulation", if (x$nsim > 1) "s", ")\n\n"))
  ##
  if (comb.fixed) {
    rownames(x$ranking.matrix.fixed) <-
      treats(x$ranking.matrix.fixed, nchar.trts)
    ##
    cat("Fixed effects model: \n\n")
    prmatrix(meta:::formatN(x$ranking.matrix.fixed, digits),
             quote = FALSE, right = TRUE, ...)
    if (comb.random)
      cat("\n")
  }
  ##
  if (comb.random) {
    rownames(x$ranking.matrix.random) <-
      treats(x$ranking.matrix.random, nchar.trts)
    ##
    cat("Random effects model: \n\n")
    prmatrix(meta:::formatN(x$ranking.matrix.random, digits),
             quote = FALSE, right = TRUE, ...)
  }
  
  
  invisible(NULL)
}
