#' Rankogram
#'
#' @description
#' The function \code{\link{rankograms}} gives the probabilities of each treatment
#' being at each possible rank and the SUCRAs (Surface Under the Cumulative RAnking curve)
#' in frequentist network meta-analysis. 
#'
#' @param x An object of class \code{\link{netmeta}}.
#' @param comb.fixed A logical indicating to compute ranking probabilities and SUCRAs for the 
#' fixed effect (common effect) model.
#' @param comb.random A logical indicating to compute ranking probabilities and SUCRAs for the 
#' random effects model.
#' @param nsim number of simulations.
#' @param small.values A character string specifying whether small treatment effects indicate a beneficial (\code{"good"}) 
#' or harmful (\code{"bad"}) effect, can be abbreviated.
#'
#' @details 
#' We derive a matrix showing the probability of each treatment being at 
#' each possible rank. To this aim, we use resampling from a multivariate normal 
#' distribution with means the estimated network effects and variance covariance matrix their 
#' estimated variance covariance matrix. We then summarise them using the ranking metric
#' SUCRAs (Surface Under the Cumulative RAnking curve).
#'
#' @return 
#' \item{ranking.matrix.fixed}{Numeric matrix giving the probability of each treatment being 
#' at each possible rank for the fixed effect model.}
#' \item{SUCRA.fixed}{SUCRA values for the fixed effect model.}
#' \item{ranking.matrix.random}{Numeric matrix giving the probability of each treatment being 
#' at each possible rank for the random effects model.}
#' \item{SUCRA.random}{SUCRA values for the random effects model.}
#'
#' @seealso \code{\link{netmeta}} \code{\link{netrank}}
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
#' p1 <- pairwise(treatment, event = r, n = N, studlab = author, data = Woods2010, sm = "OR")
#' net1 <- netmeta(p1, small.values="good")
#' 
#' ran1=rankogram(x=net1, comb.fixed=T, comb.random=T, 
#'                 nsim=1000, small.values="good")
#'                 
#'                 
#' @rdname rankogram
#' @export rankogram

rankogram = function( x
                    , comb.fixed=x$comb.fixed
                    , comb.random=x$comb.random
                    , nsim=1000
                    , small.values=x$small.values)
{

  resampling <- lapply(1:nsim,function(y) matrix(NA,nrow=nrow(x$TE.random),ncol=ncol(x$TE.random),
                       dimnames = list(rownames(x$TE.random), colnames(x$TE.random))))
  
  rearrange <- function(y,simul){
    resampling[[y]][lower.tri(resampling[[y]])] <- simul[y,]
    resampling[[y]] <- t(resampling[[y]])
    resampling[[y]][lower.tri(resampling[[y]])] <- -simul[y,]
    resampling
  }
  
  if (comb.random) {
    #random effects
    simul.random <- mvtnorm::rmvnorm(nsim, t(x$TE.random)[lower.tri(x$TE.random)], x$Cov.random)
    
    resampling.random <-lapply(1:nsim,function(y) rearrange(y,simul = simul.random)[[y]])
    
    if(small.values=="good"){rankings.random<-lapply(1:nsim,function(y) rank(rowSums(resampling.random[[y]]> 0, na.rm = T)))}
    if(small.values=="bad"){rankings.random<-lapply(1:nsim,function(y) (x$n + 1)-rank(rowSums(resampling.random[[y]]> 0, na.rm = T)))}
    
    if(small.values=="good") {sortedtreats.random=order(x$TE.random[,1])}
    if(small.values=="bad") {sortedtreats.random=order(x$TE.random[1,])}
    
    df.random=data.frame(x=mapply(function(i){(rankings.random[[i]][sortedtreats.random])},1:nsim))
    
    #ranking matrix
    ranking.matrix.random=mapply(function(x){
      mapply(function(y)
      {sum((df.random[x,]==y)/nsim)},
      1:length(rownames(df.random)))},1:length(rownames(df.random)))
    
    rownames(ranking.matrix.random) <- rownames(df.random)
    colnames(ranking.matrix.random) <- 1:x$n
    
    #cumulative ranking matrix
    rank.cum.random=t(mapply(function(i) cumsum(ranking.matrix.random[i,]), 
                             1:nrow(ranking.matrix.random)))
    dimnames(rank.cum.random) = list(rownames(ranking.matrix.random), 1:x$n)
    
    #SUCRAs
    SUCRA.random=mapply(function(i) sum(rank.cum.random[i,-length(rank.cum.random[i,])])/(nrow(rank.cum.random)-1), 
                        1:nrow(ranking.matrix.random))
    names(SUCRA.random)=rownames(ranking.matrix.random)
  }
  
  if (comb.fixed) {
     #fixed effects
     simul.fixed <- mvtnorm::rmvnorm(nsim, t(x$TE.fixed)[lower.tri(x$TE.fixed)], x$Cov.fixed)
     
     resampling.fixed <-lapply(1:nsim,function(y) rearrange(y,simul = simul.fixed)[[y]])
     
     if(small.values=="good"){rankings.fixed<-lapply(1:nsim,function(y) rank(rowSums(resampling.fixed[[y]]> 0, na.rm = T)))}
     if(small.values=="bad"){rankings.fixed<-lapply(1:nsim,function(y) (x$n + 1)-rank(rowSums(resampling.fixed[[y]]> 0, na.rm = T)))}
     
     if(small.values=="good") {sortedtreats.fixed=order(x$TE.fixed[,1])}
     if(small.values=="bad") {sortedtreats.fixed=order(x$TE.fixed[1,])}
     
     df.fixed=data.frame(x=mapply(function(i){(rankings.fixed[[i]][sortedtreats.fixed])},1:nsim))
     
     #ranking matrix
     ranking.matrix.fixed=mapply(function(x){
       mapply(function(y)
       {sum((df.fixed[x,]==y)/nsim)},
       1:length(rownames(df.fixed)))},1:length(rownames(df.fixed)))
     
     rownames(ranking.matrix.fixed) <- rownames(df.fixed)
     colnames(ranking.matrix.fixed) <- 1:x$n
     
     #cumulative ranking matrix
     rank.cum.fixed=t(mapply(function(i) cumsum(ranking.matrix.fixed[i,]), 
                             1:nrow(ranking.matrix.fixed)))
     dimnames(rank.cum.fixed) = list(rownames(ranking.matrix.fixed), 1:x$n)
     
     #SUCRAs
     SUCRA.fixed=mapply(function(i) sum(rank.cum.fixed[i,-length(rank.cum.fixed[i,])])/(nrow(rank.cum.fixed)-1), 
                        1:nrow(rank.cum.fixed))
     names(SUCRA.fixed)=rownames(ranking.matrix.fixed)
   }

  ranking.matrix.fixed = if (comb.fixed) ranking.matrix.fixed else NULL
  SUCRA.fixed = if (comb.fixed) SUCRA.fixed else NULL
  ranking.matrix.random = if (comb.random) ranking.matrix.random else NULL
  SUCRA.random = if (comb.random) SUCRA.random else NULL
  
  res <- list(ranking.matrix.fixed = ranking.matrix.fixed,
              SUCRA.fixed = SUCRA.fixed,
              ranking.matrix.random = ranking.matrix.random,
              SUCRA.random = SUCRA.random,
              x = x)
     
  class(res) <- "rankogram"
     
  res

}
