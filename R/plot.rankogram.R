#' Rankograms
#'
#' @description
#' Produce rankograms, an image plot of ranking probabilities for all treatments generated with R function rankogram.
#' 
#'
#' @param x An object of class \code{\link{rankogram}}.
#' @param type A character string specifying whether a "bar" chart or a "line" graph should be drawn.
#' 
#' @details 
#' This function produces an image plot of ranking probabilities for all treatments as a bar graph or as a line graph. 
#' Treatments are sorted according to their mean effects. 
#'
#' @return 
#' The rankogram plot for all treatments \code{plotranks}
#' 
#' #' @seealso \code{\link{netmeta}} \code{\link{netrank}} \code{\link{rankogram}}
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
#' plot(ran1,type="bar")
#' 
#' @method plot rankogram
#' @export
#' @export plot.rankogram

plot.rankogram = function(x,type="bar") 
{

meta:::chkclass(x, c("rankogram"))

mytheme=theme(axis.line.x = ggplot2::element_line(colour = "grey22", size = 1, linetype = "solid"),
                axis.line.y = ggplot2::element_line(colour = "grey22", size = 1, linetype = "solid"),
                panel.grid.major  = ggplot2::element_line(color = "grey90"),
                panel.grid.minor  = ggplot2::element_line(color = "grey90"),
                panel.background = ggplot2::element_rect(fill = "grey90"),
                plot.background = ggplot2::element_rect(
                fill = "grey90",
                colour = "white",
                size = 1)
 )

  if(!is.null(x$ranking.matrix.random)){
    rankmatrix = "ranking.matrix.random"
    sucras = "SUCRA.random"
  }else{
    rankmatrix = "ranking.matrix.fixed"
    sucras = "SUCRA.fixed"
  }

plotranks=function(treat){
  df=data.frame(pos=1:nrow(x[[rankmatrix]]), ranks=x[[rankmatrix]][treat,])
  mymaxvalue=max(x[[rankmatrix]])
  
  if(type=="bar"){
    p=ggplot2::ggplot(df, aes(pos, ranks)) +
      ggplot2::geom_col()
  }
  if(type=="line"){
    p=ggplot2::ggplot(df, aes(pos, ranks)) +
      ggplot2::geom_line()
  }
  p=p + ggplot2::scale_x_continuous(breaks = seq(1,nrow(x[[rankmatrix]]),1))
  p=p + ggplot2::labs(x = paste("Rank of", treat))
  p=p + ggplot2::labs(y = "Probability")
  p=p + ggplot2::expand_limits(x=c(1,x$n),y=c(0,mymaxvalue))
  p+mytheme
}

sortedtreats=names(sort(x[[sucras]], decreasing=T))


rankplots=do.call(gridExtra::grid.arrange, lapply(sortedtreats,plotranks))

}



