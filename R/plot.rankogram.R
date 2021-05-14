#' Plot rankograms
#'
#' @description
#' This function produces a rankogram, i.e., an image plot of ranking
#' probabilities for all treatments.
#'
#' @param x An object of class \code{\link{rankogram}}.
#' @param type A character string specifying whether a "bar" chart or
#'   a "line" graph should be drawn.
#' @param ylim The y limits (min, max) of the plot.
#' @param ylab A label for the y-axis.
#' @param \dots Additional graphical arguments (ignored at the
#'   moment).
#' 
#' @details
#' This function produces an image plot of ranking probabilities for
#' all treatments as a bar graph or as a line graph. Treatments are
#' sorted according to their mean effects.
#'
#' @author Theodoros Papakonstantinou
#'   \email{papakonstantinou@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{rankogram}}
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
#' plot(ran1, type = "l")      
#' 
#' @method plot rankogram
#' @export
#' @export plot.rankogram


plot.rankogram <- function(x, type = "bar", ylim, ylab, ...) {
  
  
  meta:::chkclass(x, c("rankogram"))
  ##
  type <- meta:::setchar(type, c("bar", "line"))
  ##
  missing.ylim <- missing(ylim)
  if (!missing.ylim)
    meta:::chknumeric(ylim, min = 0, max = 1, length = 2)
  ##
  if (!missing(ylab))
    meta:::chkchar(ylab, length = 1)
  else
    ylab <- "Probability"
  
  
  mytheme <-
    theme(axis.line.x = ggplot2::element_line(colour = "grey22", size = 1,
                                              linetype = "solid"),
          axis.line.y = ggplot2::element_line(colour = "grey22", size = 1,
                                              linetype = "solid"),
          panel.grid.major  = ggplot2::element_line(color = "grey90"),
          panel.grid.minor  = ggplot2::element_line(color = "grey90"),
          panel.background = ggplot2::element_rect(fill = "grey90"),
          plot.background = ggplot2::element_rect(
                                       fill = "grey90",
                                       colour = "white",
                                       size = 1)
          )
  
  
  if (!is.null(x$ranking.matrix.random)) {
    rankmatrix = "ranking.matrix.random"
    sucras = "ranking.random"
  }
  else {
    rankmatrix = "ranking.matrix.fixed"
    sucras = "ranking.fixed"
  }
  
  
  plotranks <- function(treat) {
    df <- data.frame(pos = 1:nrow(x[[rankmatrix]]),
                     ranks = x[[rankmatrix]][treat, ])
    mymaxvalue <- max(x[[rankmatrix]])

    if (type == "bar")
      p <- ggplot2::ggplot(mapping = aes(df$pos, df$ranks)) +
        ggplot2::geom_col()
    else if (type == "line")
      p <- ggplot2::ggplot(mapping = aes(df$pos, df$ranks)) +
        ggplot2::geom_line()
    ##
    p <- p + ggplot2::scale_x_continuous(breaks = seq(1, nrow(x[[rankmatrix]]), 1))
    p <- p + ggplot2::labs(x = paste("Rank of", treat))
    p <- p + ggplot2::labs(y = ylab)
    ##
    if (missing.ylim)
      p <- p + ggplot2::expand_limits(x = c(1, x$x$n), y = c(0, mymaxvalue))
    else
      p <- p + ggplot2::expand_limits(x = c(1, x$x$n), y = ylim)
    ##
    p + mytheme
  }
  ##
  sortedtreats <- names(sort(x[[sucras]], decreasing = TRUE))
  rankplots <- do.call(gridExtra::grid.arrange, lapply(sortedtreats, plotranks))
}
