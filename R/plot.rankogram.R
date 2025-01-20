#' Plot rankograms
#'
#' @description
#' This function produces a rankogram, i.e., an image plot of ranking
#' probabilities for all treatments.
#'
#' @param x An object of class \code{\link{rankogram}}.
#' @param type A character string specifying whether a "bar" chart, a
#'   "line" graph, or "step" functions should be drawn. Can be
#'   abbreviated.
#' @param pooled A character string indicating whether results for the
#'   common (\code{"common"}) or random effects model (\code{"random"})
#'   should be plotted. Can be abbreviated.
#' @param sort A logical indicating whether treatments should be
#'   sorted by decreasing SUCRAs.
#' @param trts Treatment(s) to show in rankogram.
#' @param cumulative.rankprob A logical indicating whether cumulative
#'   ranking probabilites should be shown.
#' @param ylim The y limits (min, max) of the plot.
#' @param ylab A label for the y-axis.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param \dots Additional graphical arguments (ignored at the
#'   moment).
#' 
#' @details
#' This function produces an image plot of (cumulative) ranking
#' probabilities for all treatments as a bar graph, a line graph or as
#' step functions (argument \code{type}).
#'
#' By default (argument \code{pooled}), results for the random effects
#' model are shown if a network meta-analysis was conducted for both
#' the common and random effects model.
#'
#' Treatments are sorted according to their mean effects if argument
#' \code{sort = TRUE} (default).  A subset of treatments can be
#' specified using argument \code{trts}.
#'
#' Cumulative ranking probabilites are shown if
#' \code{cumulative.rankprob = TRUE}. By default, step functions are
#' shown for cumulative ranking probabilites.
#'
#' @author Theodoros Papakonstantinou \email{dev@@tpapak.com}, Guido
#'   Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{rankogram}}, \code{\link[metadat]{dat.woods2010}}
#'
#' @references 
#' Salanti G, Ades AE, Ioannidis JP (2011):
#' Graphical methods and numerical summaries for presenting results
#' from multiple-treatment meta-analysis: an overview and tutorial.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{64}, 163--71
#' 
#' @examples
#' p1 <- pairwise(treatment, event = r, n = N, studlab = author,
#'   data = dat.woods2010, sm = "OR")
#' net1 <- netmeta(p1, small.values = "good")
#'
#' set.seed(1909) # get reproducible results
#' ran1 <- rankogram(net1, nsim = 100)
#' ran1
#'
#' plot(ran1)
#' plot(ran1, type = "l")
#' plot(ran1, cumulative.rankprob = TRUE)
#' 
#' @method plot rankogram
#' @export

plot.rankogram <- function(x,
                           type = if (cumulative.rankprob) "step" else "bar",
                           pooled = ifelse(x$random, "random", "common"),
                           sort = TRUE, trts,
                           cumulative.rankprob = x$cumulative.rankprob,
                           ylim, ylab,
                           nchar.trts = x$nchar.trts, ...) {  
  
  chkclass(x, c("rankogram"))
  x <- updateversion(x)
  ##
  is_installed_package("gridExtra")
  
  
  type <- setchar(type, c("bar", "line", "step"))
  ##
  pooled <- setchar(pooled, c("common", "random", "fixed"))
  pooled[pooled == "fixed"] <- "common"
  ##
  chklogical(sort)
  ##
  if (is.null(cumulative.rankprob))
    cumulative.rankprob <- FALSE
  chklogical(cumulative.rankprob)
  ##
  missing.ylim <- missing(ylim)
  if (!missing.ylim)
    chknumeric(ylim, min = 0, max = 1, length = 2)
  ##
  if (!missing(ylab))
    chkchar(ylab, length = 1)
  else {
    if (cumulative.rankprob)
      ylab <- "Cumulative\nprobability"
    else
      ylab <- "Probability"
  }
  if (missing(trts))
    trts <- NULL
  ##
  if (is.null(nchar.trts))
    nchar.trts <- 666
  chknumeric(nchar.trts, length = 1)
  
  
  mytheme <-
    theme(axis.line.x = ggplot2::element_line(colour = "black", size = 1,
                                              linetype = "solid"),
          axis.line.y = ggplot2::element_line(colour = "black", size = 1,
                                              linetype = "solid"),
          panel.grid.major = ggplot2::element_line(color = "transparent"),
          panel.grid.minor = ggplot2::element_line(color = "transparent"),
          panel.background = ggplot2::element_rect(fill = "transparent"),
          plot.background = ggplot2::element_rect(
                                       fill = "transparent",
                                       colour = "transparent",
                                       size = 1)
          )
  
  
  if (pooled == "random") {
    if (cumulative.rankprob)
      rankmatrix <- "cumrank.matrix.random"
    else 
      rankmatrix <- "ranking.matrix.random"
    ##
    sucras <- "ranking.random"
  }
  else if (pooled == "common") {
    if (cumulative.rankprob)
      rankmatrix <- "cumrank.matrix.common"
    else 
      rankmatrix <- "ranking.matrix.common"
    ##
    sucras <- "ranking.common"
  }
  ##
  if (is.null(x[[rankmatrix]]))
    stop("No results for ", pooled, " effect",
         if (pooled == "random") "s", " model available. ",
         "Rerun rankogram() with argument '", pooled, " = TRUE'.",
         call. = FALSE)
  
  
  plotranks <- function(treat) {
    df <- data.frame(pos = 1:nrow(x[[rankmatrix]]),
                     ranks = x[[rankmatrix]][treat, ])
    mymaxvalue <- max(x[[rankmatrix]])
    ##
    if (type == "bar")
      p <- ggplot2::ggplot(mapping = aes(df$pos, df$ranks)) +
        ggplot2::geom_col()
    else if (type == "line")
      p <- ggplot2::ggplot(mapping = aes(df$pos, df$ranks)) +
        ggplot2::geom_line()
    else if (type == "step")
      p <- ggplot2::ggplot(mapping = aes(df$pos, df$ranks)) +
        ggplot2::geom_step()
    ##
    p <- p +
      ggplot2::scale_x_continuous(breaks = seq(1, nrow(x[[rankmatrix]]), 1))
    p <- p + ggplot2::labs(x = paste("Rank of", treats(treat, nchar.trts)))
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
  if (sort)
    treatnames <- names(sort(x[[sucras]], decreasing = TRUE))
  else
    treatnames <- names(x[[sucras]])
  ##
  if (!is.null(trts)) {
    trts.c <- setref(trts, treatnames, varname = "trts", length = 0)
    treatnames <- treatnames[treatnames %in% trts.c]
  }
  ##
  rankplots <- do.call(gridExtra::grid.arrange, lapply(treatnames, plotranks))
}
