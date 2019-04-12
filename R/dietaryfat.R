#' Network meta-analysis of dietary fat
#' 
#' @description
#' Network meta-analysis comparing the effects of two diets to control
#' on mortality.
#' 
#' The data are rates, given as the number of deaths and person-years.
#' These data are used as an example in the supplemental material of
#' Dias et al. (2013).
#' 
#' @name dietaryfat
#' 
#' @docType data
#' 
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{treat1}}\tab treatment 1 \cr
#' \bold{\emph{treat2}}\tab treatment 2 \cr
#' \bold{\emph{treat3}}\tab treatment 3 \cr
#' \bold{\emph{years1}}\tab person years arm 1 \cr
#' \bold{\emph{years2}}\tab person years arm 2 \cr
#' \bold{\emph{years3}}\tab person years arm 3 \cr
#' \bold{\emph{d1}}\tab events (deaths) arm 1 \cr
#' \bold{\emph{d2}}\tab events (deaths) arm 2 \cr
#' \bold{\emph{d3}}\tab events (deaths) arm 3 \cr
#' \bold{\emph{ID}}\tab study ID
#' }
#' 
#' @seealso \code{\link{pairwise}}, \code{\link{metainc}},
#' \code{\link{netmeta}}, \code{\link{netgraph.netmeta}}
#' 
#' @source
#' Dias S, Sutton AJ, Ades AE and Welton NJ (2013):
#' Evidence synthesis for decision making 2: A generalized linear
#' modeling framework for pairwise and network meta-analysis of
#' randomized controlled trials.
#' \emph{Medical Decision Making},
#' \bold{33}, 607--17
#' 
#' @keywords datasets
#' 
#' @examples
#' data(dietaryfat)
#' 
#' # Transform data from arm-based format to contrast-based format
#' # Using incidence rate ratios (sm = "IRR") as effect measure.
#' # Note, the argument 'sm' is not necessary as this is the default
#' # in R function metainc() called internally
#' #
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'                list(d1, d2, d3),
#'                time = list(years1, years2, years3),
#'                studlab = ID,
#'                data = dietaryfat, sm = "IRR")
#' p1
#' 
#' # Conduct network meta-analysis
#' #
#' net1 <- netmeta(p1)
#' summary(net1)
#' 
#' # Conduct network meta-analysis using incidence rate differences
#' # (sm = "IRD")
#' #
#' p2 <- pairwise(list(treat1, treat2, treat3),
#'                list(d1, d2, d3),
#'                time = list(years1, years2, years3),
#'                studlab = ID,
#'                data = dietaryfat, sm = "IRD")
#' net2 <- netmeta(p2)
#' summary(net2)
#' 
#' # Draw network graph
#' #
#' netgraph(net1, points = TRUE, cex.points = 3, cex = 1.25)
#' tname <- c("Control","Diet", "Diet 2")
#' netgraph(net1, points = TRUE, cex.points = 3, cex = 1.25, labels = tname)


"dietaryfat"
