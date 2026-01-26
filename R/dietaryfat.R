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
#' @seealso \code{\link[meta]{pairwise}}, \code{\link[meta]{metainc}},
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
#' head(dietaryfat)
#' 
#' # Examples: example(netmeta)

NULL
