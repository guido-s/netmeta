#' Network meta-analysis of interventions for smoking cessation
#' 
#' @description
#' Network meta-analysis comparing the effects of a number of
#' interventions for smoking cessation.
#' 
#' These data are used as an example in Dias et al. (2013), page 651.
#' 
#' @name smokingcessation
#' 
#' @docType data
#' 
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{event1}}\tab number of individuals with successful
#'   smoking cessation in arm 1 \cr
#' \bold{\emph{n1}}\tab number of individuals in arm 1 \cr
#' \bold{\emph{event2}}\tab number of individuals with successful
#'   smoking cessation in arm 2 \cr
#' \bold{\emph{n2}}\tab number of individuals in arm 2 \cr
#' \bold{\emph{event3}}\tab number of individuals with successful
#'   smoking cessation in arm 3 \cr
#' \bold{\emph{n3}}\tab number of individuals in arm 3 \cr
#' \bold{\emph{treat1}}\tab treatment 1 \cr
#' \bold{\emph{treat2}}\tab treatment 2 \cr \bold{\emph{treat3}}\tab
#'   treatment 3
#' }
#' 
#' @seealso \code{\link[meta]{pairwise}}, \code{\link[meta]{metabin}},
#'   \code{\link{netmeta}}, \code{\link{netgraph.netmeta}}
#' 
#' @source
#' Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G and Ades AE (2013):
#' Evidence Synthesis for Decision Making 4: Inconsistency in networks
#' of evidence based on randomized controlled trials.
#' \emph{Medical Decision Making},
#' \bold{33}, 641--56
#' 
#' @keywords datasets
#' 
#' @examples
#' data(smokingcessation)
#' head(smokingcessation)
#' 
#' # Examples: example(netmeta)

NULL
