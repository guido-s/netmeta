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
#' @seealso \code{\link{pairwise}}, \code{\link{metabin}},
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
#' 
#' # Transform data from arm-based format to contrast-based format
#' # Argument 'sm' has to be used for odds ratio as summary measure;
#' # by default the risk ratio is used in the metabin function called
#' # internally.
#' #
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'                event = list(event1, event2, event3),
#'                n = list(n1, n2, n3),
#'                data = smokingcessation,
#'                sm = "OR")
#' p1
#' 
#' # Conduct network meta-analysis
#' #
#' net1 <- netmeta(p1)
#' net1
#' 
#' # Draw network graph
#' #
#' netgraph(net1, points = TRUE, cex.points = 3, cex = 1.25)
#' tname <- c("No intervention", "Self-help",
#'            "Individual counselling", "Group counselling")
#' netgraph(net1, points = TRUE, cex.points = 3, cex = 1.25, labels = tname)


NULL
