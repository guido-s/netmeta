#' Network meta-analysis of treatments for Parkinson's disease
#' 
#' @description
#' Network meta-analysis comparing the effects of a number of
#' treatments for Parkinson's disease.
#' 
#' The data are the mean lost work-time reduction in patients given
#' dopamine agonists as adjunct therapy in Parkinson’s disease
#' (Franchini et al. 2012). The data are given as sample size, mean
#' and standard deviation in each trial arm. Treatments are placebo
#' and four active drugs. These data are used as an example in the
#' supplemental material of Dias et al. (2013) where placebo is coded
#' as 1 and the four active drugs as 2 to 5.
#' 
#' @name Franchini2012
#' @aliases parkinson
#' 
#' @docType data
#' 
#' @format A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{Study}}\tab study label \cr
#' \bold{\emph{Treatment1}}\tab treatment 1 \cr
#' \bold{\emph{y1}}\tab treatment effect arm 1 \cr
#' \bold{\emph{sd1}}\tab Standard deviation arm 1 \cr
#' \bold{\emph{n1}}\tab Sample size arm 1 \cr
#' \bold{\emph{Treatment2}}\tab treatment 2 \cr
#' \bold{\emph{y2}}\tab treatment effect arm 2 \cr
#' \bold{\emph{sd2}}\tab Standard deviation arm 2 \cr
#' \bold{\emph{n2}}\tab Sample size arm 2 \cr
#' \bold{\emph{Treatment3}}\tab treatment 3 \cr
#' \bold{\emph{y3}}\tab treatment effect arm 3 \cr
#' \bold{\emph{sd3}}\tab Standard deviation arm 3 \cr
#' \bold{\emph{n3}}\tab Sample size arm 3
#' }
#' 
#' @note
#' The dataset Franchini2012 is identical to dataset
#' \code{\link[metadat]{dat.franchini2012}} in R package \bold{metadat}.
#' 
#' @seealso \code{\link[metadat]{dat.franchini2012}},
#'   \code{\link[meta]{pairwise}}, \code{\link[meta]{metacont}},
#'   \code{\link{netmeta}}, \code{\link{netgraph.netmeta}}
#' 
#' @source
#' Dias S, Sutton AJ, Ades AE and Welton NJ (2013):
#' Evidence synthesis for decision making 2: A generalized linear
#' modeling framework for pairwise and network meta-analysis of
#' randomized controlled trials.
#' \emph{Medical Decision Making},
#' \bold{33}, 607--17
#'
#' Franchini AJ, Dias S, Ades AE, Jansen JP, Welton NJ (2012):
#' Accounting for correlation in network meta-analysis with multi-arm
#' trials.
#' \emph{Research Synthesis Methods},
#' \bold{3}, 142--60
#' 
#' @keywords datasets
#' 
#' @examples
#' head(dat.franchini2012)
#' 
#' # Example using pairwise() and netmeta():
#' # example(dat.franchini2012, run.dontrun = TRUE)

NULL
