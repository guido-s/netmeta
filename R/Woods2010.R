#' Count statistics of survival data
#' 
#' @description
#' Count mortality statistics in randomised controlled trials of
#' treatments for chronic obstructive pulmonary disease (Woods et
#' al. (2010), Table 1).
#' 
#' @name Woods2010
#'
#' @docType data
#'
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{author}}\tab first author / study name \cr
#' \bold{\emph{treatment}}\tab treatment \cr
#' \bold{\emph{r}}\tab number of deaths in treatment arm \cr
#' \bold{\emph{N}}\tab number of patients in treatment arm
#' }
#' 
#' @note
#' The dataset Woods2010 is identical to dataset
#' \code{\link[metadat]{dat.woods2010}} in R package \bold{metadat}.
#' 
#' @seealso \code{\link[metadat]{dat.woods2010}},
#'   \code{\link[meta]{pairwise}}, \code{\link[meta]{metabin}},
#'   \code{\link{netmeta}}
#' 
#' @source
#' Woods BS, Hawkins N, Scott DA (2010):
#' Network meta-analysis on the log-hazard scale, combining count and
#' hazard ratio statistics accounting for multi-arm trials: A
#' tutorial.
#' \emph{BMC Medical Research Methodology},
#' \bold{10}, 54
#' 
#' @keywords datasets
#' 
#' @examples
#' head(dat.woods2010)
#' 
#' # Example using pairwise() and netmeta():
#' # example(dat.woods2010, run.dontrun = TRUE)

NULL
