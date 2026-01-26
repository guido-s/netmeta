#' Network meta-analysis for chronic obstructive pulmonary disease
#' 
#' @description
#' Network meta-analysis comparing inhaled medications in patients
#' with chronic obstructive pulmonary disease.
#' 
#' @name Dong2013
#' 
#' @docType data
#' 
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{id}}\tab study ID \cr
#' \bold{\emph{treatment}}\tab treatment \cr
#' \bold{\emph{death}}\tab mortality \cr
#' \bold{\emph{randomized}}\tab number of individuals in treatment arm
#' }
#' 
#' @note
#' The dataset Dong2013 is identical to dataset
#' \code{\link[metadat]{dat.dong2013}} in R package \bold{metadat}.
#' 
#' @seealso \code{\link[metadat]{dat.dong2013}},
#'   \code{\link[meta]{pairwise}}, \code{\link[meta]{metabin}},
#'   \code{\link{netmetabin}}
#' 
#' @source
#' Dong Y-H, Lin H-H, Shau W-Y, Wu Y-C, Chang C-H, Lai M-S (2013):
#' Comparative safety of inhaled medications in patients with chronic
#' obstructive pulmonary disease: systematic review and mixed
#' treatment comparison meta-analysis of randomised controlled trials.
#' \emph{Thorax},
#' \bold{68}, 48--56
#' 
#' @keywords datasets
#' 
#' @examples
#' head(dat.dong2013)
#' 
#' # Example using pairwise() and netmetabin():
#' # example(dat.dong2013, run.dontrun = TRUE)

NULL
