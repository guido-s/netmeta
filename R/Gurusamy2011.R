#' Network meta-analysis on blood loss during liver transplantation
#' 
#' @description
#' Network meta-analysis comparing the effects of a number of
#' interventions for decreasing blood loss and blood transfusion
#' requirements during liver transplantation.
#' 
#' @name Gurusamy2011
#' 
#' @docType data
#' 
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{study}}\tab study information (first author, year) \cr
#' \bold{\emph{treatment}}\tab treatment \cr
#' \bold{\emph{death}}\tab mortality at 60 days post-transplantation
#'   \cr
#' \bold{\emph{n}}\tab number of individuals in treatment arm
#' } 
#' 
#' @note
#' The dataset Gurusamy2011 is identical to dataset
#' \code{\link[metadat]{dat.gurusamy2011}} in R package \bold{metadat}.
#' 
#' @seealso \code{\link[metadat]{dat.gurusamy2011}},
#'   \code{\link[meta]{pairwise}}, \code{\link[meta]{metabin}},
#'   \code{\link{netmetabin}}
#' 
#' @source
#' Gurusamy KS, Pissanou T, Pikhart H, Vaughan J, Burroughs AK,
#' Davidson BR (2011):
#' Methods to decrease blood loss and transfusion requirements for
#' liver transplantation.
#' \emph{Cochrane Database of Systematic Reviews},
#' CD009052
#' 
#' @keywords datasets
#' 
#' @examples
#' head(dat.gurusamy2011)
#' 
#' # Example using pairwise() and netmetabin():
#' # example(dat.gurusamy2011, run.dontrun = TRUE)

NULL
