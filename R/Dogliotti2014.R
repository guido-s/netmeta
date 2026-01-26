#' Network meta-analysis of antithrombotic treatments in patients with
#' non-valvular atrial fibrillation
#' 
#' @description
#' This dataset comes from a systematic review aiming to determine
#' the effects of eight antithrombotic treatments in reducing the
#' incidence of major thrombotic events in patients with non-valvular
#' atrial fibrillation (Dogliotti et al., 2014). The review included 20
#' studies (79,808 participants), four of which were three-arm
#' studies. The primary outcome is stroke reduction.
#' 
#' @name Dogliotti2014
#' 
#' @docType data
#' 
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{study}}\tab study label \cr
#' \bold{\emph{id}}\tab study ID \cr
#' \bold{\emph{treatment}}\tab treatment \cr
#' \bold{\emph{stroke}}\tab number of strokes \cr
#' \bold{\emph{total}}\tab number of individuals in treatment arm
#' }
#' 
#' @note
#' The dataset Dogliotti2014 is identical to dataset
#' \code{\link[metadat]{dat.dogliotti2014}} in R package \bold{metadat}.
#' 
#' @seealso \code{\link[metadat]{dat.dogliotti2014}},
#'   \code{\link[meta]{pairwise}}, \code{\link[meta]{metabin}},
#'   \code{\link{netmetabin}}
#' 
#' @source
#' Dogliotti A, Paolasso E, Giugliano RP (2014):
#' Current and new oral antithrombotics in non-valvular atrial
#' fibrillation: a network meta-analysis of 79 808 patients.
#' \emph{Heart},
#' \bold{100}, 396--405
#' 
#' @keywords datasets
#' 
#' @examples
#' head(dat.dogliotti2014)
#' 
#' # Example using pairwise(), netmeta() and netmetabin():
#' # example(dat.dogliotti2014, run.dontrun = TRUE)

NULL
