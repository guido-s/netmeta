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
#' @seealso \code{\link[meta]{pairwise}}, \code{\link[meta]{metabin}},
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
#' data(Gurusamy2011)
#' 
#' # Only consider three studies (to reduce runtime of example)
#' #
#' studies <- c("Findlay 2001", "Garcia-Huete 1997", "Dalmau 2000")
#' three <- subset(Gurusamy2011, study %in% studies)
#' 
#' # Transform data from long arm-based format to contrast-based
#' # format. Argument 'sm' has to be used for odds ratio as summary
#' # measure; by default the risk ratio is used in the metabin
#' # function called internally.
#' #
#' p1 <- pairwise(treatment, death, n, studlab = study,
#'   data = three, sm = "OR")
#' 
#' # Conduct Mantel-Haenszel network meta-analysis
#' #
#' netmetabin(p1, ref = "cont")
#' 
#' \dontrun{
#' p2 <- pairwise(treatment, death, n, studlab = study,
#'   data = Gurusamy2011, sm = "OR")
#' 
#' # Conduct Mantel-Haenszel network meta-analysis
#' netmetabin(p2, ref = "cont")
#' }


NULL
