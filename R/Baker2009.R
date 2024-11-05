#' Network meta-analysis of pharmacologic treatments for chronic
#' obstructive pulmonary disease
#' 
#' @description
#' This dataset comes from a systematic review of randomized
#' controlled trials on pharmacologic treatments for chronic
#' obstructive pulmonary disease (COPD) (Baker et al., 2009).
#'
#' The primary outcome, occurrence of one or more episodes of COPD
#' exacerbation, is binary (yes / no). For this outcome, five drug
#' treatments (fluticasone, budesonide, salmeterol, formoterol,
#' tiotropium) and two combinations (fluticasone + salmeterol,
#' budesonide + formoterol) were compared to placebo. The authors
#' considered the two combinations as separate treatments instead of
#' evaluating the individual components.
#' 
#' @name Baker2009
#' 
#' @docType data
#' 
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{study}}\tab study label \cr
#' \bold{\emph{year}}\tab year of publication \cr
#' \bold{\emph{id}}\tab study ID \cr
#' \bold{\emph{treatment}}\tab treatment \cr
#' \bold{\emph{exac}}\tab one or more episodes of COPD exacerbation
#'   \cr
#' \bold{\emph{total}}\tab number of individuals in treatment arm
#' }
#' 
#' @seealso \code{\link[meta]{pairwise}}, \code{\link[meta]{metabin}},
#'   \code{\link{netmetabin}}
#' 
#' @source
#'
#' Baker WL, Baker EL, Coleman CI (2009):
#' Pharmacologic Treatments for Chronic Obstructive Pulmonary Disease:
#' A Mixed-Treatment Comparison Meta-analysis.
#' \emph{Pharmacotherapy: The Journal of Human Pharmacology and Drug
#' Therapy},
#' \bold{29}, 891--905
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Baker2009)
#' Baker2009
#'
#' \dontrun{
#' # Transform data from long arm-based format to contrast-based
#' # format. Argument 'sm' has to be used for odds ratio as summary
#' # measure; by default the risk ratio is used in the metabin
#' # function called internally.
#' #
#' p1 <- pairwise(treatment, exac, total, studlab = paste(study, year),
#'   data = Baker2009, sm = "OR")
#' 
#' # Conduct network meta-analysis
#' #
#' net1 <- netmeta(p1, ref = "plac")
#' 
#' # Conduct component network meta-analysis
#' #
#' cnet1 <- netcomb(net1)
#' cnet1
#' }


NULL
