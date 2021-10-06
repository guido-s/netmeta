#' Network meta-analysis of treatments for depression
#' 
#' @description
#' Network meta-analysis of nine classes of antidepressants including
#' placebo for the primary care setting; partly shown in Linde et
#' al. (2015), supplementary Table 2.
#' 
#' @name Linde2015
#' 
#' @docType data
#' 
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{id}}\tab Study ID \cr
#' \bold{\emph{author}}\tab First author \cr
#' \bold{\emph{year}}\tab Publication year \cr
#' \bold{\emph{treatment1}}\tab First treatment \cr
#' \bold{\emph{treatment2}}\tab Second treatment \cr
#' \bold{\emph{treatment3}}\tab Third treatment \cr
#' \bold{\emph{n1}}\tab Number of patients receiving first treatment
#'   \cr
#' \bold{\emph{resp1}}\tab Number of early responder (treatment 1) \cr
#' \bold{\emph{remi1}}\tab Number of early remissions (treatment 1)
#'   \cr
#' \bold{\emph{loss1}}\tab Number of patients loss to follow-up
#'   (treatment 1) \cr
#' \bold{\emph{loss.ae1}}\tab Number of patients loss to follow-up due
#'   to adverse events (treatment 1) \cr
#' \bold{\emph{ae1}}\tab Number of patients with adverse events
#'   (treatment 1) \cr
#' \bold{\emph{n2}}\tab Number of patients receiving second treatment
#'   \cr
#' \bold{\emph{resp2}}\tab Number of early responder (treatment 2) \cr
#' \bold{\emph{remi2}}\tab Number of early remissions (treatment 2)
#'   \cr
#' \bold{\emph{loss2}}\tab Number of patients loss to follow-up
#'   (treatment 2) \cr
#' \bold{\emph{loss.ae2}}\tab Number of patients loss to follow-up due
#'   to adverse events (treatment 2) \cr
#' \bold{\emph{ae2}}\tab Number of patients with adverse events
#'   (treatment 2) \cr
#' \bold{\emph{n3}}\tab Number of patients receiving third treatment
#'   \cr
#' \bold{\emph{resp3}}\tab Number of early responder (treatment 3) \cr
#' \bold{\emph{remi3}}\tab Number of early remissions (treatment 3)
#'   \cr
#' \bold{\emph{loss3}}\tab Number of patients loss to follow-up
#'   (treatment 3) \cr
#' \bold{\emph{loss.ae3}}\tab Number of patients loss to follow-up due
#'   to adverse events (treatment 3) \cr
#' \bold{\emph{ae3}}\tab Number of patients with adverse events
#'   (treatment 3) 
#' }
#' 
#' @seealso \code{\link{pairwise}}, \code{\link{metabin}},
#'   \code{\link{netmeta}}, \code{\link{netposet}}
#' 
#' @source
#' Linde K, Kriston L, Rücker G, et al. (2015):
#' Efficacy and acceptability of pharmacological treatments for
#' depressive disorders in primary care: Systematic review and network
#' meta-analysis.
#' \emph{Annals of Family Medicine},
#' \bold{13}, 69--79
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Linde2015)
#' 
#' # Transform data from arm-based format to contrast-based format
#' # Outcome: early response
#' p1 <- pairwise(list(treatment1, treatment2, treatment3),
#'                event = list(resp1, resp2, resp3),
#' 	       n = list(n1, n2, n3),
#'                studlab = id, data = Linde2015, sm = "OR")
#' 
#' # Define order of treatments
#' trts <- c("TCA", "SSRI", "SNRI", "NRI",
#'           "Low-dose SARI", "NaSSa", "rMAO-A", "Hypericum",
#'           "Placebo")
#' 
#' # Conduct random effects network meta-analysis
#' net1 <- netmeta(p1, comb.fixed = FALSE,
#'                 reference = "Placebo",
#' 		seq = trts)
#' print(net1, digits = 2)


NULL
