#' Network meta-analysis of primary care depression treatments
#' 
#' @description
#' Network meta-analysis of 22 treatments (including placebo and usual
#' care) for the primary care of depression.
#' 
#' @name Linde2016
#' 
#' @docType data
#' 
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{id}}\tab Study ID \cr
#' \bold{\emph{author}}\tab First author \cr
#' \bold{\emph{year}}\tab Year of publication \cr
#' \bold{\emph{treat1}}\tab First treatment (abbreviated) \cr
#' \bold{\emph{treat2}}\tab Second treatment (abbreviated) \cr
#' \bold{\emph{treat1.long}}\tab First treatment \cr
#' \bold{\emph{treat2.long}}\tab Second treatment \cr
#' \bold{\emph{lnOR}}\tab Response after treatment (log odds ratio) \cr
#' \bold{\emph{selnOR}}\tab Standard error of log odds ratio \cr
#' \bold{\emph{resp1}}\tab Responder (first treatment) \cr
#' \bold{\emph{n1}}\tab Sample size (first treatment) \cr
#' \bold{\emph{resp2}}\tab Responder (second treatment) \cr
#' \bold{\emph{n2}}\tab Sample size (second treatment) \cr
#' }
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netcomb}}
#' 
#' @source
#' Linde K, Rücker G, Schneider A et al. (2016):
#' Questionable assumptions hampered interpretation of a network
#' meta-analysis of primary care depression treatments.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{71}, 86--96
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Linde2016)
#' 
#' # Only consider studies including Face-to-face PST (to reduce
#' # runtime of example)
#' #
#' face <- subset(Linde2016, id %in% c(16, 24, 49, 118))
#' 
#' # Conduct random effects network meta-analysis
#' #
#' net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'   data = face, reference.group = "placebo",
#'   sm = "OR", common = FALSE, nchar.trts = 6)
#' #
#' net1
#' 
#' \dontrun{
#' # Conduct random effects network meta-analysis
#' #
#' net2 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'   data = Linde2016, reference.group = "placebo",
#'   sm = "OR", common = FALSE, nchar.trts = 6)
#' #
#' net2
#' }


NULL
