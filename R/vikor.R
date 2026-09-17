#' Rank treatments across all outcomes using the VišeKriterijumska Optimizacija
#' I Kompromisno Rešenje (VIKOR) multi-criteria decision analysis method
#' 
#' @description
#' This function employs the VišeKriterijumska Optimizacija I Kompromisno
#' Rešenje (VIKOR) method to analyze all outcome-specific ranking lists.
#' It provides both an amalgamated ranking list and guidance on
#' which treatments correspond to the best compromise solutions.
#' 
#' @param x An object of class \code{\link{netposet}} or a matrix.
#' @param pooled A character string indicating whether the VIKOR method should
#'   be based on the common (\code{"common"}) or random effects
#'   model (\code{"random"}). Can be abbreviated.
#' @param weights Outcome weights. The weights should always sum to 1. If not
#'   then they are standardized. If NULL, the function will assume equal outcome
#'   weights. 
#' @param v A scalar from 0 to 1 interpreted as the weight of the decision
#'   making process. Following guidance from the multi-criteria decision
#'   analysis field it is set to 0.5.
#' @param method A character string specifying the ranking metric. Either
#'   \code{"P-score"}, \code{"SUCRA"}, \code{"best"}, or
#'   \code{"ranking probabilities"}; can be abbreviated.
#' @param digits A numeric specifying the number of digits to print the
#'   ranking matrix Q.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' This function takes a single mandatory argument, which is either an object
#' of class \code{\link{netposet}}, a matrix, or a data frame. It then uses the
#' multi-criteria decision analysis method VišeKriterijumska Optimizacija I
#' Kompromisno Rešenje (VIKOR) to produce an amalgamated ranking list across
#' all outcomes (Opricovic & Tzeng, 2004).
#' 
#' The VIKOR approach is only applied when the \code{method} argument is
#' equal to \code{"P-score"}, \code{"SUCRA"}, \code{"best"}, or
#' \code{"ranking probabilities"} in \code{\link{netposet}}.
#' 
#' The final ranking list is calculated based on treatments common across all
#' outcomes. Treatments not present across all outcomes are excluded internally.
#' 
#' Using the argument 'weights' the users can specify the weight that each
#' outcome should have in the decision making process. For each outcome this
#' argument should have a value from 0 to 1 while the sum of all outcome
#' weights should be 1. If the sum of all weights is not 1, they are
#' internally normalized to sum to 1 and a warning with the normalized weight
#' values is printed. Finally, if NULL then equal weights are assumed across
#' all outcomes.
#'
#' The argument 'v' specifies the weight of the decision making process.
#' The VIKOR method is a compromise programming approach that aims to balance
#' between each treatments overall and worst performance across all outcomes.
#' The balance between these two criteria is achieved using the parameter 'v'
#' which takes values from 0 to 1. Values close to 1 will give more weight to
#' the treatment's overall performance while values close to 0 will give more
#' weight to penalize the treatment's worst performance. The most common
#' choice of 'v' is typically 0.5 (default also here), thereby allowing for a
#' balanced decision making between treatment's overall and worst performance.
#' 
#' @return
#' The function returns a 'vikor' object. This consists of three ranking lists
#' which are the following:
#' \itemize{
#' \item A ranking list Q referring to the ranking when balancing both each
#'   treatment's overall and worst performance. This is the main ranking list
#'   of the method. 
#' \item A ranking list S referring to the ranking in terms of each treatment's
#'   overall performance.
#' \item A ranking list R referring to the ranking in terms of penalising each
#'   treatment's worst performance.
#' }
#' In addition to the ranking lists, the function also evaluates the necessary
#' conditions defined by the VIKOR method and returns a message indicating the
#' set of compromise solutions.
#'
#' @references
#' Opricovic S, Tzeng GH (2004):
#' Compromise solution by MCDM methods: A comparative analysis of VIKOR and
#' TOPSIS.
#' \emph{European Journal of Operational Research},
#' \bold{156}, 445--55
#' 
#' @examples
#' \donttest{
#' # Define order of treatments in depression data set Linde2015
#' #
#' trts <- c("TCA", "SSRI", "SNRI", "NRI",
#'   "Low-dose SARI", "NaSSa", "rMAO-A", "Hypericum", "Placebo")
#' 
#' # Outcome labels
#' #
#' outcomes <- c("Early response", "Early remission")
#' 
#' # (1) Early response
#' #
#' pw1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#' #
#' nma1 <- netmeta(pw1, common = FALSE,
#'   seq = trts, ref = "Placebo", small.values = "undesirable")
#' 
#' # (2) Early remission
#' #
#' pw2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#' #
#' nma2 <- netmeta(pw2, common = FALSE,
#'   seq = trts, ref = "Placebo", small.values = "undesirable")
#' 
#' # Partial order of treatment rankings (two outcomes)
#' #
#' po12 <- netposet(netrank(nma1), netrank(nma2), outcomes = outcomes)
#' 
#' # Get the best compromise solution across the efficacy outcomes
#' vikor(po12)
#' 
#' # Use larger weight for response than remission
#' vikor(po12, weights = c(0.6, 0.4))
#'
#' 
#' # Consider five outcomes
#' #
#' # Outcome labels
#' #
#' outcomes <- c("Early response", "Early remission",
#'   "Lost to follow-up", "Lost to follow-up due to AEs",
#'    "Adverse events (AEs)")
#' 
#' # (3) Loss to follow-up
#' #
#' pw3 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss1, loss2, loss3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#' #
#' nma3 <- netmeta(pw3, common = FALSE,
#'   seq = trts, ref = "Placebo", small.values = "desirable")
#' 
#' # (4) Loss to follow-up due to adverse events
#' #
#' pw4 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss.ae1, loss.ae2, loss.ae3), n = list(n1, n2, n3),
#'   studlab = id, data = subset(Linde2015, id != 55), sm = "OR")
#' #
#' nma4 <- netmeta(pw4, common = FALSE,
#'   seq = trts, ref = "Placebo", small.values = "desirable")
#' 
#' # (5) Adverse events
#' #
#' pw5 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(ae1, ae2, ae3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#' #
#' nma5 <- netmeta(pw5, common = FALSE,
#'   seq = trts, ref = "Placebo", small.values = "desirable")
#' 
#' # Partial order of treatment rankings (based on netrank() objects)
#' #
#' po12345 <- netposet(netrank(nma1), netrank(nma2),
#'   netrank(nma3), netrank(nma4), netrank(nma5), outcomes = outcomes)
#' 
#' # Get the best compromise solution across all outcomes
#' vikor(po12345)
#' 
#' # Use larger weight for efficacy than safety outcomes
#' vikor(po12345, weights = c(0.35, 0.35, 0.1, 0.1, 0.1))
#' 
#' # Example using ranking matrix with P-scores
#' #
#' # Ribassin-Majed L, Marguet S, Lee A.W., et al. (2017):
#' # What is the best treatment of locally advanced nasopharyngeal
#' # carcinoma? An individual patient data network meta-analysis.
#' # Journal of Clinical Oncology, 35, 498-505
#' #
#' # P-scores (from Table 1)
#' #
#' pscore.os  <- c(15, 33, 63, 70, 96, 28, 45) / 100
#' pscore.pfs <- c( 4, 46, 79, 52, 94, 36, 39) / 100
#' pscore.lc  <- c( 9, 27, 47, 37, 82, 58, 90) / 100
#' pscore.dc  <- c(16, 76, 95, 48, 72, 32, 10) / 100
#' #
#' pscore.matrix <- data.frame(pscore.os, pscore.pfs, pscore.lc, pscore.dc)
#' rownames(pscore.matrix) <-
#'   c("RT", "IC-RT", "IC-CRT", "CRT", "CRT-AC", "RT-AC", "IC-RT-AC")
#' colnames(pscore.matrix) <- c("OS", "PFS", "LC", "DC")
#' pscore.matrix
#' #
#' po <- netposet(pscore.matrix)
#' vikor(po)
#' # same result
#' vikor(pscore.matrix)
#' }
#'
#' @rdname vikor
#' @method vikor netposet
#' @export

vikor.netposet <- function(x,
                           pooled = ifelse(x$random, "random", "common"),
                           weights = NULL, v = 0.5, ...) {
  
  chkclass(x, "netposet")
  x  <- updateversion(x)
  #
  pooled <- setchar(pooled, c("common", "random", "fixed"))
  pooled[pooled == "fixed"] <- "common"
  #
  chknumeric(v, min = 0, max = 1, length = 1)
  
  if (x$method %in% c("P-score", "SUCRA", "best", "ranking probabilities")) {
    if (pooled == "common")
      res <- vikor_internal(x$P.common, weights = weights, v = v)
    else
      res <- vikor_internal(x$P.random, weights = weights, v = v)
  }
  else {
    stop("VIKOR method is only available for P-scores, SUCRAs, ",
         "probabilities of being best, and ranking probabilities.",
         call. = FALSE)
  }
  #
  attr(res, "ranking.method") <- x$method
  #
  res
}


#' @rdname vikor
#' @method vikor matrix
#' @export

vikor.matrix <- function(x, weights = NULL, v = 0.5, method = "SUCRA", ...) {
  
  chkclass(x, "matrix")
  #
  chknumeric(v, min = 0, max = 1, length = 1)
  #
  method <-
    setchar(method, c("P-score", "SUCRA", "best", "ranking probabilities"))
  #
  res <- vikor_internal(x, weights = weights, v = v)
  #
  attr(res, "ranking.method") <- method
  #
  res
}


#' @rdname vikor
#' @method vikor data.frame
#' @export

vikor.data.frame <- function(x, weights = NULL, v = 0.5,
                             method = "SUCRA", ...) {
  
  chkclass(x, "data.frame")
  #
  chknumeric(v, min = 0, max = 1, length = 1)
  #
  method <-
    setchar(method, c("P-score", "SUCRA", "best", "ranking probabilities"))
  #
  res <- vikor_internal(as.matrix(x), weights = weights, v = v)
  #
  attr(res, "ranking.method") <- method
  #
  res
}


#' @rdname vikor
#' @export vikor

vikor <- function(x, ...)
  UseMethod("vikor")


#' @rdname vikor
#' @method print vikor
#' @export

print.vikor <- function(x, digits = 4, ...) {
  
  chkclass(x, "vikor")
  #
  chknumeric(digits, min = 0, length = 1)
    
  Q <- x %>% select(Q)
  S <- x %>% select(S)
  R <- x %>% select(R)
  #
  trts <- row.names(Q)
  #
  DQ <- 1 / (length(trts) - 1)
  
  cond1 <- Q$Q[2] - Q$Q[1] >= DQ
  #
  cond2_1 <- isTRUE(row.names(Q)[1] == row.names(S)[1])
  cond2_2 <- isTRUE(row.names(Q)[1] == row.names(R)[1])
  #
  cond2 <- isTRUE(cond2_1 & cond2_2)
  #
  if (cond1 & cond2) {
    solution <- row.names(Q)[1]
    #
    txt <- paste("Compromise treatment across all outcomes:", solution)
  }
  else if ((cond1) & (!cond2)) {
    solution <- paste(row.names(Q)[1:2], collapse = ", ")
    #
    txt <- paste("Compromise set of treatments across all outcomes:",
                 solution)
  }
  else if (!cond1) {
    compr <- Q$Q - Q$Q[1] < DQ
    #
    E <- which(compr)
    #
    solution <- paste(row.names(Q)[E], collapse = ", ")
    #
    txt <- paste("Compromise set of treatments across all outcomes:",
                 solution)
  }
  else if (!cond1 & !cond2)
    txt <- paste("No compromise solution was identified. Please consider",
                 "different outcome weights.")
  
  res_mat <- cbind(Q, S, R)
  
  if (attr(x, "ranking.method") %in%
      c("P-score", "SUCRA", "best", "ranking probabilities"))
    cat("VIKOR results\n\n")
  #
  prmatrix(round(res_mat, digits = digits), quote = FALSE, right = TRUE)
  #
  cat(paste0("\n", txt, "\n"))
  #
  cat(paste("Threshold for acceptable advantage:",
            round(1 / (nrow(res_mat) - 1), 3), "\n"))
  #
  drop_trts <- attr(x, "dropped_treatments")
  if (length(drop_trts) > 0) {
    cat(paste0("\nThe following treatment",
               if (length(drop_trts) > 1) "s are" else " is",
               " not considered in the VIKOR method due to ",
               "missing information: ",
               paste(drop_trts, collapse = ", "),
               "\n"))
  }
  #
  invisible(NULL)
}
