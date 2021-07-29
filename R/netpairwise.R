#' Conduct pairwise meta-analyses for all comparisons with direct
#' evidence in a network meta-analysis
#' 
#' @description
#' Conduct pairwise meta-analyses for all comparisons with direct
#' evidence in a network meta-analysis.
#' 
#' @param x An object of class \code{netmeta} or \code{netpairwise}.
#' @param object An object of class \code{netpairwise}.
#' @param separate A logical indicating whether results for
#'   pairwise comparisons should be printed as separate meta-analyses
#'   or as subgroups which is more concise.
#' @param comb.fixed A logical indicating whether a fixed effects
#'   (common effects) network meta-analysis should be conducted.
#' @param comb.random A logical indicating whether a random effects
#'   network meta-analysis should be conducted.
#' @param level The level used to calculate confidence intervals for
#'   individual comparisons.
#' @param level.comb The level used to calculate confidence intervals
#'   for pooled estimates.
#' @param prediction A logical indicating whether prediction intervals
#'   should be printed.
#' @param level.predict The level used to calculate prediction
#'   intervals for a new study.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau}. Either \code{"DL"}, \code{"REML"}, or
#'   \code{"ML"}, can be abbreviated.
#' @param ... Additional arguments (passed on to \code{metagen} or
#'   print functions).
#' 
#' @details
#' Conduct pairwise meta-analyses for all comparisons with direct
#' evidence in a network meta-analysis. In contrast to
#' \code{\link{netmeta}} and \code{\link{netsplit}}, unadjusted
#' standard errors are used in the calculations and the between-study
#' heterogeneity variance is allowed to differ between comparisons.
#' 
#' The R function \code{\link{metagen}} is called internally.
#' 
#' @note
#' This function must not be confused with \code{\link{pairwise}}
#' which can be used as a pre-processing step to convert data from
#' arm-based to contrast-based format by calculating all pairwise
#' comparisons within a study.
#' 
#' @return
#' Either a single \code{\link{metagen}} object with pairwise
#' comparisons as subgroups or a list with \code{\link{metagen}}
#' objects for each direct pairwise comparison.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netsplit}},
#'   \code{\link{pairwise}}
#' 
#' @examples
#' data(Senn2013)
#' oldsets <- settings.meta(digits = 2, digits.tau2 = 2, digits.tau = 2)
#' 
#' # Random effects model
#' #
#' net1 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'                 data = Senn2013, sm = "MD",
#'                 comb.fixed = FALSE)
#' 
#' # Calculate and print consise results for all pairwise
#' # meta-analyses
#' #
#' np1 <- netpairwise(net1)
#' np1
#' print(summary(np1), details.method = FALSE)
#'
#' \dontrun{
#' forest(np1)
#' 
#' # Print detailed information for each pairwise comparison
#' #
#' np2 <- netpairwise(net1, separate = TRUE)
#' forest(np2)
#' }
#'
#' settings.meta(oldsets)
#' 
#' @rdname netpairwise
#' @export
#' @export netpairwise


netpairwise <- function(x,
                        separate = FALSE,
                        comb.fixed = x$comb.fixed,
                        comb.random = x$comb.random,
                        level = x$level,
                        level.comb = x$level.comb,
                        prediction = x$prediction,
                        level.predict = x$level.predict,
                        method.tau = x$method.tau,
                        ...) {
  
  
  ##
  ## Check for netmeta object
  ##
  meta:::chkclass(x, "netmeta")
  
  
  ##
  ## Check other arguments
  ##
  chklevel <- meta:::chklevel
  chklogical <- meta:::chklogical
  ##
  chklogical(separate)
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
  ##
  chklevel(level)
  chklevel(level.comb)
  chklevel(level.predict)
  ##
  method.tau <- meta:::setchar(method.tau, c("DL", "ML", "REML"))
  

  if (!separate) {
    trts.abbr <- treats(x$trts, x$nchar.trts)
    trt1 <- as.character(factor(x$data$.treat1,
                                levels = x$trts, labels = trts.abbr))
    trt2 <- as.character(factor(x$data$.treat2,
                                levels = x$trts, labels = trts.abbr))
    res <-
      metagen(x$data$.TE, x$data$.seTE, studlab = x$data$.studlab,
              sm = x$sm,
              byvar = paste0(trt1, x$sep.trts, trt2),
              bylab = "comparison",
              print.byvar = FALSE,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              level = level,
              level.comb = level.comb,
              prediction = prediction,
              level.predict = level.predict,
              method.tau = method.tau,
              overall = FALSE, overall.hetstat = FALSE,
              test.subgroup = FALSE,
              ...)
    ##
    res$k.study <- x$k
    res$k <- x$m
    res$w.fixed[!is.na(res$w.fixed)] <- NA
    res$w.random[!is.na(res$w.random)] <- NA
    ##
    class(res) <- c(class(res), "netpairwise")
  }
  else {
    comps <- unique(x$data[, c(".treat1", ".treat2")])
    n.comps <- nrow(comps)
    ##
    res <- vector("list", length = n.comps)
    ##
    for (i in seq_len(n.comps)) {
      comp.i <- paste0(comps$.treat1[i], x$sep.trts, comps$.treat2[i])
      res[[i]] <-
        metagen(x$data$.TE, x$data$.seTE, studlab = x$data$.studlab,
                sm = x$sm,
                subset =
                  x$data$.treat1 == comps$.treat1[i] &
                  x$data$.treat2 == comps$.treat2[i],
                complab = comp.i,
                comb.fixed = comb.fixed,
                comb.random = comb.random,
                level = level,
                level.comb = level.comb,
                prediction = prediction,
                level.predict = level.predict,
                method.tau = method.tau,
                ...)
      attr(res[[i]], "comparison") <- comp.i
      ##
      class(res) <- "netpairwise"
    }
  }
  
  res
}





#' @rdname netpairwise
#' @method print netpairwise
#' @export
#' @export print.netpairwise


print.netpairwise <- function(x, ...) {
  
  meta:::chkclass(x, "netpairwise")
  
  n <- 1
  for (i in 1:length(x)) {
    if (n > 1)
      cat("\n*****\n\n")
    print(x[[i]], ...)
    n <- n + 1
  }
  
  invisible(NULL)
}





#' @rdname netpairwise
#' @method summary netpairwise
#' @export
#' @export summary.netpairwise


summary.netpairwise <- function(object, ...) {
  
  meta:::chkclass(object, "netpairwise")
  
  for (i in seq_len(length(object)))
    object[[i]] <- summary(object[[i]])
  ##
  class(object) <- "summary.netpairwise"
  
  object
}





#' @rdname netpairwise
#' @method print summary.netpairwise
#' @export
#' @export print.summary.netpairwise


print.summary.netpairwise <- function(x, ...) {
  
  meta:::chkclass(x, "summary.netpairwise")
  
  n <- 1
  for (i in seq_len(length(x))) {
    if (n > 1)
      cat("\n*****\n\n")
    print(x[[i]], ...)
    n <- n + 1
  }
  
  invisible(NULL)
}





#' @rdname netpairwise
#' @method forest netpairwise
#' @export
#' @export forest.netpairwise


forest.netpairwise <- function(x, ...) {
  
  meta:::chkclass(x, "netpairwise")
  
  for (i in seq_len(length(x)))
    forest(x[[i]],
           smlab = paste0("Comparison:\n", attr(x[[i]], "comparison")),
           ...)
  
  invisible(NULL)
}
