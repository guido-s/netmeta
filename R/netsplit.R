#' Split direct and indirect evidence in network meta-analysis
#' 
#' @description
#' Methods to split network estimates into the contribution of direct
#' and indirect evidence and to test for local inconsistency in
#' network meta-analysis.
#' 
#' @aliases netsplit print.netsplit
#' 
#' @param x An object of class \code{netmeta} or \code{netsplit}.
#' @param method A character string indicating which method to split
#'   direct and indirect evidence is to be used. Either
#'   \code{"Back-calculation"} or \code{"SIDDE"}, can be abbreviated.
#'   See Details.
#' @param upper A logical indicating whether treatment comparisons
#'   should be selected from the lower or upper triangle of the
#'   treatment effect matrices (see list elements \code{TE.common} and
#'   \code{TE.random} in the \code{netmeta} object). Ignored if
#'   argument \code{order} is provided.
#' @param reference.group Reference treatment. Ignored if argument
#'   \code{order} is provided.
#' @param baseline.reference A logical indicating whether results
#'   should be expressed as comparisons of other treatments versus the
#'   reference treatment or vice versa. This argument is only
#'   considered if \code{reference.group} is not equal to \code{""}
#'   and argument\code{order} is not provided.
#' @param order A optional character or numerical vector specifying
#'   the order of treatments in comparisons.
#' @param sep.trts A character string used in comparison names as
#'   separator between treatment labels, e.g., " vs ".
#' @param quote.trts A character used to print around treatment
#'   labels.
#' @param tol.direct A numeric defining the maximum deviation of the
#'   direct evidence proportion from 0 or 1 to classify a comparison
#'   as providing only indirect or direct evidence, respectively.
#' @param common A logical indicating whether results for the common
#'   effects network meta-analysis should be printed.
#' @param random A logical indicating whether results for the random
#'   effects network meta-analysis should be printed.
#' @param show A character string indicating which comparisons should
#'   be printed (see Details).
#' @param overall A logical indicating whether estimates from network
#'   meta-analyis should be printed in addition to direct and indirect
#'   estimates.
#' @param ci A logical indicating whether confidence intervals should
#'   be printed in addition to treatment estimates.
#' @param test A logical indicating whether results of a test
#'   comparing direct and indirect estimates should be printed.
#' @param only.reference A logical indicating whether only comparisons
#'   with the reference group should be printed.
#' @param sortvar An optional vector used to sort comparisons (must be
#'   of same length as the total number of comparisons).
#' @param subset An optional logical vector specifying a subset of
#'   comparisons to print (must be of same length as the total number of
#'   comparisons) .
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z-value
#'   of test of agreement between direct and indirect evidence, see
#'   \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of test of agreement between direct and indirect evidence, see
#'   \code{print.default}.
#' @param digits.prop Minimal number of significant digits for direct
#'   evidence proportions, see \code{print.default}.
#' @param text.NA A character string specifying text printed for
#'   missing values.
#' @param backtransf A logical indicating whether printed results
#'   should be back transformed. For example, if \code{backtransf =
#'   TRUE}, results for \code{sm = "OR"} are printed as odds ratios
#'   rather than log odds ratios.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param indent A logical indicating whether items in the legend
#'   should be indented.
#' @param warn A logical indicating whether warnings should be
#'   printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param verbose A logical indicating whether progress information
#'   should be printed.
#' @param \dots Additional arguments.
#' 
#' @details
#' A comparison of direct and indirect treatment estimates can serve
#' as check for consistency of network meta-analysis (Dias et al.,
#' 2010).
#' 
#' This function provides two methods to derive indirect estimates:
#' \itemize{
#' \item Separate Indirect from Direct Evidence (SIDE) using a
#'   back-calculation method (\code{method = "Back-calculation"})
#'   based on the \emph{direct evidence proportion} to calculate the
#'   indirect evidence (König et al., 2013);
#' \item Separate Indirect from Direct Design Evidence (SIDDE) as
#'   described in Efthimiou et al. (2019).
#' }
#' 
#' Note, for the back-calculation method, indirect treatment estimates
#' are already calculated in \code{\link{netmeta}} and this function
#' combines and prints these estimates in a user-friendly
#' way. Furthermore, this method is not available for the
#' Mantel-Haenszel and non-central hypergeometric distribution
#' approach implemented in \code{\link{netmetabin}}.
#' 
#' For the random-effects model, the direct treatment estimates are
#' based on the common between-study variance \eqn{\tau^2} from the
#' network meta-analysis, i.e. the square of list element
#' \code{x$tau}.
#' 
#' Argument \code{show} determines which comparisons are printed:
#' \tabular{ll}{
#' \dQuote{all} \tab All comparisons \cr
#' \dQuote{both} \tab Only comparisons contributing both direct and
#'   indirect evidence \cr
#' \dQuote{with.direct} \tab Comparisons providing direct evidence \cr
#' \dQuote{direct.only} \tab Comparisons providing only direct
#'   evidence \cr
#' \dQuote{indirect.only} \tab Comparisons providing only indirect
#'   evidence
#' }
#'
#' The node-splitting method and SIDDE can be compute-intensive in
#' large networks. Crude information on the computation progress is
#' printed if argument \code{verbose = TRUE}. In addition, computation
#' times are printed if R package \bold{tictoc} is installed.
#'
#' @return
#' An object of class \code{netsplit} with corresponding \code{print}
#' and \code{forest} functions. The object is a list containing the
#' following components:
#' \item{common, random}{As defined above.}
#' \item{comparison}{A vector with treatment comparisons.}
#' \item{prop.common, prop.random}{A vector with direct evidence
#'   proportions (common / random effects model).}
#' \item{common, random}{Results of network meta-analysis (common /
#'   random effects model), i.e., data frame with columns comparison,
#'   TE, seTE, lower, upper, z, and p.}
#' \item{direct.common, direct.random}{Network meta-analysis results
#'   based on direct evidence (common / random effects model), i.e.,
#'   data frame with columns comparison, TE, seTE, lower, upper, z,
#'   and p.}
#' \item{indirect.common, indirect.random}{Network meta-analysis
#'   results based on indirect evidence (common / random effects
#'   model), i.e., data frame with columns comparison, TE, seTE,
#'   lower, upper, z, and p.}
#' \item{compare.common, compare.random}{Comparison of direct and
#'   indirect evidence in network meta-analysis (common / random
#'   effects model), i.e., data frame with columns comparison, TE,
#'   seTE, lower, upper, z, and p.}
#' \item{sm}{A character string indicating underlying summary measure}
#' \item{level.ma}{The level used to calculate confidence intervals
#'   for pooled estimates.}
#' \item{tictoc}{Computation times for node-splitting method or SIDDE
#'   (if R package \bold{tictoc} is installed).}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}, Gerta
#'   Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}, Orestis Efthimiou
#'   \email{oremiou@@gmail.com}
#' 
#' @seealso \code{\link{forest.netsplit}}, \code{\link{netmeta}},
#'   \code{\link{netmetabin}}, \code{\link{netmeasures}}
#' 
#' @references
#' Dias S, Welton NJ, Caldwell DM, Ades AE (2010):
#' Checking consistency in mixed treatment comparison meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{29}, 932--44
#' 
#' Efthimiou O, Rücker G, Schwarzer G, Higgins J, Egger M, Salanti G
#' (2019):
#' A Mantel-Haenszel model for network meta-analysis of rare events.
#' \emph{Statistics in Medicine},
#' \bold{38}, 2992--3012
#' 
#' König J, Krahn U, Binder H (2013):
#' Visualizing the flow of evidence in network meta-analysis and
#' characterizing mixed treatment comparisons.
#' \emph{Statistics in Medicine},
#' \bold{32}, 5414--29
#' 
#' Puhan MA, Schünemann HJ, Murad MH, et al. (2014):
#' A GRADE working group approach for rating the quality of treatment
#' effect estimates from network meta-analysis.
#' \emph{British Medical Journal},
#' \bold{349}, g5630
#' 
#' @examples
#' data(Woods2010)
#' #
#' p1 <- pairwise(treatment, event = r, n = N,
#'   studlab = author, data = Woods2010, sm = "OR")
#' #
#' net1 <- netmeta(p1)
#' #
#' print(netsplit(net1), digits = 2)
#' 
#' \dontrun{
#' print(netsplit(net1), digits = 2,
#'   backtransf = FALSE, common = FALSE)
#'
#' # Sort by increasing number of studies in direct comparisons
#' print(netsplit(net1), digits = 2, sortvar = k)
#' # Sort by decreasing number of studies in direct comparisons
#' print(netsplit(net1), digits = 2, sortvar = -k)
#' 
#' # Sort by increasing evidence proportion under common effects model
#' print(netsplit(net1), digits = 2, sortvar = prop.common)
#' # Sort by decreasing evidence proportion under common effects model
#' print(netsplit(net1), digits = 2, sortvar = -prop.common)
#' 
#' # Sort by decreasing evidence proportion under common effects model
#' # and number of studies
#' print(netsplit(net1), digits = 2, sortvar = cbind(-prop.common, -k))
#' 
#' data(Senn2013)
#' #
#' net2 <- netmeta(TE, seTE, treat1.long, treat2.long,
#'   studlab, data = Senn2013)
#' #
#' print(netsplit(net2), digits = 2)
#' # Layout of Puhan et al. (2014), Table 1
#' print(netsplit(net2), digits = 2, ci = TRUE, test = FALSE)
#' 
#' data(Dong2013)
#' p3 <- pairwise(treatment, death, randomized, studlab = id,
#'   data = Dong2013, sm = "OR")
#' net3 <- netmetabin(p3)
#' netsplit(net3)
#' }
#' 
#' @rdname netsplit
#' @export netsplit


netsplit <- function(x, method,
                     upper = TRUE,
                     reference.group = x$reference.group,
                     baseline.reference = x$baseline.reference,
                     order = NULL,
                     sep.trts = x$sep.trts, quote.trts = "",
                     tol.direct = 0.0005,
                     common = x$common,
                     random = x$random,
                     backtransf = x$backtransf,
                     warn = FALSE, warn.deprecated = gs("warn.deprecated"),
                     verbose = FALSE,
                     ...) {
  
  ##
  ##
  ## (1) Check for netmeta object and upgrade object
  ##
  ##
  chkclass(x, "netmeta")
  x <- updateversion(x)
  ##
  is.bin <- inherits(x, "netmetabin")
  ##
  is.tictoc <- is.installed.package("tictoc", stop = FALSE)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  if (!missing(method))
    method <- setchar(method, c("Back-calculation", "SIDDE"))
  else {
    if (is.bin)
      method <- "SIDDE"
    else
      method <- "Back-calculation"
  }
  ##
  chklogical(upper)
  chklogical(baseline.reference)
  ##
  if (!is.null(order)) {
    order <- setseq(order, x$trts)
    baseline.reference <- FALSE
    reference.group <- ""
  }
  ##
  chkchar(sep.trts)
  chkchar(quote.trts)
  chknumeric(tol.direct, min = 0, length = 1)
  if (!is.null(backtransf))
    chklogical(backtransf)
  chklogical(warn)
  chklogical(verbose)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args <- list(...)
  chklogical(warn.deprecated)
  ##
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  ##
  x$common <- common
  x$random <- random
  
  
  ##
  ##
  ## (3) Create dat.trts
  ##
  ##
  dat.trts <- comptrts(x, upper, reference.group, baseline.reference,
                       order, sep.trts, quote.trts)
  
  
  ##
  ##
  ## (4) Change order of prop.direct.common and prop.direct.random
  ##
  ##
  if (!(is.bin & method == "SIDDE")) {
    prop.common <- sortprop(x$prop.direct.common, dat.trts, x$sep.trts)
    prop.random <- sortprop(x$prop.direct.random, dat.trts, x$sep.trts)
  }
  else
    prop.common <- prop.random <- NULL
  
  
  ##
  ##
  ## (5) Calculate / extract indirect estimates
  ##
  ##
  x.direct.indirect <- x
  ##
  if (method == "SIDDE") {
    ind <- sidde(x.direct.indirect, sep.trts, verbose, warn, is.tictoc)
    ##
    x.direct.indirect$TE.indirect.common <- ind$TE.indirect.common
    x.direct.indirect$seTE.indirect.common <- ind$seTE.indirect.common
    ##
    if (!is.bin) {
      x.direct.indirect$TE.indirect.random <- ind$TE.indirect.random
      x.direct.indirect$seTE.indirect.random <- ind$seTE.indirect.random
    }
  }
  ##
  direct.indirect <- direct.indirect(x.direct.indirect, tol.direct)
  
  
  ##
  ##
  ## (6) Transform matrices to data frames
  ##
  ##
  m2d.f <- mat2dat.split(direct.indirect, "common", dat.trts)
  m2d.r <- mat2dat.split(direct.indirect, "random", dat.trts)
  
  
  ##
  ##
  ## (7) Return results
  ##
  ##
  res <- list(comparison = dat.trts$comparison,
              ##
              k = m2d.f$k,
              ##
              prop.common = prop.common,
              ##
              common = m2d.f$nma,
              direct.common = m2d.f$direct,
              indirect.common = m2d.f$indirect,
              compare.common = m2d.f$compare,
              ##
              prop.random = prop.random,
              ##
              random = m2d.r$nma,
              direct.random = m2d.r$direct,
              indirect.random = m2d.r$indirect,
              compare.random = m2d.r$compare,
              ##
              predict = m2d.r$predict,
              ##
              method = method,
              ##
              sm = x$sm,
              level.ma = x$level.ma,
              ##
              prediction = x$prediction,
              level.predict = x$level.predict,
              tau = x$tau,
              ##
              reference.group = reference.group,
              baseline.reference = baseline.reference,
              order = order,
              sep.trts = sep.trts,
              quote.trts = quote.trts,
              nchar.trts = x$nchar.trts,
              ##
              tol.direct = tol.direct,
              backtransf = backtransf,
              ##
              x = x,
              ##
              version = packageDescription("netmeta")$Version
              )
  ##
  if (method == "SIDDE" & is.tictoc)
    res$tictoc <- ind$tictoc
  ##
  ## Backward compatibility
  ##
  res$prop.fixed <- res$prop.common
  res$fixed <- res$common
  res$direct.fixed <- res$direct.common
  res$indirect.fixed <- res$indirect.common
  res$compare.fixed <- res$compare.common
  ##
  class(res) <-
    c("netsplit",
      if (is.bin & method == "SIDDE")
        "netsplit.netmetabin")
  
  res
}





#' @rdname netsplit
#' @method print netsplit
#' @export


print.netsplit <- function(x,
                           common = x$x$common,
                           random = x$x$random,
                           ##
                           show = "all",
                           overall = TRUE,
                           ci = FALSE,
                           test = show %in% c("all", "with.direct", "both"),
                           only.reference = FALSE,
                           ##
                           sortvar = NULL,
                           subset = NULL,
                           ##
                           nchar.trts = x$nchar.trts,
                           ##
                           digits = gs("digits"),
                           digits.stat = gs("digits.stat"),
                           digits.pval = gs("digits.pval"),
                           digits.prop = max(gs("digits.pval") - 2, 2),
                           ##
                           text.NA = ".",
                           backtransf = x$backtransf,
                           scientific.pval = gs("scientific.pval"),
                           big.mark = gs("big.mark"),
                           legend = TRUE,
                           ##
                           indent = TRUE,
                           warn.deprecated = gs("warn.deprecated"),
                           ##
                           ...) {
  
  ##
  ##
  ## (1) Check for netsplit object and upgrade object
  ##
  ##
  chkclass(x, "netsplit")
  x <- updateversion(x)
  ##
  is.bin <- inherits(x, "netsplit.netmetabin")
  ##  
  ## All individual results in a single row - be on the save side:
  ##
  oldopts <- options(width = 200)
  on.exit(options(oldopts))
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  ##
  args <- list(...)
  ##
  show <-
    deprecated(show, missing(show), args, "showall")
  if (is.logical(show))
    if (show)
      show <- "all"
    else
      show <- "both"
  ##
  show <- setchar(show, c("all", "both", "with.direct",
                          "direct.only", "indirect.only",
                          "reference.only"))
  ##
  if (show == "reference.only") {
    warning("Argument 'show = \"reference.only\" replaced with ",
            "'only.reference = TRUE'.",
            call. = FALSE)
    show <- "both"
    if (missing.only.reference)
      only.reference <- TRUE
  }
  ##
  chklogical(overall)
  chklogical(ci)
  chklogical(test)
  ##
  missing.only.reference <- missing(only.reference)
  if (!missing.only.reference)
    chklogical(only.reference)
  #
  # Catch sortvar and subset from data:
  #
  mc <- match.call()
  sfsp <- sys.frame(sys.parent())
  #
  error <- try(sortvar.x <- catch("sortvar", mc, x, sfsp), silent = TRUE)
  if (!any(class(error) == "try-error"))
    sortvar <- sortvar.x
  ##
  if (!is.null(sortvar)) {
    if (length(dim(sortvar)) == 2) {
      if (dim(sortvar)[1] != length(x$comparison))
        stop("Argument 'sortvar' must be of length ",
             length(x$comparison), ".",
             call. = FALSE)
      ##
      ## Set proportions to 0 or 1
      ##
      if (is.numeric(sortvar)) {
        sortvar[is.zero(abs(sortvar), n = 1000)] <- 0
        sortvar[is.zero(1 - abs(sortvar), n = 1000)] <-
          1 * sign(sortvar)[is.zero(1 - abs(sortvar), n = 1000)]
      }
      sortvar <- order(do.call(order, as.list(as.data.frame(sortvar))))
    }
    else
      chklength(sortvar, length(x$comparison),
                text = paste0("Argument 'sortvar' must be of length ",
                              length(x$comparison), "."))
    ##
    if (!is.numeric(sortvar))
      sortvar <- setchar(sortvar, x$comparison)
  }
  #
  error <- try(subset.x <- catch("subset", mc, x, sfsp), silent = TRUE)
  if (!any(class(error) == "try-error"))
    subset <- subset.x
  ##
  if (!is.null(subset)) {
    chklength(subset, length(x$comparison),
              text = paste0("Argument 'subset' must be of length ",
                            length(x$comparison), "."))
    if (!is.logical(subset))
      stop("Argument 'subset' must be a logical vector.")
  }
  ##
  if (is.null(nchar.trts))
    nchar.trts <- 666
  chknumeric(nchar.trts, length = 1)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)
  ##
  if (is.null(backtransf))
    backtransf <- TRUE
  chklogical(backtransf)
  chklogical(scientific.pval)
  chklogical(legend)
  chklogical(indent)
  ##
  ## Check for deprecated arguments in '...'
  ##
  fun <- "print.netmeta"
  ##
  chklogical(warn.deprecated)
  ##
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  common.logical <- common
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  random.logical <- random
  
  
  ##
  ##
  ## (3) Some additional settings and checks
  ##
  ##
  sm <- x$sm
  sm.lab <- sm
  ##
  relative <- is.relative.effect(sm)
  ##
  if (!backtransf & relative)
    sm.lab <- paste("log", sm, sep = "")
  ##
  if (!(sm.lab == "" | sm.lab == "log"))
    sm.lab <- paste("(", sm.lab, ") ", sep = "")
  else
    sm.lab <- ""
  ##
  level.ma <- x$level.ma
  ci.lab <- paste(100 * level.ma, "%-CI", sep ="")
  ##  
  random.available <- !is.null(x$random)
  ##
  if (!random.available & random) {
    warning("No results for random effects model available. ",
            "Argument 'random' set to FALSE.",
            call. = FALSE)
    ##
    random <- FALSE
  }
  
  
  if (show == "all")
    sel <- rep_len(TRUE, length(x$direct.common$TE))
  else if (show == "with.direct")
    sel <- !is.na(x$direct.common$TE)
  else if (show == "both")
    sel <- !is.na(x$direct.common$TE) & !is.na(x$indirect.common$TE)
  else if (show == "direct.only")
    sel <- !is.na(x$direct.common$TE) & is.na(x$indirect.common$TE)
  else if (show == "indirect.only")
    sel <- is.na(x$direct.common$TE) & !is.na(x$common$TE)
  ##
  if (sum(sel) == 0) {
    warning("No results for argument 'show = ", show, "'.",
            call. = FALSE)
    return(invisible(NULL))
  }
  ##
  if (only.reference) {
    if (x$reference.group == "") {
      warning("First treatment used as reference as argument ",
              "'reference.group' was unspecified in netsplit().",
              call. = FALSE)
      x$reference.group <-
        compsplit(x$comparison, x$sep.trts)[[1]][1]
    }
    ##
    sel.ref <-
      apply(!is.na(sapply(compsplit(x$comparison, x$sep.trts),
                          match, x$reference.group)), 2, sum) >= 1
    ##
    sel <- sel & sel.ref
  }
  ##
  comp <- x$comparison[sel]
  ##
  k <- x$k[sel]
  ##
  prop.common <- x$prop.common[sel]
  ##
  TE.common <- x$common$TE[sel]
  lower.common <- x$common$lower[sel]
  upper.common <- x$common$upper[sel]
  ##
  TE.direct.common <- x$direct.common$TE[sel]
  lower.direct.common <- x$direct.common$lower[sel]
  upper.direct.common <- x$direct.common$upper[sel]
  ##
  TE.indirect.common <- x$indirect.common$TE[sel]
  lower.indirect.common <- x$indirect.common$lower[sel]
  upper.indirect.common <- x$indirect.common$upper[sel]
  ##
  TE.compare.common <- x$compare.common$TE[sel]
  lower.compare.common <- x$compare.common$lower[sel]
  upper.compare.common <- x$compare.common$upper[sel]
  statistic.compare.common <- x$compare.common$statistic[sel]
  pval.compare.common <- x$compare.common$p[sel]
  ##
  if (random.available) {
    prop.random <- x$prop.random[sel]
    ##
    TE.random <- x$random$TE[sel]
    lower.random <- x$random$lower[sel]
    upper.random <- x$random$upper[sel]
    ##
    TE.direct.random <- x$direct.random$TE[sel]
    lower.direct.random <- x$direct.random$lower[sel]
    upper.direct.random <- x$direct.random$upper[sel]
    ##
    TE.indirect.random <- x$indirect.random$TE[sel]
    lower.indirect.random <- x$indirect.random$lower[sel]
    upper.indirect.random <- x$indirect.random$upper[sel]
    ##
    TE.compare.random <- x$compare.random$TE[sel]
    lower.compare.random <- x$compare.random$lower[sel]
    upper.compare.random <- x$compare.random$upper[sel]
    statistic.compare.random <- x$compare.random$statistic[sel]
    pval.compare.random <- x$compare.random$p[sel]
  }
  
  
  if (backtransf & relative) {
    TE.common <- exp(TE.common)
    lower.common <- exp(lower.common)
    upper.common <- exp(upper.common)
    ##
    TE.direct.common <- exp(TE.direct.common)
    lower.direct.common <- exp(lower.direct.common)
    upper.direct.common <- exp(upper.direct.common)
    ##
    TE.indirect.common <- exp(TE.indirect.common)
    lower.indirect.common <- exp(lower.indirect.common)
    upper.indirect.common <- exp(upper.indirect.common)
    ##
    TE.compare.common <- exp(TE.compare.common)
    lower.compare.common <- exp(lower.compare.common)
    upper.compare.common <- exp(upper.compare.common)
    ##
    if (random.available) {
      TE.random <- exp(TE.random)
      lower.random <- exp(lower.random)
      upper.random <- exp(upper.random)
      ##
      TE.direct.random <- exp(TE.direct.random)
      lower.direct.random <- exp(lower.direct.random)
      upper.direct.random <- exp(upper.direct.random)
      ##
      TE.indirect.random <- exp(TE.indirect.random)
      lower.indirect.random <- exp(lower.indirect.random)
      upper.indirect.random <- exp(upper.indirect.random)
      ##
      TE.compare.random <- exp(TE.compare.random)
      lower.compare.random <- exp(lower.compare.random)
      upper.compare.random <- exp(upper.compare.random)
    }
  }
  
  
  common <- list(comp = comp,
                 k = k,
                 prop = formatPT(prop.common, digits = digits.prop))
  names.common <- c("comparison", "k", "prop")
  ##
  if (overall) {
    common$TE.common <- formatN(TE.common, digits, text.NA = text.NA,
                                big.mark = big.mark)
    names.common <- c(names.common, "nma")
    if (ci) {
      common$ci.common <- formatCI(round(lower.common, digits),
                                   round(upper.common, digits))
      common$ci.common[is.na(common$ci.common)] <- text.NA
      names.common <- c(names.common, ci.lab)
    }
  }
  ##
  common$TE.direct.common <-
    formatN(TE.direct.common, digits, text.NA = text.NA,
            big.mark = big.mark)
  names.common <- c(names.common, "direct")
  if (ci) {
    common$ci.direct.common <-
      formatCI(round(lower.direct.common, digits),
               round(upper.direct.common, digits))
    common$ci.direct.common[is.na(common$ci.direct.common)] <- text.NA
    names.common <- c(names.common, ci.lab)
  }
  ##
  common$TE.indirect.common <-
    formatN(TE.indirect.common, digits,
            text.NA = text.NA, big.mark = big.mark)
  names.common <- c(names.common, "indir.")
  ##
  if (ci) {
    common$ci.indirect.common <-
      formatCI(round(lower.indirect.common, digits),
               round(upper.indirect.common, digits))
    common$ci.indirect.common[is.na(common$ci.indirect.common)] <- text.NA
    names.common <- c(names.common, ci.lab)
  }
  ##
  if (test) {
    common$diff <- formatN(TE.compare.common, digits, text.NA = text.NA,
                           big.mark = big.mark)
    names.common <-
      c(names.common, if (backtransf & relative) "RoR" else "Diff")
    if (ci) {
      common$ci.diff <- formatCI(round(lower.compare.common, digits),
                                 round(upper.compare.common, digits))
      common$ci.diff[is.na(common$ci.diff)] <- text.NA
      names.common <- c(names.common, ci.lab)
    }
    ##
    common$statistic <- formatN(statistic.compare.common, digits.stat,
                                big.mark = big.mark)
    common$statistic[common$statistic == "--"] <- text.NA
    common$p <- formatPT(pval.compare.common, digits = digits.pval,
                         scientific = scientific.pval)
    common$p[rmSpace(common$p) == "--"] <- text.NA
    names.common <- c(names.common, c("z", "p-value"))
  }
  common <- as.data.frame(common)
  names(common) <- names.common
  
  
  if (random.available) {
    random <- list(comp = comp,
                   k = k,
                   prop = formatPT(prop.random, digits = digits.prop))
    names.random <- c("comparison", "k", "prop")
    ##
    if (overall) {
      random$TE.random <- formatN(TE.random, digits, text.NA = text.NA,
                                  big.mark = big.mark)
      names.random <- c(names.random, "nma")
      if (ci) {
        random$ci.random <- formatCI(round(lower.random, digits),
                                     round(upper.random, digits))
        random$ci.random[is.na(random$ci.random)] <- text.NA
        names.random <- c(names.random, ci.lab)
      }
    }
    ##
    random$TE.direct.random <- formatN(TE.direct.random, digits,
                                       text.NA = text.NA,
                                       big.mark = big.mark)
    names.random <- c(names.random, "direct")
    if (ci) {
      random$ci.direct.random <- formatCI(round(lower.direct.random, digits),
                                          round(upper.direct.random, digits))
      random$ci.direct.random[is.na(random$ci.direct.random)] <- text.NA
      names.random <- c(names.random, ci.lab)
    }
    ##
    random$TE.indirect.random <- formatN(TE.indirect.random, digits,
                                         text.NA = text.NA,
                                         big.mark = big.mark)
    names.random <- c(names.random, "indir.")
    if (ci) {
      random$ci.indirect.random <-
        formatCI(round(lower.indirect.random, digits),
                 round(upper.indirect.random, digits))
      random$ci.indirect.random[is.na(random$ci.indirect.random)] <- text.NA
      names.random <- c(names.random, ci.lab)
    }
    ##
    if (test) {
      random$diff <- formatN(TE.compare.random, digits, text.NA = text.NA,
                             big.mark = big.mark)
      names.random <- c(names.random,
                        if (backtransf & relative) "RoR" else "Diff")
      if (ci) {
        random$ci.diff <- formatCI(round(lower.compare.random, digits),
                                   round(upper.compare.random, digits))
        random$ci.diff[is.na(random$ci.diff)] <- text.NA
        names.random <- c(names.random, ci.lab)
      }
      ##
      random$statistic <- formatN(statistic.compare.random, digits.stat,
                                  big.mark = big.mark)
      random$statistic[random$statistic == "--"] <- text.NA
      random$p <- formatPT(pval.compare.random, digits = digits.pval,
                           scientific = scientific.pval)
      random$p[rmSpace(random$p) == "--"] <- text.NA
      names.random <- c(names.random, c("z", "p-value"))
    }
    random <- as.data.frame(random)
    names(random) <- names.random
  }
  
  
  ## Do not print direct evidence proportion for node-splitting method
  ## or SIDDE
  ## 
  ##
  noprop <- is.bin | x$method == "SIDDE" | all(common$prop == "")
  ##
  if (noprop) {
    common <- common[, !(names(common) %in% "prop")]
    if (random.available)
      random <- random[, !(names(random) %in% "prop")]
  }
  

  if (!is.null(sortvar))
    sortvar <- sortvar[sel]
  #
  if (!is.null(subset))
    subset <- subset[sel]
  #
  if (!is.null(sortvar)) {
    o <- order(sortvar)
    #
    if (common.logical)
      common <- common[o, ]
    if (random.logical)
      random <- random[o, ]
    #
    if (!is.null(subset))
      subset <- subset[o]
  }
  #
  if (!is.null(subset)) {
    if (common.logical)
      common <- common[subset, ]
    if (random.logical)
      random <- random[subset, ]
  }
  
  
  if (common.logical | random.logical) {
    if (x$method == "SIDDE")
      cat("Separate indirect from direct design evidence (SIDDE)\n\n")
    else
      cat(paste("Separate indirect from direct evidence (SIDE)",
                "using back-calculation method\n\n"))
  }
  else
    legend <- FALSE
  
  
  if (common.logical) {
    cat("Common effects model: \n\n")
    common[is.na(common)] <- text.NA
    trts <- unique(sort(unlist(compsplit(common$comparison, x$sep.trts))))
    common$comparison <- comps(common$comparison, trts, x$sep.trts, nchar.trts)
    prmatrix(common, quote = FALSE, right = TRUE,
             rowlab = rep("", dim(common)[1]))
    if (random.logical)
      cat("\n")
  }
  ##
  if (random.logical) {
    cat("Random effects model: \n\n")
    random[is.na(random)] <- text.NA
    trts <- unique(sort(unlist(compsplit(random$comparison, x$sep.trts))))
    random$comparison <- comps(random$comparison, trts, x$sep.trts, nchar.trts)
    prmatrix(random, quote = FALSE, right = TRUE,
             rowlab = rep("", dim(random)[1]))
  }
  ##
  if (legend) {
    cat("\nLegend:\n")
    cat(" comparison - Treatment comparison\n")
    cat(paste0(" k", if (indent) "          " else " ",
               "- Number of studies providing direct evidence\n"))
    if (!noprop)
      cat(paste0(" prop", if (indent) "       " else " ",
                 "- Direct evidence proportion\n"))
    if (overall)
      cat(paste0(" nma", if (indent) "        " else " ",
                 "- Estimated treatment effect ", sm.lab,
                 "in network meta-analysis\n", sep = ""))
    cat(paste0(" direct", if (indent) "     " else " ",
               "- Estimated treatment effect ", sm.lab,
               "derived from direct evidence\n", sep = ""))
    cat(paste0(" indir.", if (indent) "     " else " ",
               "- Estimated treatment effect ", sm.lab,
               "derived from indirect evidence\n", sep = ""))
    if (test) {
      if (backtransf & relative)
        cat(paste0(" RoR", if (indent) "        " else " ",
                   "- Ratio of Ratios ",
                   "(direct versus indirect)\n"))
      else
        cat(paste0(" Diff", if (indent) "       " else " ",
                   "- Difference between direct and ",
                   "indirect treatment estimates\n"))
      ##
      cat(paste0(" z", if (indent) "          " else " ",
                 "- z-value of test for disagreement ",
                 "(direct versus indirect)\n"))
      cat(paste0(" p-value", if (indent) "    " else " ",
                 "- p-value of test for disagreement ",
                 "(direct versus indirect)\n"))
    }
    ##
    ## Add legend with abbreviated treatment labels
    ##
    legendabbr(trts, treats(trts, nchar.trts), TRUE, header = "\n")
  }
  
  
  invisible(NULL)
}
