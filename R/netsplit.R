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
#'   \code{"Back-calculation"} or \code{"SIDDE"}, can be
#'   abbreviated. See Details.
#' @param upper A logical indicating whether treatment comparisons
#'   should be selected from the lower or upper triangle of the
#'   treatment effect matrices (see list elements \code{TE.fixed} and
#'   \code{TE.random} in the \code{netmeta} object).
#' @param reference.group Reference treatment.
#' @param baseline.reference A logical indicating whether results
#'   should be expressed as comparisons of other treatments versus the
#'   reference treatment or vice versa. This argument is only
#'   considered if \code{reference.group} is not equal to \code{""}.
#' @param sep.trts A character string used in comparison names as
#'   separator between treatment labels, e.g., " vs ".
#' @param quote.trts A character used to print around treatment
#'   labels.
#' @param tol.direct A numeric defining the maximum deviation of the
#'   direct evidence proportion from 0 or 1 to classify a comparison
#'   as providing only indirect or direct evidence, respectively.
#' @param comb.fixed A logical indicating whether results for the
#'   fixed effects (common effects) network meta-analysis should be
#'   printed.
#' @param comb.random A logical indicating whether results for the
#'   random effects network meta-analysis should be printed.
#' @param show A character string indicating which comparisons should
#'   be printed (see Details).
#' @param overall A logical indicating whether estimates from network
#'   meta-analyis should be printed in addition to direct and indirect
#'   estimates.
#' @param ci A logical indicating whether confidence intervals should
#'   be printed in addition to treatment estimates.
#' @param test A logical indicating whether results of a test
#'   comparing direct and indirect estimates should be printed.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.zval Minimal number of significant digits for z-value
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
#' @param legend A logical indicating whether a legend show be
#'   printed.
#' @param warn A logical indicating whether warnings should be
#'   printed.
#' @param ... Additional arguments (ignored at the moment)
#' 
#' @details
#' A comparison of direct and indirect treatment estimates can serve
#' as check for consistency of network meta-analysis (Dias et al.,
#' 2010).
#' 
#' This function provides two methods to derive indirect estimates:
#' \itemize{
#' \item Separate Indirect from Direct Evidence (SIDE) using a
#'   back-calculation method. The \emph{direct evidence proportion} as
#'   described in König et al. (2013) is used in the calculation of
#'   the indirect evidence;
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
#' @return
#' An object of class \code{netsplit} with corresponding \code{print}
#' and \code{forest} functions. The object is a list containing the
#' following components:
#' \item{comb.fixed, comb.random}{As defined above.}
#' \item{comparison}{A vector with treatment comparisons.}
#' \item{prop.fixed, prop.random}{A vector with direct evidence
#'   proportions (fixed / random effects model).}
#' \item{fixed, random}{Results of network meta-analysis (fixed /
#'   random effects model), i.e., data frame with columns comparison,
#'   TE, seTE, lower, upper, z, and p.}
#' \item{direct.fixed, direct.random}{Network meta-analysis results
#'   based on direct evidence (fixed / random effects model), i.e.,
#'   data frame with columns comparison, TE, seTE, lower, upper, z,
#'   and p.}
#' \item{indirect.fixed, indirect.random}{Network meta-analysis
#'   results based on indirect evidence (fixed / random effects
#'   model), i.e., data frame with columns comparison, TE, seTE,
#'   lower, upper, z, and p.}
#' \item{compare.fixed, compare.random}{Comparison of direct and
#'   indirect evidence in network meta-analysis (fixed / random
#'   effects model), i.e., data frame with columns comparison, TE,
#'   seTE, lower, upper, z, and p.}
#' \item{sm}{A character string indicating underlying summary measure}
#' \item{level.comb}{The level used to calculate confidence intervals
#'   for pooled estimates.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}, Gerta
#'   Rücker \email{ruecker@@imbi.uni-freiburg.de}, Orestis Efthimiou
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
#' 1--21, https://doi.org/10.1002/sim.8158
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
#'                studlab = author, data = Woods2010, sm = "OR")
#' #
#' net1 <- netmeta(p1)
#' #
#' print(netsplit(net1), digits = 2)
#' print(netsplit(net1), digits = 2,
#'       backtransf = FALSE, comb.fixed = FALSE)
#' 
#' data(Senn2013)
#' #
#' net2 <- netmeta(TE, seTE, treat1.long, treat2.long,
#'                 studlab, data = Senn2013)
#' #
#' print(netsplit(net2), digits = 2)
#' # Layout of Puhan et al. (2014), Table 1
#' print(netsplit(net2), digits = 2, ci = TRUE, test = FALSE)
#' 
#' \dontrun{
#' data(Dong2013)
#' p3 <- pairwise(treatment, death, randomized, studlab = id,
#'                data = Dong2013, sm = "OR")
#' net3  <- netmetabin(p3)
#' netsplit(net3)
#' }
#' 
#' @rdname netsplit
#' @export netsplit


netsplit <- function(x, method,
                     upper = TRUE,
                     reference.group = x$reference.group,
                     baseline.reference = x$baseline.reference,
                     sep.trts = x$sep.trts, quote.trts = "",
                     tol.direct = 0.0005,
                     comb.fixed = x$comb.fixed,
                     comb.random = x$comb.random,
                     backtransf = x$backtransf,
                     warn = FALSE) {
  
  
  meta:::chkclass(x, "netmeta")
  ##
  is.bin <- inherits(x, "netmetabin")
  ##
  if (!missing(method))
    method <- meta:::setchar(method, c("Back-calculation", "SIDDE"))
  else {
    if (is.bin)
      method <- "SIDDE"
    else
      method <- "Back-calculation"
  }
  ##
  meta:::chklogical(upper)
  meta:::chklogical(baseline.reference)
  meta:::chkchar(sep.trts)
  meta:::chkchar(quote.trts)
  meta:::chknumeric(tol.direct, min = 0, length = 1)
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  if (!is.null(backtransf))
    meta:::chklogical(backtransf)
  meta:::chklogical(warn)
  
  
  seq.comps <- rownames(x$Cov.fixed)
  ##
  trts <- matrix(unlist(compsplit(seq.comps, x$sep.trts)),
                 ncol = 2, byrow = TRUE)
  trts <- as.data.frame(trts, stringsAsFactors = FALSE)
  names(trts) <- c("treat1", "treat2")
  ##
  if (!upper) {
    ##
    ## Comparison names are column:row (and must be switched)
    ##
    t1 <- trts$treat1
    trts$treat1 <- trts$treat2
    trts$treat2 <- t1
  }
  
  
  ##
  ## Number of studies with direct comparisons
  ##
  k <- lowertri(x$A.matrix)
  
  
  ##
  ## Change treatment order if
  ## - reference group is specified, i.e., unequal to ""
  ## - reference group is first treatment (argument 'baseline.reference' is TRUE)
  ## - reference group is second treatment (argument 'baseline.reference' is FALSE)
  ##
  wo <- rep_len(FALSE, length(seq.comps))
  ##
  if (reference.group != "") {
    reference.group <- setref(reference.group, colnames(x$TE.fixed))
    ##
    if (baseline.reference)
      wo <- trts$treat1 == reference.group
    else
      wo <- trts$treat2 == reference.group
  }
  else
    if (!missing(baseline.reference))
      warning("Argument 'baseline.reference' ignored as reference group is not defined (argument 'reference.group').")
  ##
  if (any(wo)) {
    t1.wo <- trts$treat1[wo]
    trts$treat1[wo] <- trts$treat2[wo]
    trts$treat2[wo] <- t1.wo
  }
  ##
  comparison <- as.character(interaction(paste(quote.trts, trts$treat1,
                                               quote.trts, sep = ""),
                                         paste(quote.trts, trts$treat2,
                                               quote.trts, sep = ""),
                                         sep = sep.trts))
  
  
  ##
  ##
  ## Back-calculation method
  ## - based on direct evidence proportion (König et al. (2013)
  ##
  ##
  if (method == "Back-calculation") {
    ##
    ## Indirect estimate is NA if only direct evidence is available
    ##
    sel.one.fixed <- abs(x$P.fixed - 1) < tol.direct
    ##
    TE.indirect.fixed <- x$TE.indirect.fixed
    seTE.indirect.fixed <- x$seTE.indirect.fixed
    lower.indirect.fixed <- x$lower.indirect.fixed
    upper.indirect.fixed <- x$upper.indirect.fixed
    zval.indirect.fixed <- x$zval.indirect.fixed
    pval.indirect.fixed <- x$pval.indirect.fixed
    ##
    TE.indirect.fixed[sel.one.fixed] <- NA
    seTE.indirect.fixed[sel.one.fixed] <- NA
    lower.indirect.fixed[sel.one.fixed] <- NA
    upper.indirect.fixed[sel.one.fixed] <- NA
    zval.indirect.fixed[sel.one.fixed] <- NA
    pval.indirect.fixed[sel.one.fixed] <- NA
    ##
    sel.one.random <- abs(x$P.random - 1) < tol.direct
    ##
    TE.indirect.random <- x$TE.indirect.random
    seTE.indirect.random <- x$seTE.indirect.random
    lower.indirect.random <- x$lower.indirect.random
    upper.indirect.random <- x$upper.indirect.random
    zval.indirect.random <- x$zval.indirect.random
    pval.indirect.random <- x$pval.indirect.random
    ##
    TE.indirect.random[sel.one.random] <- NA
    seTE.indirect.random[sel.one.random] <- NA
    lower.indirect.random[sel.one.random] <- NA
    upper.indirect.random[sel.one.random] <- NA
    zval.indirect.random[sel.one.random] <- NA
    pval.indirect.random[sel.one.random] <- NA
  }
  
  
  ##
  ##
  ## Separate Indirect from Direct Design Evidence (SIDDE)
  ##
  ##
  if (method == "SIDDE") {
    ##
    if (is.null(x$data))
      stop("SIDDE method only available for network meta-analysis objects ",
           "created with argument 'keepdata' equal to TRUE.")
    ##
    dat <- x$data
    dat <- dat[order(dat$.studlab, dat$.treat1, dat$.treat2), ]
    ##
    if (!is.null(dat$.subset))
      dat <- dat[dat$.subset, , drop = FALSE]
    ##
    if (!is.null(dat$.drop))
      dat <- dat[!dat$.drop, , drop = FALSE]
    ##
    ## Determine comparisons with direct evidence
    ##
    idx.d <- which(!is.na(x$TE.direct.fixed), arr.ind = TRUE)
    idx.d <- idx.d[idx.d[, 1] < idx.d[, 2], , drop = FALSE]
    ##
    rownames(idx.d) <- seq_len(nrow(idx.d))
    idx1 <- idx.d[, 1]
    idx2 <- idx.d[, 2]
    ##
    n.comps <- nrow(idx.d)
    ##
    trts <- x$trts
    ##
    ## Perform network meta-analyses for indirect evidence
    ## (by dropping one direct comparison at a time)
    ##
    TE.indirect.fixed <- x$TE.direct.fixed
    TE.indirect.fixed[!is.na(TE.indirect.fixed)] <- NA
    seTE.indirect.fixed <- TE.indirect.fixed
    ##
    TE.indirect.random <- x$TE.direct.random
    TE.indirect.random[!is.na(TE.indirect.random)] <- NA
    seTE.indirect.random <- TE.indirect.random
    ##
    for (i in seq_len(n.comps)) {
      ##
      idx1.i <- idx1[i]
      idx2.i <- idx2[i]
      ##
      drop.i <-
        (dat$.treat1 == trts[idx1.i] & dat$.treat2 == trts[idx2.i]) |
        (dat$.treat2 == trts[idx1.i] & dat$.treat1 == trts[idx2.i])
      ##
      ## Studies (potentially with multi-arm studies) to drop from
      ## calculation of indirect estimate
      ##
      drop.studies <- unique(dat$studlab[drop.i])
      ##
      ## Drop studies
      ##
      dat.i <- dat[!(dat$studlab %in% drop.studies), , drop = FALSE]
      dat.i$.design <- NULL
      ##
      con <- netconnection(dat.i$treat1, dat.i$treat2, dat.i$studlab)
      ##
      if (con$n.subnets == 1) {
        ##
        if (is.bin)
          net.i <- netmetabin(dat.i$.event1, dat.i$.n1,
                              dat.i$.event2, dat.i$.n2,
                              dat.i$.treat1, dat.i$.treat2,
                              dat.i$.studlab,
                              data = dat.i,
                              sm = x$sm, method = x$method,
                              comb.fixed = comb.fixed,
                              comb.random = comb.random,
                              warn = warn)
        else
          net.i <- netmeta(dat.i$.TE, dat.i$.seTE,
                           dat.i$.treat1, dat.i$.treat2,
                           dat.i$.studlab,
                           data = dat.i,
                           comb.fixed = comb.fixed,
                           comb.random = comb.random,
                           warn = warn)
        ##
        if (trts[idx1.i] %in% rownames(net.i$TE.fixed) &
            trts[idx2.i] %in% colnames(net.i$TE.fixed)) {
          TE.indirect.fixed[idx1.i, idx2.i] <-
            net.i$TE.fixed[trts[idx1.i], trts[idx2.i]]
          TE.indirect.fixed[idx2.i, idx1.i] <-
            net.i$TE.fixed[trts[idx2.i], trts[idx1.i]]
          ##
          seTE.indirect.fixed[idx1.i, idx2.i] <-
            seTE.indirect.fixed[idx2.i, idx1.i] <-
            net.i$seTE.fixed[trts[idx1.i], trts[idx2.i]]
        }
        ##
        if (!is.bin) {
          if (trts[idx1.i] %in% rownames(net.i$TE.random) &
              trts[idx2.i] %in% colnames(net.i$TE.random)) {
            TE.indirect.random[idx1.i, idx2.i] <-
              net.i$TE.random[trts[idx1.i], trts[idx2.i]]
            TE.indirect.random[idx2.i, idx1.i] <-
              net.i$TE.random[trts[idx2.i], trts[idx1.i]]
            ##
            seTE.indirect.random[idx1.i, idx2.i] <-
              seTE.indirect.random[idx2.i, idx1.i] <-
              net.i$seTE.random[trts[idx1.i], trts[idx2.i]]
          }
        }
      }
    }
    ##
    ci.if <- ci(TE.indirect.fixed, seTE.indirect.fixed, x$level.comb)
    ##
    lower.indirect.fixed <- ci.if$lower
    upper.indirect.fixed <- ci.if$upper
    zval.indirect.fixed <- ci.if$z
    pval.indirect.fixed <- ci.if$p
    ##
    if (!is.bin) {
      ci.ir <- ci(TE.indirect.random, seTE.indirect.random, x$level.comb)
      ##
      lower.indirect.random <- ci.ir$lower
      upper.indirect.random <- ci.ir$upper
      zval.indirect.random <- ci.ir$z
      pval.indirect.random <- ci.ir$p
    }
  }
  
  
  ##
  ## Set direct evidence estimates to NA if only indirect evidence is
  ## available
  ##
  TE.direct.fixed <- x$TE.direct.fixed
  seTE.direct.fixed <- x$seTE.direct.fixed
  lower.direct.fixed <- x$lower.direct.fixed
  upper.direct.fixed <- x$upper.direct.fixed
  zval.direct.fixed <- x$zval.direct.fixed
  pval.direct.fixed <- x$pval.direct.fixed
  ##
  if (!is.null(x$P.fixed)) {
    sel.zero.fixed <- abs(x$P.fixed) < tol.direct
    ##
    TE.direct.fixed[sel.zero.fixed]    <- NA
    seTE.direct.fixed[sel.zero.fixed]  <- NA
    lower.direct.fixed[sel.zero.fixed] <- NA
    upper.direct.fixed[sel.zero.fixed] <- NA
    zval.direct.fixed[sel.zero.fixed]  <- NA
    pval.direct.fixed[sel.zero.fixed] <- NA
  }
  ##
  TE.direct.random <- x$TE.direct.random
  seTE.direct.random <- x$seTE.direct.random
  lower.direct.random <- x$lower.direct.random
  upper.direct.random <- x$upper.direct.random
  zval.direct.random <- x$zval.direct.random
  pval.direct.random <- x$pval.direct.random
  ##
  if (!is.null(x$P.random)) {
    sel.zero.random <- abs(x$P.random) < tol.direct
    ##
    TE.direct.random[sel.zero.random] <- NA
    seTE.direct.random[sel.zero.random] <- NA
    lower.direct.random[sel.zero.random] <- NA
    upper.direct.random[sel.zero.random] <- NA
    zval.direct.random[sel.zero.random] <- NA
    pval.direct.random[sel.zero.random] <- NA
  }
  
  
  ##
  ## Fixed effects model
  ##
  fixed.low <- data.frame(comparison,
                          TE = lowertri(x$TE.fixed),
                          seTE = lowertri(x$seTE.fixed),
                          lower = lowertri(x$lower.fixed),
                          upper = lowertri(x$upper.fixed),
                          z = lowertri(x$zval.fixed),
                          p = lowertri(x$pval.fixed),
                          stringsAsFactors = FALSE)
  ##
  direct.fixed.low <- data.frame(comparison,
                                 TE = lowertri(TE.direct.fixed),
                                 seTE = lowertri(seTE.direct.fixed),
                                 lower = lowertri(lower.direct.fixed),
                                 upper = lowertri(upper.direct.fixed),
                                 z = lowertri(zval.direct.fixed),
                                 p = lowertri(pval.direct.fixed),
                                 stringsAsFactors = FALSE)
  ##
  indirect.fixed.low <- data.frame(comparison,
                                   TE = lowertri(TE.indirect.fixed),
                                   seTE = lowertri(seTE.indirect.fixed),
                                   lower = lowertri(lower.indirect.fixed),
                                   upper = lowertri(upper.indirect.fixed),
                                   z = lowertri(zval.indirect.fixed),
                                   p = lowertri(pval.indirect.fixed),
                                   stringsAsFactors = FALSE)
  ##
  fixed.upp <- data.frame(comparison,
                          TE = uppertri(x$TE.fixed),
                          seTE = uppertri(x$seTE.fixed),
                          lower = uppertri(x$lower.fixed),
                          upper = uppertri(x$upper.fixed),
                          z = uppertri(x$zval.fixed),
                          p = uppertri(x$pval.fixed),
                          stringsAsFactors = FALSE)
  ##
  direct.fixed.upp <- data.frame(comparison,
                                 TE = uppertri(TE.direct.fixed),
                                 seTE = uppertri(seTE.direct.fixed),
                                 lower = uppertri(lower.direct.fixed),
                                 upper = uppertri(upper.direct.fixed),
                                 z = uppertri(zval.direct.fixed),
                                 p = uppertri(pval.direct.fixed),
                                 stringsAsFactors = FALSE)
  ##
  indirect.fixed.upp <- data.frame(comparison,
                                   TE = uppertri(TE.indirect.fixed),
                                   seTE = uppertri(seTE.indirect.fixed),
                                   lower = uppertri(lower.indirect.fixed),
                                   upper = uppertri(upper.indirect.fixed),
                                   z = uppertri(zval.indirect.fixed),
                                   p = uppertri(pval.indirect.fixed),
                                   stringsAsFactors = FALSE)
  ##
  if (!upper) {
    fixed <- fixed.low
    direct.fixed <- direct.fixed.low
    indirect.fixed <- indirect.fixed.low
    ##
    if (any(wo)) {
      fixed[wo, ] <- fixed.upp[wo, ]
      direct.fixed[wo, ] <- direct.fixed.upp[wo, ]
      indirect.fixed[wo, ] <- indirect.fixed.upp[wo, ]
    }
  }
  else {
    fixed <- fixed.upp
    direct.fixed <- direct.fixed.upp
    indirect.fixed <- indirect.fixed.upp
    ##
    if (any(wo)) {
      fixed[wo, ] <- fixed.low[wo, ]
      direct.fixed[wo, ] <- direct.fixed.low[wo, ]
      indirect.fixed[wo, ] <- indirect.fixed.low[wo, ]
    }
  }
  ##
  m.fixed <- metagen(direct.fixed$TE - indirect.fixed$TE,
                     sqrt(direct.fixed$seTE^2 + indirect.fixed$seTE^2),
                     level = x$level.comb)
  ##
  compare.fixed <- data.frame(comparison,
                              TE = m.fixed$TE,
                              seTE = m.fixed$seTE,
                              lower = m.fixed$lower,
                              upper = m.fixed$upper,
                              z = m.fixed$zval,
                              p = m.fixed$pval,
                              stringsAsFactors = FALSE)
  ##
  sel.k0 <- k == 0
  vars <- c("TE", "seTE", "lower", "upper", "z", "p")
  ##
  indirect.fixed[sel.k0, vars] <- fixed[sel.k0, vars]
  
  
  if (!is.bin) {
    ##
    ## Random effects model
    ##
    random.low <- data.frame(comparison,
                             TE = lowertri(x$TE.random),
                             seTE = lowertri(x$seTE.random),
                             lower = lowertri(x$lower.random),
                             upper = lowertri(x$upper.random),
                             z = lowertri(x$zval.random),
                             p = lowertri(x$pval.random),
                             stringsAsFactors = FALSE)
    ##
    direct.random.low <- data.frame(comparison,
                                    TE = lowertri(TE.direct.random),
                                    seTE = lowertri(seTE.direct.random),
                                    lower = lowertri(lower.direct.random),
                                    upper = lowertri(upper.direct.random),
                                    z = lowertri(zval.direct.random),
                                    p = lowertri(pval.direct.random),
                                    stringsAsFactors = FALSE)
    ##
    indirect.random.low <- data.frame(comparison,
                                      TE = lowertri(TE.indirect.random),
                                      seTE = lowertri(seTE.indirect.random),
                                      lower = lowertri(lower.indirect.random),
                                      upper = lowertri(upper.indirect.random),
                                      z = lowertri(zval.indirect.random),
                                      p = lowertri(pval.indirect.random),
                                      stringsAsFactors = FALSE)
    ##
    predict.low <- data.frame(comparison,
                              lower = lowertri(x$lower.predict),
                              upper = lowertri(x$upper.predict),
                              stringsAsFactors = FALSE)
    ##
    random.upp <- data.frame(comparison,
                             TE = uppertri(x$TE.random),
                             seTE = uppertri(x$seTE.random),
                             lower = uppertri(x$lower.random),
                             upper = uppertri(x$upper.random),
                             z = uppertri(x$zval.random),
                             p = uppertri(x$pval.random),
                             stringsAsFactors = FALSE)
    ##
    direct.random.upp <- data.frame(comparison,
                                    TE = uppertri(TE.direct.random),
                                    seTE = uppertri(seTE.direct.random),
                                    lower = uppertri(lower.direct.random),
                                    upper = uppertri(upper.direct.random),
                                    z = uppertri(zval.direct.random),
                                    p = uppertri(pval.direct.random),
                                    stringsAsFactors = FALSE)
    ##
    indirect.random.upp <- data.frame(comparison,
                                      TE = uppertri(TE.indirect.random),
                                      seTE = uppertri(seTE.indirect.random),
                                      lower = uppertri(lower.indirect.random),
                                      upper = uppertri(upper.indirect.random),
                                      z = uppertri(zval.indirect.random),
                                      p = uppertri(pval.indirect.random),
                                      stringsAsFactors = FALSE)
    ##
    predict.upp <- data.frame(comparison,
                              lower = uppertri(x$lower.predict),
                              upper = uppertri(x$upper.predict),
                              stringsAsFactors = FALSE)
    ##
    if (!upper) {
      random <- random.low
      direct.random <- direct.random.low
      indirect.random <- indirect.random.low
      predict <- predict.low
      ##
      if (any(wo)) {
        random[wo, ] <- random.upp[wo, ]
        direct.random[wo, ] <- direct.random.upp[wo, ]
        indirect.random[wo, ] <- indirect.random.upp[wo, ]
        predict[wo, ] <- predict.upp[wo, ]
      }
    }
    else {
      random <- random.upp
      direct.random <- direct.random.upp
      indirect.random <- indirect.random.upp
      predict <- predict.upp
      ##
      if (any(wo)) {
        random[wo, ] <- random.low[wo, ]
        direct.random[wo, ] <- direct.random.low[wo, ]
        indirect.random[wo, ] <- indirect.random.low[wo, ]
        predict[wo, ] <- predict.low[wo, ]
      }
    }
    ##
    m.random <- metagen(direct.random$TE - indirect.random$TE,
                        sqrt(direct.random$seTE^2 + indirect.random$seTE^2),
                        level = x$level.comb)
    ##
    compare.random <- data.frame(comparison,
                                 TE = m.random$TE,
                                 seTE = m.random$seTE,
                                 lower = m.random$lower,
                                 upper = m.random$upper,
                                 z = m.random$zval,
                                 p = m.random$pval,
                                 stringsAsFactors = FALSE)
    ##
    indirect.random[sel.k0, vars] <- random[sel.k0, vars]
  }
  else {
    if (!upper)
      random <- fixed.low
    else
      random <- fixed.upp
    ##
    random[!is.na(random)] <- NA
    random$comparison <- comparison
    direct.random <- indirect.random <- compare.random <- random
    predict <- random[, c("comparison", "lower", "upper")]
  }
  
  
  res <- list(comparison = comparison,
              ##
              k = k,
              ##
              prop.fixed =
                if (is.bin) NULL else x$prop.direct.fixed[seq.comps],
              fixed = fixed,
              direct.fixed = direct.fixed,
              indirect.fixed = indirect.fixed,
              compare.fixed = compare.fixed,
              ##
              prop.random =
                if (is.bin) NULL else x$prop.direct.random[seq.comps],
              random = random,
              direct.random = direct.random,
              indirect.random = indirect.random,
              compare.random = compare.random,
              predict = predict,
              ##
              method = method,
              ##
              sm = x$sm,
              level.comb = x$level.comb,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              ##
              prediction = x$prediction,
              level.predict = x$level.predict,
              tau = x$tau,
              ##
              reference.group = reference.group,
              baseline.reference = baseline.reference,
              sep.trts = sep.trts,
              quote.trts = quote.trts,
              ##
              tol.direct = tol.direct,
              backtransf = backtransf,
              ##
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- c("netsplit",
                  if (is.bin & method == "SIDDE") "netsplit.netmetabin")
  
  res
}





#' @rdname netsplit
#' @method print netsplit
#' @export
#' @export print.netsplit


print.netsplit <- function(x,
                           comb.fixed = x$comb.fixed,
                           comb.random = x$comb.random,
                           show = "all",
                           overall = TRUE,
                           ci = FALSE,
                           test = show %in% c("all", "with.direct", "both"),
                           digits = gs("digits"),
                           digits.zval = gs("digits.zval"),
                           digits.pval = gs("digits.pval"),
                           digits.prop = max(gs("digits.pval") - 2, 2),
                           text.NA = ".",
                           backtransf = x$backtransf,
                           scientific.pval = gs("scientific.pval"),
                           big.mark = gs("big.mark"),
                           legend = TRUE,
                           ...) {


  meta:::chkclass(x, "netsplit")
  ##
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  formatCI <- meta:::formatCI
  formatN <- meta:::formatN
  formatPT <- meta:::formatPT
  is.relative.effect <- meta:::is.relative.effect
  rmSpace <- meta:::rmSpace
  setchar <- meta:::setchar
  
  
  is.bin <- inherits(x, "netsplit.netmetabin")
  
  
  ## All individual results in a single row - be on the save side:
  ##
  oldopts <- options(width = 200)
  on.exit(options(oldopts))


  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(overall)
  chklogical(ci)
  chklogical(test)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.zval, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)
  ##
  if (is.null(backtransf))
    backtransf <- TRUE
  chklogical(backtransf)
  chklogical(scientific.pval)
  chklogical(legend)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  ## Check whether first argument is a list. In this case only use
  ## this list as input.
  if (length(args) > 0 && is.list(args[[1]]))
    args <- args[[1]]
  ##
  additional.arguments <- names(args)
  ##
  if (length(additional.arguments) > 0) {
    if (!is.na(charmatch("showa", additional.arguments)))
      if (!missing(show))
        warning("Deprecated argument 'showall' ignored as argument 'show' is also provided.")
      else {
        warning("Deprecated argument 'showall' has been replaced by argument 'show'.")
        show <- args[[charmatch("showa", additional.arguments)]]
        if (show)
          show <- "all"
        else
          show <- "both"
      }
  }
  ##
  show <- setchar(show, c("all", "both", "with.direct", "direct.only", "indirect.only"))


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


  level.comb <- x$level.comb
  ci.lab <- paste(100 * level.comb, "%-CI", sep ="")
  
  
  random.available <- !is.null(x$random)
  ##
  if (!random.available & comb.random) {
    warning("No results for random effects model available. ",
            "Argument 'comb.random' set to FALSE.",
            call. = FALSE)
    ##
    comb.random <- FALSE
  }
  
  
  if (show == "all")
    sel <- rep_len(TRUE, length(x$direct.fixed$TE))
  else if (show == "with.direct")
    sel <- !is.na(x$direct.fixed$TE)
  else if (show == "both")
    sel <- !is.na(x$direct.fixed$TE) & !is.na(x$indirect.fixed$TE)
  else if (show == "direct.only")
    sel <- !is.na(x$direct.fixed$TE) & is.na(x$indirect.fixed$TE)
  else if (show == "indirect.only")
    sel <- is.na(x$direct.fixed$TE) & !is.na(x$fixed$TE)
  ##
  comp <- x$comparison[sel]
  ##
  k <- x$k[sel]
  ##
  prop.fixed <- x$prop.fixed[sel]
  ##
  TE.fixed <- x$fixed$TE[sel]
  lower.fixed <- x$fixed$lower[sel]
  upper.fixed <- x$fixed$upper[sel]
  ##
  TE.direct.fixed <- x$direct.fixed$TE[sel]
  lower.direct.fixed <- x$direct.fixed$lower[sel]
  upper.direct.fixed <- x$direct.fixed$upper[sel]
  ##
  TE.indirect.fixed <- x$indirect.fixed$TE[sel]
  lower.indirect.fixed <- x$indirect.fixed$lower[sel]
  upper.indirect.fixed <- x$indirect.fixed$upper[sel]
  ##
  TE.compare.fixed <- x$compare.fixed$TE[sel]
  lower.compare.fixed <- x$compare.fixed$lower[sel]
  upper.compare.fixed <- x$compare.fixed$upper[sel]
  zval.compare.fixed <- x$compare.fixed$z[sel]
  pval.compare.fixed <- x$compare.fixed$p[sel]
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
    zval.compare.random <- x$compare.random$z[sel]
    pval.compare.random <- x$compare.random$p[sel]
  }
  
  
  if (backtransf & relative) {
    TE.fixed <- exp(TE.fixed)
    lower.fixed <- exp(lower.fixed)
    upper.fixed <- exp(upper.fixed)
    ##
    TE.direct.fixed <- exp(TE.direct.fixed)
    lower.direct.fixed <- exp(lower.direct.fixed)
    upper.direct.fixed <- exp(upper.direct.fixed)
    ##
    TE.indirect.fixed <- exp(TE.indirect.fixed)
    lower.indirect.fixed <- exp(lower.indirect.fixed)
    upper.indirect.fixed <- exp(upper.indirect.fixed)
    ##
    TE.compare.fixed <- exp(TE.compare.fixed)
    lower.compare.fixed <- exp(lower.compare.fixed)
    upper.compare.fixed <- exp(upper.compare.fixed)
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
  
  
  fixed <- list(comp = comp,
                k = k,
                prop = formatPT(prop.fixed, digits = digits.prop))
  names.fixed <- c("comparison", "k", "prop")
  ##
  if (overall) {
    fixed$TE.fixed <- formatN(TE.fixed, digits, text.NA = text.NA,
                              big.mark = big.mark)
    names.fixed <- c(names.fixed, "nma")
    if (ci) {
      fixed$ci.fixed <- formatCI(round(lower.fixed, digits),
                                 round(upper.fixed, digits))
      fixed$ci.fixed[is.na(fixed$ci.fixed)] <- text.NA
      names.fixed <- c(names.fixed, ci.lab)
    }
  }
  ##
  fixed$TE.direct.fixed <- formatN(TE.direct.fixed, digits, text.NA = text.NA,
                                   big.mark = big.mark)
  names.fixed <- c(names.fixed, "direct")
  if (ci) {
    fixed$ci.direct.fixed <- formatCI(round(lower.direct.fixed, digits),
                                      round(upper.direct.fixed, digits))
    fixed$ci.direct.fixed[is.na(fixed$ci.direct.fixed)] <- text.NA
    names.fixed <- c(names.fixed, ci.lab)
  }
  ##
  fixed$TE.indirect.fixed <- formatN(TE.indirect.fixed, digits,
                                     text.NA = text.NA, big.mark = big.mark)
  names.fixed <- c(names.fixed, "indir.")
  ##
  if (ci) {
    fixed$ci.indirect.fixed <- formatCI(round(lower.indirect.fixed, digits),
                                        round(upper.indirect.fixed, digits))
    fixed$ci.indirect.fixed[is.na(fixed$ci.indirect.fixed)] <- text.NA
    names.fixed <- c(names.fixed, ci.lab)
  }
  ##
  if (test) {
    fixed$diff <- formatN(TE.compare.fixed, digits, text.NA = text.NA,
                          big.mark = big.mark)
    names.fixed <- c(names.fixed, if (backtransf & relative) "RoR" else "Diff")
    if (ci) {
      fixed$ci.diff <- formatCI(round(lower.compare.fixed, digits),
                                round(upper.compare.fixed, digits))
      fixed$ci.diff[is.na(fixed$ci.diff)] <- text.NA
      names.fixed <- c(names.fixed, ci.lab)
    }
    ##
    fixed$z <- formatN(zval.compare.fixed, digits.zval,
                       big.mark = big.mark)
    fixed$z[fixed$z == "--"] <- text.NA
    fixed$p <- formatPT(pval.compare.fixed, digits = digits.pval,
                        scientific = scientific.pval)
    fixed$p[rmSpace(fixed$p) == "--"] <- text.NA
    names.fixed <- c(names.fixed, c("z", "p-value"))
  }
  fixed <- as.data.frame(fixed)
  names(fixed) <- names.fixed
  
  
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
      random$ci.indirect.random <- formatCI(round(lower.indirect.random, digits),
                                            round(upper.indirect.random, digits))
      random$ci.indirect.random[is.na(random$ci.indirect.random)] <- text.NA
      names.random <- c(names.random, ci.lab)
    }
    ##
    if (test) {
      random$diff <- formatN(TE.compare.random, digits, text.NA = text.NA,
                             big.mark = big.mark)
      names.random <- c(names.random, if (backtransf & relative) "RoR" else "Diff")
      if (ci) {
        random$ci.diff <- formatCI(round(lower.compare.random, digits),
                                   round(upper.compare.random, digits))
        random$ci.diff[is.na(random$ci.diff)] <- text.NA
        names.random <- c(names.random, ci.lab)
      }
      ##
      random$z <- formatN(zval.compare.random, digits.zval,
                          big.mark = big.mark)
      random$z[random$z == "--"] <- text.NA
      random$p <- formatPT(pval.compare.random, digits = digits.pval,
                           scientific = scientific.pval)
      random$p[rmSpace(random$p) == "--"] <- text.NA
      names.random <- c(names.random, c("z", "p-value"))
    }
    random <- as.data.frame(random)
    names(random) <- names.random
  }
  
  
  ## Do not print direct evidence proportion for SIDDE
  ##
  noprop <- is.bin | x$method == "SIDDE" | all(fixed$prop == "")
  if (noprop) {
    fixed <- fixed[, !(names(fixed) %in% "prop")]
    if (random.available)
      random <- random[, !(names(random) %in% "prop")]
  }
  
  
  if (x$method == "SIDDE")
    cat("Separate indirect from direct design evidence (SIDDE)\n\n")
  else
    cat(paste("Separate indirect from direct evidence (SIDE)",
              "using back-calculation method\n\n"))
  
  
  if (comb.fixed) {
    cat("Fixed effects model: \n\n")
    fixed[is.na(fixed)] <- text.NA
    prmatrix(fixed, quote = FALSE, right = TRUE,
             rowlab = rep("", dim(fixed)[1]))
    if (comb.random)
      cat("\n")
  }
  ##
  if (comb.random) {
    cat("Random effects model: \n\n")
    random[is.na(random)] <- text.NA
    prmatrix(random, quote = FALSE, right = TRUE,
             rowlab = rep("", dim(random)[1]))
  }
  ##
  if (legend) {
    cat("\nLegend:\n")
    cat(" comparison - Treatment comparison\n")
    cat(" k          - Number of studies providing direct evidence\n")
    if (!noprop)
      cat(" prop       - Direct evidence proportion\n")
    if (overall)
      cat(paste(" nma        - Estimated treatment effect ", sm.lab,
                "in network meta-analysis\n", sep = ""))
    cat(paste(" direct     - Estimated treatment effect ", sm.lab,
              "derived from direct evidence\n", sep = ""))
    cat(paste(" indir.     - Estimated treatment effect ", sm.lab,
              "derived from indirect evidence\n", sep = ""))
    if (test) {
      if (backtransf & relative)
        cat(" RoR        - Ratio of Ratios (direct versus indirect)\n")
      else
        cat(" Diff       - Difference between direct and indirect treatment estimates\n")
      cat(" z          - z-value of test for disagreement (direct versus indirect)\n")
      cat(" p-value    - p-value of test for disagreement (direct versus indirect)\n")
    }
  }


  invisible(NULL)
}
