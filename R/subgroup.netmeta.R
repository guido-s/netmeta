#' Subgroup analysis for network meta-analysis
#'
#' @description
#' Subgroup analysis for objects of class \code{netmeta}.
#'
#' @param x An object of class \code{netmeta} (or \code{subgroup.netmeta}).
#' @param subgroup A vector defining the subgroups considered in the
#'   network meta-analysis.
#' @param only.connected A logical indicating whether networks of subgroups
#'   must be connected.
#' @param common A logical indicating whether results for common
#'   effect subgroup network meta-analysis should be printed.
#' @param random A logical indicating whether results for random
#'   effects subgroup network meta-analysis should be printed.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau}. Either \code{"DL"}, \code{"REML"}, or
#'   \code{"ML"}, can be abbreviated.
#' @param level.ma The level used to calculate confidence intervals
#'   for network estimates.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names
#'   (see \code{\link{netmeta}}).
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard
#'   errors.
#' @param digits.Q Minimal number of significant digits for
#'   heterogeneity statistic Q, see \code{print.default}.
#' @param digits.pval.Q Minimal number of significant digits for
#'   p-value of heterogeneity test, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance \eqn{\tau^2}, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for
#'   \eqn{\tau}, the square root of the between-study variance
#'   \eqn{\tau^2}.
#' @param big.mark A character used as thousands separator.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
#' @param print.se A logical specifying whether standard errors should be
#'   printed.
#' @param print.tau2 A logical specifying whether between-study
#'   variance \eqn{\tau^2} should be printed.
#' @param print.tau A logical specifying whether \eqn{\tau}, the
#'   square root of the between-study variance \eqn{\tau^2}, should be
#'   printed.
#' @param print.Q A logical value indicating whether to print the
#'   results of the test of heterogeneity.
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param details.methods A logical specifying whether details on statistical
#'   methods should be printed.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param \dots Additional arguments.
#'
#' @return
#' An object of class \code{"subgroup.netmeta"} with corresponding \code{print}
#'   and \code{forest} function.
#'
#' @seealso \code{\link{forest.subgroup.netmeta}}, \code{\link{forest.netmeta}}
#'
#' @examples
#' \dontrun{
#' data("Senn2013")
#' # Add variable with (fictitious) risk of bias values
#' Senn2013$rob <- NA
#' set.seed(1909)
#' for (i in unique(Senn2013$studlab))
#'   Senn2013$rob[Senn2013$studlab == i] <- sample(1:3, 1)
#' Senn2013$rob <- factor(Senn2013$rob, levels = 1:3,
#'   labels = c("low", "moderate", "high"))
#' # Conduct network meta-analysis
#' net <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'   data = Senn2013, sm = "MD", reference = "plac", nchar.trts = 4)
#' # Conduct subgroup network meta-analysis
#' subgroup(net, rob, common = FALSE)
#' }
#'
#' @method subgroup netmeta
#' @export

subgroup.netmeta <- function(x, subgroup, only.connected = FALSE,
                             common = x$common, random = x$random,
                             method.tau = x$method.tau,
                             level.ma = x$level.ma,
                             backtransf = x$backtransf,
                             nchar.trts = x$nchar.trts,
                             ...) {
  
  chkclass(x, "netmeta")
  #
  chklogical(only.connected)
  chklogical(common)
  chklogical(random)
  chklevel(level.ma)
  chklogical(backtransf)
  #
  method.tau <- setchar(method.tau, c("DL", "REML", "ML"))
  #
  chknumeric(nchar.trts, min = 1, length = 1)
  
  
  # Get rid of warning 'Undefined global functions or variables'
  .studlab <- .treat1 <- .treat2 <- .TE <- .seTE <-
    comparison <- subnet <- NULL
    
  trts <- x$seq
  subgroup.name <- deparse(substitute(subgroup))
  #
  subgroup <- catch("subgroup", match.call(), x$data, sys.frame(sys.parent()))
  chknull(subgroup)
  if (is.factor(subgroup))
    levs <- levels(subgroup)
  else
    levs <- sort(unique(subgroup))
  #
  matNA <- x$TE.common[x$seq, x$seq]
  matNA[!is.na(matNA)] <- NA
  #
  TE.common <- seTE.common <- lower.common <- upper.common <-
    TE.random <- seTE.random <- lower.random <- upper.random <-
    tau2.matrix <- tau.matrix <- A.matrix <-
    array(NA, dim = c(nrow(matNA), ncol(matNA), length(levs)),
          dimnames = list(rownames(matNA), colnames(matNA), levs))
  #
  networks <- vector("list", length(levs))
  names(networks) <- levs
  #
  res.c <- res.r <- data.frame()
  #
  for (k in seq_along(levs)) {
    dat.k <- subset(x$data, subgroup == levs[k])
    #
    nc.k <- netconnection(.treat1, .treat2, .studlab, data = dat.k)
    #
    if (only.connected & nc.k$n.subnets > 1)
      stop("Subgroup '", subgroup.name, " = ", levs[k],
           "' does not have a connected network.",
           call. = FALSE)
    #
    dat.nc.k <- as.data.frame(nc.k)
    #
    for (s in seq_len(nc.k$n.subnets)) {
      studies.s <- subset(dat.nc.k, subnet == s)$studlab
      dat.ks <- subset(dat.k, .studlab %in% studies.s)
      #
      net.ks <- netmeta(.TE, .seTE, .treat1, .treat2, .studlab, data = dat.ks,
                        sm = x$sm, method.tau = method.tau,
                        level.ma = level.ma, backtransf = backtransf)
      if (nc.k$n.subnets > 1)
        networks[[k]][[s]] <- net.ks
      else
        networks[[k]] <- net.ks
      #
      for (i in net.ks$trts) {
        for (j in net.ks$trts) {
          TE.common[i, j, k] <- net.ks$TE.common[i, j]
          seTE.common[i, j, k] <- net.ks$seTE.common[i, j]
          lower.common[i, j, k] <- net.ks$lower.common[i, j]
          upper.common[i, j, k] <- net.ks$upper.common[i, j]
          #
          TE.random[i, j, k] <- net.ks$TE.random[i, j]
          seTE.random[i, j, k] <- net.ks$seTE.random[i, j]
          lower.random[i, j, k] <- net.ks$lower.random[i, j]
          upper.random[i, j, k] <- net.ks$upper.random[i, j]
          #
          tau2.matrix[i, j, k] <- net.ks$tau2
          tau.matrix[i, j, k] <- net.ks$tau
          #
          A.matrix[i, j, k] <- net.ks$A.matrix[i, j]
        }
      }
    }
    #
    diag(TE.common[, , k]) <- NA
    diag(seTE.common[, , k]) <- NA
    diag(lower.common[, , k]) <- NA
    diag(upper.common[, , k]) <- NA
    #
    diag(TE.random[, , k]) <- NA
    diag(seTE.random[, , k]) <- NA
    diag(lower.random[, , k]) <- NA
    diag(upper.random[, , k]) <- NA
    #
    diag(tau2.matrix[, , k]) <- NA
    diag(tau.matrix[, , k]) <- NA
    diag(A.matrix[, , k]) <- NA
  }
  #
  for (i in seq_len(length(trts) - 1)) {
    for (j in (i + 1):length(trts)) {
      m.ij <- metagen(TE.common[trts[i], trts[j], ],
                      seTE.common[trts[i], trts[j], ],
                      method.tau = "DL", method.tau.ci = "")
      #
      for (k in levs) {
        res.ijk <-
          data.frame(treat1 = trts[i], treat2 = trts[j],
                     subgroup = k,
                     k = A.matrix[trts[i], trts[j], k],
                     TE = TE.common[trts[i], trts[j], k],
                     seTE = seTE.common[trts[i], trts[j], k],
                     lower = lower.common[trts[i], trts[j], k],
                     upper = upper.common[trts[i], trts[j], k],
                     Q = m.ij$Q, df.Q = m.ij$df.Q, pval.Q = m.ij$pval.Q)
        #
        if (nrow(res.c) == 0)
          res.c <- res.ijk
        else
          res.c <- rbind(res.c, res.ijk)
      }
    }
  }
  #
  for (i in seq_len(length(trts) - 1)) {
    for (j in (i + 1):length(trts)) {
      m.ij <- metagen(TE.random[trts[i], trts[j], ],
                      seTE.random[trts[i], trts[j], ],
                      method.tau = method.tau, method.tau.ci = "")
      #
      for (k in levs) {
        res.ijk <-
          data.frame(treat1 = trts[i], treat2 = trts[j],
                     subgroup = k,
                     k = A.matrix[trts[i], trts[j], k],
                     TE = TE.random[trts[i], trts[j], k],
                     seTE = seTE.random[trts[i], trts[j], k],
                     lower = lower.random[trts[i], trts[j], k],
                     upper = upper.random[trts[i], trts[j], k],
                     tau2 = tau2.matrix[trts[i], trts[j], k],
                     tau = tau.matrix[trts[i], trts[j], k],
                     Q = m.ij$Q, df.Q = m.ij$df.Q, pval.Q = m.ij$pval.Q)
        #
        if (nrow(res.r) == 0)
          res.r <- res.ijk
        else
          res.r <- rbind(res.r, res.ijk)
      }
    }
  }
  #
  res.c <- subset(res.c, !is.na(res.c$TE))
  #
  sel.c <- duplicated(res.c[, c("Q", "df.Q", "pval.Q")]) & res.c$df.Q != 0
  res.c$Q[sel.c] <- NA
  res.c$df.Q[sel.c] <- NA
  res.c$pval.Q[sel.c] <- NA
  #
  res.r <- subset(res.r, !is.na(res.r$TE))
  #
  sel.r <- duplicated(res.r[, c("Q", "df.Q", "pval.Q")]) & res.r$df.Q != 0
  res.r$Q[sel.r] <- NA
  res.r$df.Q[sel.r] <- NA
  res.r$pval.Q[sel.r] <- NA
  #
  res.r <- subset(res.r, !is.na(res.r$TE))
  #
  rownames(res.c) <- seq_len(nrow(res.c))
  rownames(res.r) <- seq_len(nrow(res.r))
  #
  x$common <- common
  x$random <- random
  x$method.tau <- method.tau
  x$level.ma <- level.ma
  x$backtransf <- backtransf
  x$nchar.trts <- nchar.trts
  #
  res <- list(common = res.c, random = res.r,
              #
              TE.common = TE.common, seTE.common = seTE.common,
              lower.common = lower.common, upper.common = upper.common,
              #
              TE.random = TE.random, seTE.random = seTE.random,
              lower.random = lower.random, upper.random = upper.random,
              #
              tau2.matrix = tau2.matrix, tau.matrix = tau.matrix,
              A.matrix = A.matrix,
              networks = networks,
              #
              tau.preset = x$tau.preset,
              method.tau = x$method.tau,
              x = x,
              call = match.call(),
              version = packageDescription("netmeta")$Version)
  #
  class(res) <- "subgroup.netmeta"
  #
  res
}


#' @rdname subgroup.netmeta
#' @export subgroup

subgroup <- function(x, ...)
  UseMethod("subgroup")


#' @rdname subgroup.netmeta
#' @method print subgroup.netmeta
#' @export

print.subgroup.netmeta <- function(x,
                                   common = x$x$common,
                                   random = x$x$random,
                                   backtransf = x$x$backtransf,
                                   #
                                   nchar.trts = x$x$nchar.trts,
                                   #
                                   digits = gs("digits"),
                                   digits.se = gs("digits.se"),
                                   digits.Q = gs("digits.Q"),
                                   digits.pval.Q = gs("digits.pval.Q"),
                                   digits.tau2 = gs("digits.tau2"),
                                   digits.tau = gs("digits.tau"),
                                   #
                                   big.mark = gs("big.mark"),
                                   scientific.pval = gs("scientific.pval"),
                                   zero.pval = gs("zero.pval"),
                                   JAMA.pval = gs("JAMA.pval"),
                                   #
                                   print.se = !backtransf,
                                   print.tau2 = gs("print.tau2"),
                                   print.tau = gs("print.tau"),
                                   print.Q = gs("print.Q"),
                                   #
                                   text.tau2 = gs("text.tau2"),
                                   text.tau = gs("text.tau"),
                                   #
                                   details.methods = gs("details"),
                                   legend = gs("legend"),
                                   #
                                   ...) {
  
  chkclass(x, "subgroup.netmeta")
  #
  chklogical(common)
  chklogical(random)
  chklogical(backtransf)
  #
  if (is.null(nchar.trts))
    nchar.trts <- 666
  else
    chknumeric(nchar.trts, min = 1, length = 1)
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.pval.Q, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  #
  chkchar(big.mark, length = 1)
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  #
  chklogical(print.se)
  chklogical(print.tau2)
  chklogical(print.tau)
  chklogical(print.Q)
  #
  chkchar(text.tau2)
  chkchar(text.tau)
  #
  chklogical(details.methods)
  chklogical(legend)
  
  # Get rid of warning 'Undefined global functions or variables'
  .seTE <- .studlab <- .TE <- .treat1 <- .treat2 <- comparison <-
    df.Q <- subnet <- treat1 <- treat2 <- k <- TE <- NULL
  #
  trts <- x$x$trts
  trts.abbr <- treats(x$x$trts, nchar.trts)
  sm <- sm.lab <- x$x$sm
  #
  if (!backtransf & (is_relative_effect(sm) | sm == "VE"))
    sm.lab <- paste0("log", if (sm == "VE") "VR" else sm)
  ##
  ci.lab <- paste0(round(100 * x$x$level.ma, 1), "%-CI")
  
  if (common) {
    cat("Common effects model: \n\n")
    #
    dat.common <- x$common %>% filter(is.na(df.Q) | df.Q > 0)
    #
    if (backtransf) {
      dat.common$TE    <- backtransf(dat.common$TE, sm)
      dat.common$lower <- backtransf(dat.common$lower, sm)
      dat.common$upper <- backtransf(dat.common$upper, sm)
      #
      # Switch lower and upper limit for VE if results have been
      # backtransformed
      #
      if (sm == "VE") {
        tmp.l <- dat.common$lower
        dat.common$lower <- dat.common$upper
        dat.common$upper <- tmp.l
      }
    }
    #
    dat.common$TE <-
      formatN(dat.common$TE, digits = digits, text.NA = ".")
    #
    if (print.se)
      dat.common$seTE <-
        formatN(dat.common$seTE, digits = digits.se, text.NA = ".")
    else
      dat.common$seTE <- NULL
    #
    dat.common$ci <-
      formatCI(formatN(round(dat.common$lower, digits),
                       digits, ".", big.mark = big.mark),
               formatN(round(dat.common$upper, digits),
                       digits, ".", big.mark = big.mark))
    #
    if (print.Q) {
      dat.common$Q <-
        formatN(dat.common$Q, digits = digits.Q, text.NA = ".")
      dat.common$df.Q <-
        formatN(dat.common$df.Q, digits = 0, text.NA = ".")
      dat.common$pval.Q <-
        formatPT(dat.common$pval.Q,
                 digits = digits.pval.Q, lab.NA = ".",
                 scientific = scientific.pval,
                 zero = zero.pval, JAMA = JAMA.pval)
    }
    #
    if (any(trts != trts.abbr)) {
      dat.common$treat1 <-
        factor(dat.common$treat1, levels = trts, labels = trts.abbr)
      dat.common$treat2 <-
        factor(dat.common$treat2, levels = trts, labels = trts.abbr)
    }
    #
    rowlabs <- paste(dat.common$treat1, dat.common$treat2, sep = x$x$sep.trts)
    dat.common %<>% select(-treat1, -treat2) %>%
      select(subgroup, k, TE, if (print.se) "seTE", ci,
             if (print.Q) "Q", if (print.Q) "df.Q", if (print.Q) "pval.Q")
    #
    names(dat.common)[names(dat.common) == "TE"] <- sm.lab
    names(dat.common)[names(dat.common) == "ci"] <- ci.lab
    #
    prmatrix(dat.common, quote = FALSE, right = TRUE, rowlab = rowlabs)
    #
    if (random)
      cat("\n")
  }
  #
  # Print results for random effects model
  #
  if (random) {
    cat("Random effects model: \n\n")
    #
    dat.random <- x$random %>% filter(is.na(df.Q) | df.Q > 0)
    #
    if (backtransf) {
      dat.random$TE    <- backtransf(dat.random$TE, sm)
      dat.random$lower <- backtransf(dat.random$lower, sm)
      dat.random$upper <- backtransf(dat.random$upper, sm)
      #
      # Switch lower and upper limit for VE if results have been
      # backtransformed
      #
      if (sm == "VE") {
        tmp.l <- dat.random$lower
        dat.random$lower <- dat.random$upper
        dat.random$upper <- tmp.l
      }
    }
    #
    dat.random$TE <-
      formatN(dat.random$TE, digits = digits, text.NA = ".",
              big.mark = big.mark)
    #
    if (print.se)
      dat.random$seTE <-
        formatN(dat.random$seTE, digits = digits.se, text.NA = ".")
    else
      dat.random$seTE <- NULL
    #
    dat.random$ci <-
      formatCI(formatN(round(dat.random$lower, digits),
                     digits, ".", big.mark = big.mark),
             formatN(round(dat.random$upper, digits),
                     digits, ".", big.mark = big.mark))
    #
    if (print.Q) {
      dat.random$Q <-
        formatN(dat.random$Q, digits = digits.Q, text.NA = ".",
                big.mark = big.mark)
      dat.random$df.Q <-
        formatN(dat.random$df.Q, digits = 0, text.NA = ".",
                big.mark = big.mark)
      dat.random$pval.Q <-
        formatPT(dat.random$pval.Q,
                 digits = digits.pval.Q, lab.NA = ".",
                 scientific = scientific.pval,
                 zero = zero.pval, JAMA = JAMA.pval)
    }
    #
    if (print.tau2)
      dat.random$tau2 <-
      formatPT(dat.random$tau2, digits = digits.tau2, lab.NA = ".",
               big.mark = big.mark)
    if (print.tau2)
      dat.random$tau <-
      formatPT(dat.random$tau, digits = digits.tau, lab.NA = ".",
               big.mark = big.mark)
    #
    if (any(trts != trts.abbr)) {
      dat.random$treat1 <-
        factor(dat.random$treat1, levels = trts, labels = trts.abbr)
      dat.random$treat2 <-
        factor(dat.random$treat2, levels = trts, labels = trts.abbr)
    }
    #
    rowlabs <- paste(dat.random$treat1, dat.random$treat2, sep = x$x$sep.trts)
    dat.random %<>% select(-treat1, -treat2) %>%
      select(subgroup, k, TE, if (print.se) "seTE", ci,
             if (print.tau2) "tau2", if (print.tau) "tau",
             if (print.Q) "Q", if (print.Q) "df.Q", if (print.Q) "pval.Q")
    #
    names(dat.random)[names(dat.random) == "TE"] <- sm.lab
    names(dat.random)[names(dat.random) == "ci"] <- ci.lab
    #
    prmatrix(dat.random, quote = FALSE, right = TRUE, rowlab = rowlabs)
  }
  #
  # Print details of network meta-analysis methods
  #
  if (details.methods) {
    text.details <-
      textmeth(x, random, print.tau2, print.tau,
               text.tau2, text.tau, digits.tau2, digits.tau,
               big.mark = big.mark)
    #
    cat(text.details)
  }
  #
  # Add legend with abbreviated treatment labels
  #
  legendabbr(trts, trts.abbr, legend)
  #
  invisible(NULL)
}
