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
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
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
#' @seealso \code{\link{forest.subgroup.netmeta}}, \code{\link{netmetareg}},
#'   \code{\link{forest.netmeta}}
#'
#' @examples
#' # Add examples
#'
#' @method subgroup netmeta
#' @export


subgroup.netmeta <- function(x, subgroup, only.connected = FALSE,
                             common = x$common, random = x$random,
                             method.tau = x$method.tau,
                             nchar.trts = x$nchar.trts,
                             ...) {
  
  chkclass(x, "netmeta")
  #
  chklogical(only.connected)
  chklogical(common)
  chklogical(random)
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
  levs <- sort(unique(subgroup))
  #
  matNA <- x$TE.common[x$seq, x$seq]
  matNA[!is.na(matNA)] <- NA
  #
  TE.common <- seTE.common <- TE.random <- seTE.random <-
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
                        sm = x$sm, method.tau = x$method.tau)
      if (nc.k$n.subnets > 1)
        networks[[k]][[s]] <- net.ks
      else
        networks[[k]] <- net.ks
      #
      for (i in net.ks$trts) {
        for (j in net.ks$trts) {
          TE.common[i, j, k] <- net.ks$TE.common[i, j]
          seTE.common[i, j, k] <- net.ks$seTE.common[i, j]
          TE.random[i, j, k] <- net.ks$TE.random[i, j]
          seTE.random[i, j, k] <- net.ks$seTE.random[i, j]
          tau2.matrix[i, j, k] <- net.ks$tau2
          tau.matrix[i, j, k] <- net.ks$tau
          A.matrix[i, j, k] <- net.ks$A.matrix[i, j]
        }
      }
    }
    #
    diag(TE.common[, , k]) <- NA
    diag(seTE.common[, , k]) <- NA
    diag(TE.random[, , k]) <- NA
    diag(seTE.random[, , k]) <- NA
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
                      method.tau = x$method.tau, method.tau.ci = "")
      #
      for (k in levs) {
        res.ijk <-
          data.frame(treat1 = trts[i], treat2 = trts[j],
                     subgroup = k,
                     k = A.matrix[trts[i], trts[j], k],
                     TE = TE.random[trts[i], trts[j], k],
                     seTE = seTE.random[trts[i], trts[j], k],
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
  x$nchar.trts <- nchar.trts
  #
  res <- list(common = res.c, random = res.r,
              TE.common = TE.common, seTE.common = seTE.common,
              TE.random = TE.random, seTE.random = seTE.random,
              tau2.matrix = tau2.matrix, tau.matrix = tau.matrix,
              A.matrix = A.matrix,
              networks = networks,
              x = x)
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
                                   scientific.pval = gs("scientific.pval"),
                                   big.mark = gs("big.mark"),
                                   #
                                   text.tau2 = gs("text.tau2"),
                                   text.tau = gs("text.tau"),
                                   #
                                   details.methods = gs("details.netmeta"),
                                   legend = gs("legend.netmeta"),
                                   #
                                   ...) {
  
  chkclass(x, "subgroup.netmeta")
  #
  chklogical(common)
  chklogical(random)
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
  chklogical(scientific.pval)
  #
  chkchar(text.tau2)
  chkchar(text.tau)
  #
  chklogical(details.methods)
  chklogical(legend)
  
  # Get rid of warning 'Undefined global functions or variables'
  .seTE <- .studlab <- .TE <- .treat1 <- .treat2 <- comparison <-
    df.Q <- subnet <- treat1 <- treat2 <- NULL
  #
  trts <- x$x$trts
  trts.abbr <- treats(x$x$trts, nchar.trts)
  
  if (common) {
    cat("Common effects model: \n\n")
    #
    common <- x$common %>% filter(is.na(df.Q) | df.Q > 0)
    #
    common$TE <- formatN(common$TE, digits = digits, text.NA = ".")
    common$seTE <- formatN(common$seTE, digits = digits.se, text.NA = ".")
    common$Q <- formatN(common$Q, digits = digits.Q, text.NA = ".")
    common$df.Q <- formatN(common$df.Q, digits = 0, text.NA = ".")
    common$pval.Q <- formatPT(common$pval.Q, digits = digits.pval.Q)
    #
    if (any(trts != trts.abbr)) {
      common$treat1 <- factor(common$treat1, levels = trts, labels = trts.abbr)
      common$treat2 <- factor(common$treat2, levels = trts, labels = trts.abbr)
    }
    #
    prmatrix(common %>% select(-treat1, -treat2), quote = FALSE, right = TRUE,
             rowlab = paste(common$treat1, common$treat2, sep = x$x$sep.trts))
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
    random <- x$random %>% filter(is.na(df.Q) | df.Q > 0)
    #
    random$TE <- formatN(random$TE, digits = digits, text.NA = ".",
                         big.mark = big.mark)
    random$seTE <- formatN(random$seTE, digits = digits.se, text.NA = ".",
                           big.mark = big.mark)
    random$Q <- formatN(random$Q, digits = digits.Q, text.NA = ".",
                        big.mark = big.mark)
    random$df.Q <- formatN(random$df.Q, digits = 0, text.NA = ".",
                           big.mark = big.mark)
    random$pval.Q <- formatPT(random$pval.Q, digits = digits.pval.Q,
                              scientific = scientific.pval)
    random$tau2 <- formatPT(random$tau2, digits = digits.tau2,
                            big.mark = big.mark)
    random$tau <- formatPT(random$tau, digits = digits.tau,
                           big.mark = big.mark)
    #
    if (any(trts != trts.abbr)) {
      random$treat1 <- factor(random$treat1, levels = trts, labels = trts.abbr)
      random$treat2 <- factor(random$treat2, levels = trts, labels = trts.abbr)
    }
    #
    prmatrix(random %>% select(-treat1, -treat2), quote = FALSE, right = TRUE,
             rowlab = paste(random$treat1,
                            random$treat2,
                            sep = x$x$sep.trts))
  }
  #
  # Print details of network meta-analysis methods
  #
  if (details.methods) {
    text.details <- catmeth(x, random, text.tau2, digits.tau2, big.mark)
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
