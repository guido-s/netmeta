#' Create a netmeta object from a crossnma object
#' 
#' @description
#' Auxiliary function to create a netmeta object from a crossnma object
#' 
#' @param x A \code{crossnma} object created with R package \bold{crossnma}.
#' @param keep.samples A logical indicating whether to keep the generated
#'   samples.
#' @param level The level used to calculate confidence intervals for
#'   individual comparisons.
#' @param level.ma The level used to calculate confidence intervals
#'   for network estimates.
#' @param reference.group Reference treatment (first treatment is used
#'   if argument is missing).
#' @param baseline.reference A logical indicating whether results
#'   should be expressed as comparisons of other treatments versus the
#'   reference treatment (default) or vice versa. This argument is
#'   only considered if \code{reference.group} has been specified.
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"desirable"}) or
#'   harmful (\code{"undesirable"}) effect (passed on to
#'   \code{\link{netrank}}, can be abbreviated.
#' @param all.treatments A logical or \code{"NULL"}. If \code{TRUE},
#'   matrices with all treatment effects, and confidence limits will
#'   be printed.
#' @param seq A character or numerical vector specifying the sequence
#'   of treatments in printouts.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
#' @param nchar.studlab A numeric defining the minimum number of
#'   characters used to create unique study labels.
#'
#' @return
#' A netmeta object with additional class "netmeta.crossnma".
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}
#' 
#' @examples
#' \dontrun{
#' if (requireNamespace("crossnma", quietly = TRUE)) {
#' library("crossnma")
#' # Create a JAGS model
#' mod <- crossnma.model(treat, id, relapse, n, design,
#'   prt.data = ipddata, std.data = stddata,
#'   reference = "A", trt.effect = "random", method.bias = "naive")
#' # Fit JAGS model
#' set.seed(1909)
#' fit <- crossnma(mod)
#'
#' # Print results in netmeta layout
#' cro <- crossnma2netmeta(fit)
#' cro
#' # SUCRAs
#' netrank(cro)
#' # Rankogram
#' rankogram(cro)
#' }
#' }
#' 
#' @export crossnma2netmeta

crossnma2netmeta <- function(x,
                             keep.samples = TRUE,
                             level = gs("level"), level.ma = x$model$level.ma,
                             reference.group = x$model$reference,
                             baseline.reference = gs("baseline.reference"),
                             small.values = x$model$small.values,
                             all.treatments = gs("all.treatments"),
                             seq = gs("seq"),
                             backtransf = x$model$backtransf,
                             nchar.trts = gs("nchar.trts"),
                             nchar.studlab = gs("nchar.studlab")) {
  
  #
  # (1) Check for crossnma object and extract data / settings
  #
  
  chkclass(x, "crossnma")
  #
  if (!is.null(x$model$covariate))
    stop("R function crossnma2netmeta() not suitable for ",
         "network meta-regression.",
         call. = FALSE)
  #
  chklogical(keep.samples)
  chklevel(level)
  chklevel(level.ma)
  chklogical(baseline.reference)
  small.values <- setsv(replaceNULL(small.values, gs("small.values")))
  if (!is.null(all.treatments))
    chklogical(all.treatments)
  chklogical(backtransf)
  chknumeric(nchar.trts, min = 1, length = 1)
  chknumeric(nchar.studlab, min = 1, length = 1)
  #
  common <- x$model$trt.effect == "common"
  random <- x$model$trt.effect == "random"
  sm <- x$model$sm
  #
  trts <- as.character(x$trt.key$trt.ini)
  labels <- sort(trts)
  #
  if (!is.null(seq))
    seq <- setseq(seq, labels)
  else {
    seq <- labels
    if (is.numeric(seq))
      seq <- as.character(seq)
  }
  #
  dat <- x$model$all.data.ad
  trt <- dat$trt
  n <- dat$n
  outcome <- dat$outcome
  study <- dat$study
  se <- dat$se
  
  
  #
  # (2) Calculate matrices with treatment effects, standard errors, etc.
  #
  
  dmat <- do.call(rbind, x$samples) %>% data.frame() %>%
    select(starts_with("d."))
  names(dmat) <- trts
  #
  if (random) {
    tau <- NULL
    tmat <- do.call(rbind, x$samples) %>% data.frame() %>% select(tau)
  }
  #
  TE.matrix <- seTE.matrix <-
    lower.matrix <- upper.matrix <-
    matrix(NA, nrow = ncol(dmat), ncol = ncol(dmat),
           dimnames = list(trts, trts))
  #
  for (i in seq_len(ncol(dmat)))
    for (j in seq_len(ncol(dmat)))
      if (i != j)
        TE.matrix[i, j] <- mean(dmat[, i] - dmat[, j])
  #
  for (i in seq_len(ncol(dmat)))
    for (j in seq_len(ncol(dmat)))
      if (i != j)
        seTE.matrix[i, j] <- sd(dmat[, i] - dmat[, j])
  #
  for (i in seq_len(ncol(dmat)))
    for (j in seq_len(ncol(dmat)))
      if (i != j)
        lower.matrix[i, j] <- quantile(dmat[, i] - dmat[, j],
                                       (1 - level.ma) / 2)
  #
  for (i in seq_len(ncol(dmat)))
    for (j in seq_len(ncol(dmat)))
      if (i != j)
        upper.matrix[i, j] <- quantile(dmat[, i] - dmat[, j],
                                       1 - (1 - level.ma) / 2)
  #
  diag(TE.matrix) <- diag(seTE.matrix) <-
    diag(lower.matrix) <- diag(upper.matrix) <- 0
  
  
  #
  # (3) Create pairwise object from crossnma data and run network meta-analysis
  #
  
  if (sm %in% c("OR", "RR")) {
    suppressWarnings(
      dat <-
        pairwise(treat = trt, event = outcome, n = n,
                 studlab = study, sm = sm, warn = FALSE)
    )
  }
  else if (sm %in% c("MD", "SMD")) {
    suppressWarnings(
      dat <-
        pairwise(treat = trt,
                 mean = outcome, n = n, sd = se,
                 studlab = study, sm = sm, warn = FALSE)
    )
  }
  #
  res <- suppressWarnings(netmeta(dat, common = common, random = random,
                                  warn = FALSE))
  
  
  #
  # (4) Add crossnma results to netmeta object
  #
  
  res <- setNA_nma(res)
  #
  res$keep.samples <- keep.samples
  #
  if (keep.samples) {
    res$samples <- list(d = dmat)
    #
    if (random)
      res$samples$tau <- tmat
  }
  #
  res$method <- "crossnma"
  res$method.tau <- ""
  res$m <- NA
  #
  res$level <- gs("level")
  res$level.ma <- level.ma
  res$reference.group <- reference.group
  res$baseline.reference <- baseline.reference
  res$small.values <- small.values
  res$all.treatments <- all.treatments
  res$seq <- seq
  res$backtransf <- backtransf
  res$nchar.trts <- nchar.trts
  res$nchar.studlab <- nchar.studlab
  #
  if (common) {
    res$TE.common <- TE.matrix[labels, labels]
    res$seTE.common <- seTE.matrix[labels, labels]
    res$lower.common <- lower.matrix[labels, labels]
    res$upper.common <- upper.matrix[labels, labels]
    #
    for (i in seq_along(res$treat1)) {
      res$TE.nma.common[i] <- res$TE.common[res$treat1[i], res$treat2[i]]
      res$seTE.nma.common[i] <- res$seTE.common[res$treat1[i], res$treat2[i]]
      res$lower.nma.common[i] <- res$lower.common[res$treat1[i], res$treat2[i]]
      res$upper.nma.common[i] <- res$upper.common[res$treat1[i], res$treat2[i]]
    }
    #
    res$tau2 <- res$tau <- NA
    #
    res$TE.random[!is.na(res$TE.random)] <- NA
    res$seTE.random <- res$lower.random <- res$upper.random <- res$TE.random
  }
  else if (random) {
    res$TE.random <- TE.matrix[labels, labels]
    res$seTE.random <- seTE.matrix[labels, labels]
    res$lower.random <- lower.matrix[labels, labels]
    res$upper.random <- upper.matrix[labels, labels]
    #
    for (i in seq_along(res$treat1)) {
      res$TE.nma.random[i] <- res$TE.random[res$treat1[i], res$treat2[i]]
      res$seTE.nma.random[i] <- res$seTE.random[res$treat1[i], res$treat2[i]]
      res$lower.nma.random[i] <- res$lower.random[res$treat1[i], res$treat2[i]]
      res$upper.nma.random[i] <- res$upper.random[res$treat1[i], res$treat2[i]]
    }
    #
    res$tau2 <- mean(tmat %>% unlist() %>% unname())
    res$tau <- sqrt(res$tau2)
    #
    res$TE.common[!is.na(res$TE.common)] <- NA
    res$seTE.common <- res$lower.common <- res$upper.common <- res$TE.common
  }
  
  
  #
  # (5) Return netmeta object
  #
  
  class(res) <- c("netmeta.crossnma", class(res))
  #
  res
}
