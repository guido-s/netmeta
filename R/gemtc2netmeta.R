#' Create a netmeta object from a gemtc object
#' 
#' @description
#' Auxiliary function to create a netmeta object from a gemtc object
#' (experimental feature).
#' 
#' @param x A \code{gemtc} object created with R package \bold{gemtc}.
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
#' A netmeta object with additional class "netmeta.gemtc".
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}
#' 
#' @examples
#' \donttest{
#' if (requireNamespace("gemtc", quietly = TRUE)) {
#' library("gemtc")
#' # Load the example network and generate a consistency model:
#' model <- mtc.model(smoking, type="consistency")
#' results <- mtc.run(model, thin = 10)
#' detach("package:gemtc")
#'
#' # Print results in netmeta layout
#' mtc <- gemtc2netmeta(results)
#' mtc
#' # SUCRAs
#' netrank(mtc)
#' # Rankogram
#' rankogram(mtc)
#' }
#' }
#' 
#' @export gemtc2netmeta

gemtc2netmeta <- function(x,
                          keep.samples = TRUE,
                          level = gs("level"),
                          level.ma = gs("level.ma"),
                          reference.group =
                            levels(x$model$network$treatments$id)[1],
                          baseline.reference = gs("baseline.reference"),
                          small.values = gs("small.values"),
                          all.treatments = gs("all.treatments"),
                          seq = gs("seq"),
                          backtransf = gs("backtransf"),
                          nchar.trts = gs("nchar.trts"),
                          nchar.studlab = gs("nchar.studlab")) {
  
  #
  #
  # (1) Check for gemtc object and extract data / settings
  #
  #
  
  chkclass(x, "mtc.result")
  #
  if (!is.null(x$model$regressor))
    stop("R function gemtc2netmeta() not suitable for network meta-regression.",
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
  random <- x$model$linearModel == "random"
  common <- !random
  #
  likelihood <- x$model$likelihood
  link <- x$model$link
  #
  # Determine summary measure
  #
  assign_rules <- data.frame(
    lik = rep(c("normal", "binom", "poisson"), c(2, 3, 1)),
    lin = c("identity", "smd", "logit", "log", "cloglog", "log"),
    sm = c("MD", "SMD", "OR", "RR", "HR", "IRR")
  )
  #
  # Get rid of warning 'Undefined global functions or variables'
  #
  lik <- lin <- NULL
  #
  sm <- assign_rules %>% filter(lik == likelihood, lin == link) %>%
    select(sm) %>% pull()
  #
  if (sm == "HR")
    stop("Not yet implemented.")
  #
  trts <- levels(x$model$network$treatments$id)
  labels <- sort(trts)
  #
  trts.long <- x$model$network$treatments$description
  #
  if (!is.null(seq))
    seq <- setseq(seq, labels)
  else {
    seq <- labels
    if (is.numeric(seq))
      seq <- as.character(seq)
  }
  
  
  #
  #
  # (2) Calculate matrices with treatment effects, standard errors, etc.
  #
  #
  
  dmat <- do.call(rbind, x$samples) %>% data.frame() %>%
    select(starts_with("d."))
  dmat <- cbind(0, dmat)
  names(dmat) <- trts
  #
  tau <- NULL
  tmat <- do.call(rbind, x$samples) %>%
    data.frame() %>% select(starts_with("sd."))
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
  #
  # (3) Create pairwise object from gemtc data and run network meta-analysis
  #
  #
  
  #
  # Arm-level data
  #
  dat.gemtc <- x$model$network$data.ab
  #
  if (!is.null(dat.gemtc)) {
    #
    # Get rid of warning 'Undefined global functions or variables'
    #
    treatment <- study <- responders <- sampleSize <- exposure <-
      mean <- std.err <- std.dev <- NULL
    #
    if (sm %in% c("OR", "RR")) {
      suppressWarnings(
        dat <-
          pairwise(treat = treatment, studlab = study,
                   event = responders, n = sampleSize,
                   data = dat.gemtc,
                   sm = sm, warn = FALSE))
    }
    else if (sm == "MD") {
      #
      # Only means and standard errors provided
      #
      if (is.null(dat.gemtc$std.dev)) {
        suppressWarnings(
          dat <-
            pairwise(treat = treatment,
                     TE = mean, seTE = std.err,
                     studlab = study,
                     data = dat.gemtc,
                     sm = sm, warn = FALSE))
      }
      else {
        suppressWarnings(
          dat <-
            pairwise(treat = treatment,
                     n = sampleSize, mean = mean, sd = std.dev,
                     studlab = study,
                     data = dat.gemtc,
                     sm = sm, warn = FALSE))
      }
    }
    else if (sm == "SMD") {
      suppressWarnings(
        dat <-
          pairwise(treat = treatment,
                   n = sampleSize, mean = mean, sd = std.dev,
                   studlab = study,
                   data = dat.gemtc,
                   sm = sm, warn = FALSE))
    }
    else if (sm == "IRR") {
      suppressWarnings(
        dat <-
          pairwise(treat = treatment, studlab = study,
                   event = responders, time = exposure,
                   data = dat.gemtc,
                   sm = sm, warn = FALSE))
    }
  }
  else
    stop("R function gemtc2netmeta() not suitable for relative effect data.",
         call. = FALSE)
  #
  res <- suppressWarnings(
    netmeta(dat, common = common, random = random, warn = FALSE))
  
  
  #
  #
  # (4) Add gemtc results to netmeta object
  #
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
  res$method <- "gemtc"
  res$method.tau <- ""
  res$m <- NA
  #
  res$level <- gs("level")
  res$level.ma <- level.ma
  res$reference.group <- setchar(reference.group, trts)
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
  #
  # (5) Return netmeta object
  #
  #
  
  class(res) <- c(class(res), "netmeta.gemtc")
  #
  res
}
