#' Conduct pairwise meta-analyses for all comparisons with direct
#' evidence in a network meta-analysis
#' 
#' @description
#' Conduct pairwise meta-analyses for all comparisons with direct
#' evidence in a network meta-analysis.
#' 
#' @param x An object of class \code{netmeta} or \code{netpairwise}.
#' @param object An object of class \code{netpairwise}.
#' @param separate A logical indicating whether results for pairwise
#'   comparisons should be printed as separate meta-analyses or as
#'   subgroups which is more concise.
#' @param common A logical indicating whether a common effects network
#'   meta-analysis should be conducted.
#' @param random A logical indicating whether a random effects network
#'   meta-analysis should be conducted.
#' @param level The level used to calculate confidence intervals for
#'   individual comparisons.
#' @param level.ma The level used to calculate confidence intervals
#'   for pooled estimates.
#' @param prediction A logical indicating whether prediction intervals
#'   should be printed.
#' @param level.predict The level used to calculate prediction
#'   intervals for a new study.
#' @param reference.group Reference treatment.
#' @param baseline.reference A logical indicating whether results
#'   should be expressed as comparisons of other treatments versus the
#'   reference treatment (default) or vice versa. This argument is
#'   only considered if \code{reference.group} has been specified.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau}. Either \code{"DL"}, \code{"REML"}, or
#'   \code{"ML"}, can be abbreviated.
#' @param method A character string indicating which method is to be
#'   used for pooling of studies, see \code{\link[meta]{metabin}}.
#' @param incr A numerical value which is added to cell counts,
#'   see \code{\link[meta]{metabin}}.
#' @param method.incr A character string indicating which continuity
#'   correction method should be used (\code{"only0"},
#'   \code{"if0all"}, or \code{"all"}), see \code{\link[meta]{metabin}}.
#' @param allstudies A logical indicating whether studies with zero
#'   events or non-events in all treatment arms should be included in
#'   the meta-analysis, see \code{\link[meta]{metabin}}.
#' @param order An optional character or numerical vector specifying
#'   the order of treatments.
#' @param sep.trts A character used in comparison names as separator
#'   between treatment labels.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param k.min Minimum number of studies in pairwise comparison to
#'   show funnel plot, radial plot or conduct test for funnel plot
#'   asymmetry.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (passed on to \code{metagen} or
#'   print functions and to catch deprecated arguments).
#' 
#' @details
#' Conduct pairwise meta-analyses for all comparisons with direct
#' evidence in a network meta-analysis. In contrast to
#' \code{\link{netmeta}} and \code{\link{netsplit}}, unadjusted
#' standard errors are used in the calculations and the between-study
#' heterogeneity variance is allowed to differ between comparisons.
#' 
#' The R function \code{\link[meta]{metagen}} is called internally.
#' 
#' @note
#' This function must not be confused with \code{\link[meta]{pairwise}}
#' which can be used as a pre-processing step to convert data from
#' arm-based to contrast-based format by calculating all pairwise
#' comparisons within a study.
#' 
#' @return
#' Either a single \code{\link[meta]{metagen}} object with pairwise
#' comparisons as subgroups or a list with \code{\link[meta]{metagen}}
#' objects for each direct pairwise comparison.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netsplit}},
#'   \code{\link[meta]{pairwise}}
#' 
#' @examples
#' oldsets <- settings.meta(digits = 2, digits.tau2 = 2, digits.tau = 2)
#' 
#' data(smokingcessation)
#' 
#' # Transform data from arm-based format to contrast-based format
#' #
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' 
#' # Conduct random effects network meta-analysis
#' #
#' net1 <- netmeta(p1, common = FALSE)
#' 
#' # Calculate and print concise results for all pairwise
#' # meta-analyses
#' #
#' np1 <- netpairwise(net1)
#' np1
#' print(np1, details.method = FALSE)
#'
#' \dontrun{
#' data(Senn2013)
#' 
#' # Random effects model
#' #
#' net2 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'   data = Senn2013, sm = "MD", common = FALSE, reference = "plac")
#' 
#' # Calculate and print concise results for all pairwise
#' # meta-analyses
#' #
#' np2 <- netpairwise(net2)
#' np2
#' print(np2, details.method = FALSE)
#'
#' forest(np2)
#' 
#' # Print detailed information for each pairwise comparison
#' #
#' np3 <- netpairwise(net2, separate = TRUE)
#' forest(np3)
#' funnel(np3)
#' radial(np3)
#' funnel(np3, k.min = 1)
#' }
#'
#' settings.meta(oldsets)
#' 
#' @rdname netpairwise
#' @export netpairwise

netpairwise <- function(x, ...)
  UseMethod("netpairwise")


#' @rdname netpairwise
#' @method netpairwise netmeta
#' @export

netpairwise.netmeta <- function(x,
                                separate = FALSE,
                                common = x$common,
                                random = x$random,
                                level = x$level,
                                level.ma = x$level.ma,
                                prediction = x$prediction,
                                level.predict = x$level.predict,
                                reference.group =
                                  if (missing(order)) x$reference.group else "",
                                baseline.reference = x$baseline.reference,
                                method.tau = x$method.tau,
                                order = NULL,
                                sep.trts = x$sep.trts,
                                nchar.trts = x$nchar.trts,
                                backtransf = x$backtransf,
                                warn.deprecated = gs("warn.deprecated"),
                                ...) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  chkclass(x, "netmeta")
  x <- updateversion(x)
  ##
  chklogical(separate)
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  common <- deprecated(common, missing(common), args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  chklogical(random)
  chklogical(prediction)
  ##
  chklevel(level)
  chklevel(level.ma)
  chklevel(level.predict)
  ##
  reference.group <- setref(reference.group, c(x$trts, ""))
  chklogical(baseline.reference)
  ##
  method.tau <- setchar(method.tau, c("DL", "ML", "REML"))
  ##
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  if (!missing(order)) {
    order <- catch("order", mc, x, sfsp)
    ##
    if (length(order) != length(x$trts))
      order <- setchar(order, x$trts)
    else
      order <- setseq(order, x$trts)
  }
  ##
  if (is.null(order))
    order <- x$trts
  ##
  chkchar(sep.trts)
  chknumeric(nchar.trts, min = 1, length = 1)
  chklogical(backtransf)
  
  
  ##
  ##
  ## (2) Get data
  ##
  ##
  TE <- x$data$.TE
  seTE <- x$data$.seTE
  studlab <- x$data$.studlab
  ##
  n1 <- x$data$.n1
  n2 <- x$data$.n2
  ##
  treat1 <- x$data$.treat1
  treat2 <- x$data$.treat2
  ##
  comparison <- paste(treat1, treat2, sep = sep.trts)
  comparison21 <- paste(treat2, treat1, sep = sep.trts)
  ##
  treat1.pos <- as.numeric(factor(treat1, levels = order))
  treat2.pos <- as.numeric(factor(treat2, levels = order))
  ##
  trts.abbr <- treats(x$trts, nchar.trts)
  trt1 <- as.character(factor(treat1, levels = x$trts, labels = trts.abbr))
  trt2 <- as.character(factor(treat2, levels = x$trts, labels = trts.abbr))
  ##
  comp <- paste(trt1, trt2, sep = sep.trts)
  comp21 <- paste(trt2, trt1, sep = sep.trts)
  ##
  wo <- treat1.pos > treat2.pos
  ##
  if (any(wo)) {
    ttreat1.pos <- treat1.pos
    treat1.pos[wo] <- treat2.pos[wo]
    treat2.pos[wo] <- ttreat1.pos[wo]
    ##
    TE[wo] <- -TE[wo]
    ##
    tn1 <- n1
    n1[wo] <- n2[wo]
    n2[wo] <- tn1[wo]
    ##
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
    ##
    comparison[wo] <- comparison21[wo]
    ##
    ttrt1 <- trt1
    trt1[wo] <- trt2[wo]
    trt2[wo] <- ttrt1[wo]
    ##
    comp[wo] <- comp21[wo]
  }
  ##
  o <- order(treat1.pos, treat2.pos)
  ##
  TE <- TE[o]
  seTE <- seTE[o]
  studlab <- studlab[o]
  ##
  n1 <- n1[o]
  n2 <- n2[o]
  ##
  treat1 <- treat1[o]
  treat2 <- treat2[o]
  comparison <- comparison[o]
  ##
  trt1 <- trt1[o]
  trt2 <- trt2[o]
  comp <- comp[o]
  ##
  if (reference.group != "") {
    if (baseline.reference) {
      wo1 <- trt1 == reference.group
      if (any(wo1)) {
        TE[wo1] <- -TE[wo1]
        ttrt1 <- trt1
        trt1[wo1] <- trt2[wo1]
        trt2[wo1] <- ttrt1[wo1]
        tn1 <- n1
        n1[wo1] <- n2[wo1]
        n2[wo1] <- tn1[wo1]
      }
    }
    else {
      wo2 <- trt2 == reference.group
      if (any(wo2)) {
        TE[wo2] <- -TE[wo2]
        ttrt1 <- trt1
        trt1[wo2] <- trt2[wo2]
        trt2[wo2] <- ttrt1[wo2]
        tn1 <- n1
        n1[wo2] <- n2[wo2]
        n2[wo2] <- tn1[wo2]
      }
    }
  }
  
  
  ##
  ##
  ## (3) Run pairwise meta-analyses
  ##
  ##
  if (!separate) {
    res <- metagen(TE, seTE, studlab = studlab,
                   n.e = n1, n.c = n2,
                   sm = x$sm,
                   subgroup = paste0(trt1, sep.trts, trt2),
                   subgroup.name = "comparison",
                   print.subgroup.name = FALSE,
                   common = common,
                   random = random,
                   level = level,
                   level.ma = level.ma,
                   prediction = prediction,
                   level.predict = level.predict,
                   method.tau = method.tau,
                   overall = FALSE, overall.hetstat = FALSE,
                   test.subgroup = FALSE,
                   warn.deprecated = FALSE,
                   ...)
    ##
    res$k.study <- x$k
    res$k <- x$m
    res$w.common[!is.na(res$w.common)] <- NA
    res$w.random[!is.na(res$w.random)] <- NA
    ##
    res$order <- order
    ##
    class(res) <- c(class(res), "netpairwise")
  }
  else {
    comps <- unique(data.frame(trt1, trt2))
    comps <- comps[order(comps$trt1, comps$trt2), ]
    n.comps <- nrow(comps)
    ##
    res <- vector("list", length = n.comps)
    ##
    for (i in seq_len(n.comps)) {
      comp.i <- paste0(comps$trt1[i], sep.trts, comps$trt2[i])
      res[[i]] <-
        metagen(TE, seTE, studlab = studlab,
                n.e = n1, n.c = n2,
                sm = x$sm,
                subset = trt1 == comps$trt1[i] & trt2 == comps$trt2[i],
                complab = comp.i,
                common = common,
                random = random,
                level = level,
                level.ma = level.ma,
                prediction = prediction,
                level.predict = level.predict,
                method.tau = method.tau,
                label.e = comps$trt1[i], label.c = comps$trt2[i],
                warn.deprecated = FALSE,
                ...)
    }
    ##
    attr(res, "order") <- order
    attr(res, "version") <- packageDescription("netmeta")$Version
    ##
    class(res) <- "netpairwise"
  }
  
  res
}


#' @rdname netpairwise
#' @method netpairwise netmetabin
#' @export

netpairwise.netmetabin <- function(x,
                                   separate = FALSE,
                                   common = x$common,
                                   random = x$random,
                                   level = x$level,
                                   level.ma = x$level.ma,
                                   prediction = x$prediction,
                                   level.predict = x$level.predict,
                                   reference.group =
                                     if (missing(order))
                                       x$reference.group else "",
                                   baseline.reference = x$baseline.reference,
                                   #
                                   method = x$method,
                                   incr = x$incr,
                                   method.incr = x$method.incr,
                                   allstudies = x$allstudies,
                                   #
                                   method.tau = x$method.tau,
                                   #
                                   order = NULL,
                                   sep.trts = x$sep.trts,
                                   nchar.trts = x$nchar.trts,
                                   backtransf = x$backtransf,
                                   warn.deprecated = gs("warn.deprecated"),
                                   ...) {
  
  
  #
  #
  # (1) Check and set arguments
  #
  #
  
  chkclass(x, "netmetabin")
  x <- updateversion(x)
  #
  chklogical(separate)
  #
  args  <- list(...)
  chklogical(warn.deprecated)
  common <- deprecated(common, missing(common), args, "fixed",
                       warn.deprecated)
  chklogical(common)
  #
  chklogical(random)
  chklogical(prediction)
  #
  chklevel(level)
  chklevel(level.ma)
  chklevel(level.predict)
  #
  reference.group <- setref(reference.group, c(x$trts, ""))
  chklogical(baseline.reference)
  #
  method <- setchar(method, c("MH", "NCH", "LRP"))
  if (method == "NCH")
    method <- "MH"
  #
  method.tau <- setchar(method.tau, c("DL", "ML", "REML"))
  #
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  if (!missing(order)) {
    order <- catch("order", mc, x, sfsp)
    #
    if (length(order) != length(x$trts))
      order <- setchar(order, x$trts)
    else
      order <- setseq(order, x$trts)
  }
  #
  if (is.null(order))
    order <- x$trts
  #
  chkchar(sep.trts)
  chknumeric(nchar.trts, min = 1, length = 1)
  chklogical(backtransf)
  
  
  #
  #
  # (2) Get data
  #
  #
  
  studlab <- x$data$.studlab
  #
  event1 <- x$data$.event1
  event2 <- x$data$.event2
  #
  n1 <- x$data$.n1
  n2 <- x$data$.n2
  #
  treat1 <- x$data$.treat1
  treat2 <- x$data$.treat2
  #
  comparison <- paste(treat1, treat2, sep = sep.trts)
  comparison21 <- paste(treat2, treat1, sep = sep.trts)
  #
  treat1.pos <- as.numeric(factor(treat1, levels = order))
  treat2.pos <- as.numeric(factor(treat2, levels = order))
  #
  trts.abbr <- treats(x$trts, nchar.trts)
  trt1 <- as.character(factor(treat1, levels = x$trts, labels = trts.abbr))
  trt2 <- as.character(factor(treat2, levels = x$trts, labels = trts.abbr))
  #
  comp <- paste(trt1, trt2, sep = sep.trts)
  comp21 <- paste(trt2, trt1, sep = sep.trts)
  #
  wo <- treat1.pos > treat2.pos
  #
  if (any(wo)) {
    ttreat1.pos <- treat1.pos
    treat1.pos[wo] <- treat2.pos[wo]
    treat2.pos[wo] <- ttreat1.pos[wo]
    #
    tevent1 <- event1
    event1[wo] <- event2[wo]
    event2[wo] <- tevent1[wo]
    #
    tn1 <- n1
    n1[wo] <- n2[wo]
    n2[wo] <- tn1[wo]
    #
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
    #
    comparison[wo] <- comparison21[wo]
    #
    ttrt1 <- trt1
    trt1[wo] <- trt2[wo]
    trt2[wo] <- ttrt1[wo]
    #
    comp[wo] <- comp21[wo]
  }
  #
  o <- order(treat1.pos, treat2.pos)
  #
  studlab <- studlab[o]
  #
  event1 <- event1[o]
  event2 <- event2[o]
  #
  n1 <- n1[o]
  n2 <- n2[o]
  #
  treat1 <- treat1[o]
  treat2 <- treat2[o]
  comparison <- comparison[o]
  #
  trt1 <- trt1[o]
  trt2 <- trt2[o]
  comp <- comp[o]
  #
  if (reference.group != "") {
    if (baseline.reference) {
      wo1 <- trt1 == reference.group
      if (any(wo1)) {
        ttrt1 <- trt1
        trt1[wo1] <- trt2[wo1]
        trt2[wo1] <- ttrt1[wo1]
        #
        tevent1 <- event1
        event1[wo1] <- event2[wo1]
        event2[wo1] <- tevent1[wo1]
        #
        tn1 <- n1
        n1[wo1] <- n2[wo1]
        n2[wo1] <- tn1[wo1]
      }
    }
    else {
      wo2 <- trt2 == reference.group
      if (any(wo2)) {
        ttrt1 <- trt1
        trt1[wo2] <- trt2[wo2]
        trt2[wo2] <- ttrt1[wo2]
        #
        tevent1 <- event1
        event1[wo2] <- event2[wo2]
        event2[wo2] <- tevent1[wo2]
        #
        tn1 <- n1
        n1[wo2] <- n2[wo2]
        n2[wo2] <- tn1[wo2]
      }
    }
  }
  
  
  #
  #
  # (3) Run pairwise meta-analyses
  #
  #
  if (!separate) {
    res <- metabin(event1, n1, event2, n2, studlab = studlab,
                   method = method, sm = x$sm,
                   incr = incr, method.incr = method.incr,
                   allstudies = allstudies,
                   subgroup = paste0(trt1, sep.trts, trt2),
                   subgroup.name = "comparison",
                   print.subgroup.name = FALSE,
                   common = common,
                   random = random,
                   level = level,
                   level.ma = level.ma,
                   prediction = prediction,
                   level.predict = level.predict,
                   method.tau = method.tau,
                   overall = FALSE, overall.hetstat = FALSE,
                   test.subgroup = FALSE,
                   warn = FALSE,
                   warn.deprecated = FALSE,
                   ...)
    #
    res$k.study <- x$k
    res$k <- x$m
    res$w.common[!is.na(res$w.common)] <- NA
    res$w.random[!is.na(res$w.random)] <- NA
    #
    res$order <- order
    #
    class(res) <- c(class(res), "netpairwise")
  }
  else {
    comps <- unique(data.frame(trt1, trt2))
    comps <- comps[order(comps$trt1, comps$trt2), ]
    n.comps <- nrow(comps)
    #
    res <- vector("list", length = n.comps)
    #
    for (i in seq_len(n.comps)) {
      comp.i <- paste0(comps$trt1[i], sep.trts, comps$trt2[i])
      res[[i]] <-
        metabin(event1, n1, event2, n2, studlab = studlab,
                method = method, sm = x$sm,
                incr = incr, method.incr = method.incr,
                allstudies = allstudies,
                subset = trt1 == comps$trt1[i] & trt2 == comps$trt2[i],
                complab = comp.i,
                common = common,
                random = random,
                level = level,
                level.ma = level.ma,
                prediction = prediction,
                level.predict = level.predict,
                method.tau = method.tau,
                label.e = comps$trt1[i], label.c = comps$trt2[i],
                warn = FALSE,
                warn.deprecated = FALSE,
                ...)
    }
    #
    attr(res, "order") <- order
    attr(res, "version") <- packageDescription("netmeta")$Version
    #
    class(res) <- "netpairwise"
  }
  
  res
}


#' @rdname netpairwise
#' @method print netpairwise
#' @export

print.netpairwise <- function(x, ...) {
  
  chkclass(x, "netpairwise")
  
  if (inherits(x, "metagen")) {
    print(x, ...)
  }
  else {
    n <- 1
    for (i in 1:length(x)) {
      if (n > 1)
        cat("\n*****\n\n")
      print(x[[i]], ...)
      n <- n + 1
    }
  }
  
  invisible(NULL)
}


#' @rdname netpairwise
#' @method summary netpairwise
#' @export

summary.netpairwise <- function(object, ...) {
  
  chkclass(object, "netpairwise")
  
  if (inherits(object, "metagen")) {
    res <- summary(object)
  }
  else {
    res <- object
    for (i in seq_len(length(object)))
      res[[i]] <- summary(object[[i]])
    ##
    class(res) <- "summary.netpairwise"
  }
  
  res
}


#' @rdname netpairwise
#' @method print summary.netpairwise
#' @export

print.summary.netpairwise <- function(x, ...) {
  
  chkclass(x, "summary.netpairwise")
  
  if (inherits(x, "metagen")) {
    print(x, ...)
  }
  else {
    n <- 1
    for (i in seq_len(length(x))) {
      if (n > 1)
        cat("\n*****\n\n")
      print(x[[i]], ...)
      n <- n + 1
    }
  }
  
  invisible(NULL)
}


#' @rdname netpairwise
#' @method forest netpairwise
#' @export

forest.netpairwise <- function(x, ...) {
  
  chkclass(x, "netpairwise")
  
  if (inherits(x, "metagen")) {
    res <- forest(x, ...)
    return(invisible(res))
  }
  else {
    for (i in seq_len(length(x))) {
      m.i <- x[[i]]
      forest(m.i, smlab = paste0("Comparison:\n", m.i$complab), ...)
    }
  }
  
  invisible(NULL)
}


#' @rdname netpairwise
#' @method plot netpairwise
#' @export

plot.netpairwise <- function(x, ...)
  forest(x, ...)


#' @rdname netpairwise
#' @method funnel netpairwise
#' @export

funnel.netpairwise <- function(x, k.min = 3, ...) {
  
  chkclass(x, "netpairwise")

  if (inherits(x, "metagen")) {
     stop("Funnel plot not suitable for object of class \"netpairwise\" ",
          "without argument 'separate = TRUE'")
  }
  else {
    for (i in seq_len(length(x))) {
      m.i <- x[[i]]
      if (m.i$k >= k.min) {
        funnel(m.i, ...)
        title(m.i$complab)
      }
    }
  }
  
  invisible(NULL)
}


#' @rdname netpairwise
#' @method radial netpairwise
#' @export

radial.netpairwise <- function(x, k.min = 3, ...) {
  
  chkclass(x, "netpairwise")
  
  if (inherits(x, "metagen")) {
     stop("Radial plot not suitable for object of class \"netpairwise\" ",
          "without argument 'separate = TRUE'")
  }
  else {
    for (i in seq_len(length(x))) {
      m.i <- x[[i]]
      if (m.i$k >= k.min) {
        radial(m.i, ...)
        title(m.i$complab)
      }
    }
  }
  
  invisible(NULL)
}


#' @rdname netpairwise
#' @method baujat netpairwise
#' @export

baujat.netpairwise <- function(x, k.min = 3, ...) {
  
  chkclass(x, "netpairwise")

  if (inherits(x, "metagen")) {
     stop("Baujat plot not suitable for object of class \"netpairwise\" ",
          "without argument 'separate = TRUE'")
  }
  else {
    for (i in seq_len(length(x))) {
      m.i <- x[[i]]
      if (m.i$k >= k.min) {
        baujat(m.i, ...)
        title(m.i$complab)
      }
    }
  }
  
  invisible(NULL)
}


#' @rdname netpairwise
#' @method metabias netpairwise
#' @export

metabias.netpairwise <- function(x, k.min = 10, ...) {
  
  chkclass(x, "netpairwise")
  res <- list()

  if (inherits(x, "metagen")) {
     stop("Test for funnel plot asymmetry not suitable for object of ",
          "class \"netpairwise\" without argument 'separate = TRUE'.")
  }
  else {
    n <- 0
    for (i in seq_len(length(x))) {
      m.i <- x[[i]]
      if (m.i$k >= k.min) {
        n <- n + 1
        res[[n]] <- metabias(m.i, k.min = k.min, ...)
      }
    }
    ##
    if (n > 0)
      class(res) <- "metabias.netpairwise"
    else {
      warning("No pairwise comparison with at least ", k.min, " studies.",
              call. = FALSE)
      return(invisible(NULL))
    }
  }
  
  res
}


#' @rdname netpairwise
#' @method print metabias.netpairwise
#' @export

print.metabias.netpairwise <- function(x, ...) {
  
  chkclass(x, "metabias.netpairwise")
  
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
#' @method trimfill netpairwise
#' @export

trimfill.netpairwise <- function(x, k.min = 3, ...) {
  
  chkclass(x, "netpairwise")
  res <- list()
  
  if (inherits(x, "metagen")) {
     stop("Trim-and-fill method not suitable for object of ",
          "class \"netpairwise\" without argument 'separate = TRUE'.")
  }
  else {
    n <- 0
    for (i in seq_len(length(x))) {
      m.i <- x[[i]]
      if (m.i$k >= k.min) {
        n <- n + 1
        res[[n]] <- trimfill(m.i, ...)
      }
    }
    ##
    if (n > 0)
      class(res) <- c("trimfill.netpairwise", "netpairwise")
    else {
      warning("No pairwise comparison with least ", k.min, " studies.",
              call. = FALSE)
      res <- NULL
    }
  }
  
  res
}


#' @rdname netpairwise
#' @method print trimfill.netpairwise
#' @export

print.trimfill.netpairwise <- function(x, ...) {
  
  chkclass(x, "trimfill.netpairwise")
  
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
#' @method metainf netpairwise
#' @export

metainf.netpairwise <- function(x, k.min = 2, ...) {
  
  chkclass(x, "netpairwise")
  res <- list()
  
  if (inherits(x, "metagen")) {
     stop("Trim-and-fill method not suitable for object of ",
          "class \"netpairwise\" without argument 'separate = TRUE'.")
  }
  else {
    n <- 0
    for (i in seq_len(length(x))) {
      m.i <- x[[i]]
      if (m.i$k >= k.min) {
        n <- n + 1
        res[[n]] <- metainf(m.i, ...)
      }
    }
    ##
    if (n > 0)
      class(res) <- c("metainf.netpairwise", "netpairwise")
    else {
      warning("No pairwise comparison with least ", k.min, " studies.",
              call. = FALSE)
      res <- NULL
    }
  }
  
  ##
  res
}


#' @rdname netpairwise
#' @method print metainf.netpairwise
#' @export

print.metainf.netpairwise <- function(x, ...) {
  
  chkclass(x, "metainf.netpairwise")
  
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
#' @method metacum netpairwise
#' @export

metacum.netpairwise <- function(x, k.min = 2, ...) {
  
  chkclass(x, "netpairwise")
  res <- list()
  
  if (inherits(x, "metagen")) {
     stop("Trim-and-fill method not suitable for object of ",
          "class \"netpairwise\" without argument 'separate = TRUE'.")
  }
  else {
    n <- 0
    for (i in seq_len(length(x))) {
      m.i <- x[[i]]
      if (m.i$k >= k.min) {
        n <- n + 1
        res[[n]] <- metacum(m.i, ...)
      }
    }
    ##
    if (n > 0)
      class(res) <- c("metacum.netpairwise", "netpairwise")
    else {
      warning("No pairwise comparison with least ", k.min, " studies.",
              call. = FALSE)
      res <- NULL
    }
  }
  
  res
}


#' @rdname netpairwise
#' @method print metacum.netpairwise
#' @export

print.metacum.netpairwise <- function(x, ...) {
  
  chkclass(x, "metacum.netpairwise")
  
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
#' @method metareg netpairwise
#' @export

metareg.netpairwise <- function(x, ..., k.min = 2) {
  
  chkclass(x, "netpairwise")
  res <- list()
  
  if (inherits(x, "metagen")) {
     stop("Trim-and-fill method not suitable for object of ",
          "class \"netpairwise\" without argument 'separate = TRUE'.")
  }
  else {
    n <- 0
    for (i in seq_len(length(x))) {
      m.i <- x[[i]]
      if (m.i$k >= k.min) {
        n <- n + 1
        res[[n]] <- metareg(m.i, ...)
      }
    }
    ##
    if (n > 0)
      class(res) <- c("metareg.netpairwise", "netpairwise")
    else {
      warning("No pairwise comparison with least ", k.min, " studies.",
              call. = FALSE)
      res <- NULL
    }
  }
  
  res
}


#' @rdname netpairwise
#' @method print metareg.netpairwise
#' @export

print.metareg.netpairwise <- function(x, ...) {
  
  chkclass(x, "metareg.netpairwise")
  
  n <- 1
  for (i in seq_len(length(x))) {
    if (n > 1)
      cat("\n*****\n\n")
    print(x[[i]], ...)
    n <- n + 1
  }
  
  invisible(NULL)
}
