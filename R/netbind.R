#' Combine network meta-analysis objects
#' 
#' @description
#' This function can be used to combine network meta-analysis objects
#' which is especially useful to generate a forest plot with results
#' of several network meta-analyses.
#' 
#' @param ... Any number of meta-analysis objects or a single list
#'   with network meta-analyses.
#' @param name An optional character vector providing descriptive
#'   names for the network meta-analysis objects.
#' @param comb.fixed A logical indicating whether results for the
#'   fixed effects (common effects) model should be reported.
#' @param comb.random A logical indicating whether results for the
#'   random effects model should be reported.
#' @param col.study The colour for network estimates and confidence
#'   limits.
#' @param col.inside The colour for network estimates and confidence
#'   limits if confidence limits are completely within squares.
#' @param col.square The colour for squares.
#' @param col.square.lines The colour for the outer lines of squares.
#' @param backtransf A logical indicating whether results should be
#'   back transformed. If \code{backtransf = TRUE} (default), results
#'   for \code{sm = "OR"} are printed as odds ratios rather than log
#'   odds ratios, for example.
#' @param reference.group Reference treatment.
#' @param baseline.reference A logical indicating whether results
#'   should be expressed as comparisons of other treatments versus the
#'   reference treatment (default) or vice versa. This argument is
#'   only considered if \code{reference.group} has been specified.
#' 
#' @return
#' An object of class \code{"netbind"} with corresponding
#' \code{forest} function. The object is a list containing the
#' following components:
#' \item{fixed}{A data frame with results for the fixed effects
#'   model.}
#' \item{random}{A data frame with results for the random effects
#'   model.}
#' \item{sm}{Summary measure used in network meta-analyses.}
#' \item{level.comb}{Level for confidence intervals.}
#' \item{comb.fixed, comb.random, backtransf}{As defined above.}
#' \item{reference.group, baseline.reference}{As defined above.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netcomb}},
#'   \code{\link{discomb}}, \code{\link{forest.netbind}}
#' 
#' @examples
#' data(Linde2016)
#' 
#' # Only consider studies including Face-to-face PST (to reduce
#' # runtime of example)
#' #
#' face <- subset(Linde2016, id %in% c(16, 24, 49, 118))
#' 
#' # Standard random effects NMA model (with placebo as reference
#' # treatment)
#' #
#' net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'                 data = face, reference.group = "placebo",
#'                 sm = "OR", comb.fixed = FALSE)
#' 
#' # Additive CNMA model with placebo as inactive component and
#' # reference
#' #
#' nc1 <- netcomb(net1, inactive = "placebo")
#' 
#' # Combine results of standard NMA and CNMA
#' #
#' nb1 <- netbind(nc1, net1,
#'                name = c("Additive CNMA", "Standard NMA"),
#'                col.study = c("red", "black"),
#'                col.square = c("red", "black"))
#' forest(nb1,
#'        col.by = "black", addrow.subgroups = FALSE,
#'        fontsize = 10, spacing = 0.7, squaresize = 0.9,
#'        label.left = "Favours Placebo",
#'        label.right = "Favours other")
#' 
#' @export netbind


netbind <- function(..., name,
                    comb.fixed, comb.random,
                    ##
                    col.study = "black",
                    col.inside = "white",
                    col.square = "gray",
                    col.square.lines = col.square,
                    ##
                    backtransf,
                    reference.group, baseline.reference
                    ) {
  
  
  ##
  ##
  ## (1) Extract list elements and basic checks
  ##
  ##
  is.nma <- function(x)
    inherits(x, "netmeta") |
      inherits(x, "netcomb") |
      inherits(x, "discomb")
  ##
  chklogical <- meta:::chklogical
  ##
  args <- list(...)
  ##
  n.netmeta <- length(args)
  n.i <- seq_len(n.netmeta)
  ##
  if (length(args) == 1) {
    if (!is.list(args[[1]]))
      stop("All elements of argument '...' must be of classes ",
           "'netmeta', 'netcomb', or 'discomb'.",
           call. = FALSE)
    ##
    if (!is.nma(args[[1]])) {
      n.netmeta <- length(args[[1]])
      n.i <- seq_len(n.netmeta)
      ##
      args2 <- list()
      for (i in n.i)
        args2[[i]] <- args[[1]][[i]]
    }
    args <- args2
  }
  ##  
  for (i in n.i) {
    if (!is.nma(args[[i]]))
      stop("All elements of argument '...' must be of classes ",
           "'netmeta', 'netcomb', or 'discomb'.",
           call. = FALSE)
  }
  ##
  levs <- numeric(0)
  for (i in n.i)
    levs[i] <- args[[i]]$level.comb
  ##
  if (length(unique(levs)) != 1)
    stop("Different confidence levels used in network meta-analyses ",
         "(see list element 'level.comb').",
         call. = FALSE)
  ##
  sms <- character(0)
  for (i in n.i)
    sms[i] <- args[[i]]$sm
  ##
  if (length(unique(sms)) != 1)
    stop("Different summary measure used in network meta-analyses ",
         "(see list element 'sm').",
         call. = FALSE)
  ##
  if (missing(col.study))
    col.study <- rep(col.study, n.netmeta)
  else
    if (length(col.study) != n.netmeta) 
      stop("Length of argument 'col.study' must be the same as the ",
           "number of network meta-analyses.",
           call. = FALSE)
  ##
  if (missing(col.square))
    col.square <- rep(col.square, n.netmeta)
  else
    if (length(col.square) != n.netmeta) 
      stop("Length of argument 'col.square' must be the same as the ",
           "number of network meta-analyses.",
           call. = FALSE)
  ##
  if (missing(col.square.lines))
    col.square.lines <- rep(col.square.lines, n.netmeta)
  else
    if (length(col.square.lines) != n.netmeta) 
      stop("Length of argument 'col.square.lines' must be the same as the ",
           "number of network meta-analyses.",
           call. = FALSE)
  ##
  if (missing(col.inside))
    col.inside <- rep(col.inside, n.netmeta)
  else
    if (length(col.inside) != n.netmeta) 
      stop("Length of argument 'col.inside' must be the same as the ",
           "number of network meta-analyses.",
           call. = FALSE)
  ##
  print.warning1 <- FALSE
  print.warning2 <- FALSE
  print.warning3 <- FALSE
  print.warning4 <- FALSE
  print.warning5 <- FALSE
  print.warning6 <- FALSE


  ##
  ##
  ## (2) Set names of network meta-analysis objects
  ##
  ##
  if (missing(name))
    name <- paste0("netmeta", n.i)
  else {
    if (length(name) != n.netmeta) 
      stop("Number of network meta-analyses and ",
           "names provided in argument 'name' differ.",
           call. = FALSE)
    ##
    if (length(unique(name)) != length(name)) {
      warning1 <- paste0("Network meta-analyses are labelled 'netmeta1' to 'netmeta",
                         n.netmeta,
                         "' as values of argument 'name' are not all disparate.")
      print.warning1 <- TRUE
      ##
      name <- paste0("netmeta", n.i)
    }
  }
  
  
  ##
  ## (3) Determine comb.fixed
  ##
  if (missing(comb.fixed)) {
    cfs <- logical(0)
    ##
    for (i in n.i)
      cfs[i] <- args[[i]]$comb.fixed
    ##
    cfs <- unique(cfs)
    ##
    if (length(cfs) != 1) {
      comb.fixed <- TRUE
      warning2 <- paste0("Argument 'comb.fixed' set to TRUE ",
                         "(as it is not unique in network meta-analyses).")
      print.warning2 <- TRUE
    }
    else
      comb.fixed <- cfs
  }
  else
    chklogical(comb.fixed)
  
  
  ##
  ## (4) Determine comb.random
  ##
  if (missing(comb.random)) {
    crs <- logical(0)
    ##
    for (i in n.i)
      crs[i] <- args[[i]]$comb.random
    ##
    crs <- unique(crs)
    ##
    if (length(crs) != 1) {
      comb.random <- TRUE
      warning3 <- paste0("Argument 'comb.random' set to TRUE ",
                         "(as it is not unique in network meta-analyses).")
      print.warning3 <- TRUE
    }
    else
      comb.random <- crs
  }
  else
    chklogical(comb.random)
  
  
  ##
  ## (5) Determine backtransf
  ##
  if (missing(backtransf)) {
    backt <- logical(0)
    ##
    for (i in n.i)
      backt[i] <- args[[i]]$backtransf
    ##
    backt <- unique(backt)
    ##
    if (length(backt) != 1) {
      backtransf <- TRUE
      warning4 <- paste0("Argument 'backtransf' set to TRUE ",
                         "(as it is not unique in network meta-analyses).")
      print.warning4 <- TRUE
    }
    else
      backtransf <- backt
  }
  else
    chklogical(backtransf)
  
  
  ##
  ##
  ## (4) Determine reference group
  ##
  ##
  if (missing(reference.group)) {
    refs <- character(0)
    ##
    for (i in n.i)
      refs[i] <- args[[i]]$reference.group
    ##
    refs <- unique(refs)
    refs <- refs[refs != ""]
    ##
    if (length(refs) == 0) {
      reference.group <- args[[i]]$trts[1]
      warning5 <- paste0("Unspecified argument 'reference.group' is set to '",
                         reference.group,
                         "'.")
      print.warning5 <- TRUE
    }
    else if (length(refs) == 1)
      reference.group <- refs
    else {
      reference.group <- refs[1]
      warning5 <- paste0("Argument 'reference.group' set to '",
                         reference.group,
                         "' (as it is not unique in network meta-analyses).")
      print.warning5 <- TRUE
    }
  }
  ##
  trts.all <- character(0)
  ##
  for (i in n.i)
    trts.all <- c(trts.all, args[[i]]$trts)
  ##
  trts.all <- unique(trts.all)
  ##
  reference.group <- setref(reference.group, trts.all)
  ##
  for (i in n.i)
    if (!(reference.group %in% args[[i]]$trts))
      stop("Reference treatment '", reference.group,
           "' not included in all network meta-analyses.")
  
  
  ##
  ## (6) Determine baseline reference
  ##
  if (missing(baseline.reference)) {
    bref <- logical(0)
    ##
    for (i in n.i)
      bref[i] <- args[[i]]$baseline.reference
    ##
    bref <- unique(bref)
    ##
    if (length(bref) != 1) {
      bref <- TRUE
      warning6 <- paste0("Argument 'baseline.reference' set to TRUE ",
                         "(as it is not unique in network meta-analyses).")
      print.warning6 <- TRUE
    }
    ##
    baseline.reference <- bref
  }
  ##
  chklogical(baseline.reference)
  
  
  fixed <- data.frame(name = character(0),
                      treat = character(0),
                      TE = numeric(0), seTE = numeric(0),
                      lower = numeric(0), upper = numeric(0),
                      zval = numeric(0), pval = numeric(0),
                      ##
                      col.study = character(0),
                      col.square = character(0),
                      col.square.lines = character(0),
                      col.inside = character(0),
                      ##
                      stringsAsFactors = FALSE)
  ##
  for (i in n.i) {
    ##
    rn <- rownames(args[[i]]$TE.fixed)
    seq1 <- charmatch(setseq(args[[i]]$seq, rn), rn)
    ##
    TE.i    <- args[[i]]$TE.fixed[seq1, seq1]
    seTE.i  <- args[[i]]$seTE.fixed[seq1, seq1]
    lower.i <- args[[i]]$lower.fixed[seq1, seq1]
    upper.i <- args[[i]]$upper.fixed[seq1, seq1]
    zval.i  <- args[[i]]$zval.fixed[seq1, seq1]
    pval.i  <- args[[i]]$pval.fixed[seq1, seq1]
    ##
    cnam <- colnames(TE.i)
    rnam <- rownames(TE.i)
    selc <- cnam == reference.group
    selr <- rnam == reference.group
    ##
    if (baseline.reference) {
      fixed <- rbind(fixed,
                     data.frame(name = name[i],
                                treat = cnam,
                                ##
                                TE = TE.i[, selc],
                                seTE = seTE.i[, selc],
                                lower = lower.i[, selc],
                                upper = upper.i[, selc],
                                zval = zval.i[, selc],
                                pval = pval.i[, selc],
                                ##
                                col.study = col.study[i],
                                col.square = col.square[i],
                                col.square.lines = col.square.lines[i],
                                col.inside = col.inside[i],
                                ##
                                stringsAsFactors = FALSE)
                     )
    }
    else {
      fixed <- rbind(fixed,
                     data.frame(name = name[i],
                                treat = cnam,
                                ##
                                TE = TE.i[selr, ],
                                seTE = seTE.i[selr, ],
                                lower = lower.i[selr, ],
                                upper = upper.i[selr, ],
                                zval = zval.i[selr, ],
                                pval = pval.i[selr, ],
                                ##
                                col.study = col.study[i],
                                col.square = col.square[i],
                                col.square.lines = col.square.lines[i],
                                col.inside = col.inside[i],
                                ##
                                stringsAsFactors = FALSE)
                     )
    }
  }
  ##
  rownames(fixed) <- seq_len(nrow(fixed))
  
  
  random <- data.frame(name = character(0),
                       treat = character(0),
                       ##
                       TE = numeric(0), seTE = numeric(0),
                       lower = numeric(0), upper = numeric(0),
                       zval = numeric(0), pval = numeric(0),
                       ##
                       col.study = character(0),
                       col.square = character(0),
                       col.square.lines = character(0),
                       col.inside = character(0),
                       ##
                       stringsAsFactors = FALSE)
  ##
  for (i in n.i) {
    ##
    rn <- rownames(args[[i]]$TE.random)
    seq1 <- charmatch(setseq(args[[i]]$seq, rn), rn)
    ##
    TE.i    <- args[[i]]$TE.random[seq1, seq1]
    seTE.i  <- args[[i]]$seTE.random[seq1, seq1]
    lower.i <- args[[i]]$lower.random[seq1, seq1]
    upper.i <- args[[i]]$upper.random[seq1, seq1]
    zval.i  <- args[[i]]$zval.random[seq1, seq1]
    pval.i  <- args[[i]]$pval.random[seq1, seq1]
    ##
    cnam <- colnames(TE.i)
    rnam <- rownames(TE.i)
    selc <- cnam == reference.group
    selr <- rnam == reference.group
    ##
    if (baseline.reference) {
      random <- rbind(random,
                      data.frame(name = name[i],
                                 treat = cnam,
                                 ##
                                 TE = TE.i[, selc],
                                 seTE = seTE.i[, selc],
                                 lower = lower.i[, selc],
                                 upper = upper.i[, selc],
                                 zval = zval.i[, selc],
                                 pval = pval.i[, selc],
                                 ##
                                 col.study = col.study[i],
                                 col.square = col.square[i],
                                 col.square.lines = col.square.lines[i],
                                 col.inside = col.inside[i],
                                 ##
                                 stringsAsFactors = FALSE)
                      )
    }
    else {
      random <- rbind(random,
                      data.frame(name = name[i],
                                 treat = cnam,
                                 ##
                                 TE = TE.i[selr, ],
                                 seTE = seTE.i[selr, ],
                                 lower = lower.i[selr, ],
                                 upper = upper.i[selr, ],
                                 zval = zval.i[selr, ],
                                 pval = pval.i[selr, ],
                                 ##
                                 col.study = col.study[i],
                                 col.square = col.square[i],
                                 col.square.lines = col.square.lines[i],
                                 col.inside = col.inside[i],
                                 ##
                                 stringsAsFactors = FALSE)
                      )
    }
  }
  ##
  rownames(random) <- seq_len(nrow(random))
  
  
  res <- list(fixed = fixed,
              random = random,
              sm = sms[1],
              level.comb = levs[1],
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              backtransf = backtransf,
              reference.group = reference.group,
              baseline.reference = baseline.reference)
  ##
  class(res) <- "netbind"


  ##
  ## Print warnings
  ##
  if (print.warning1)
    warning(warning1, call. = FALSE)
  if (print.warning2)
    warning(warning2, call. = FALSE)
  if (print.warning3)
    warning(warning3, call. = FALSE)
  if (print.warning4)
    warning(warning4, call. = FALSE)
  if (print.warning5)
    warning(warning5, call. = FALSE)
  if (print.warning6)
    warning(warning6, call. = FALSE)


  res
}
