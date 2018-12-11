netmetabin <- function(event1, n1, event2, n2,
                       treat1, treat2, studlab,
                       data = NULL, subset = NULL,
                       sm,
                       method = "MH",
                       cc.pooled = FALSE,
                       incr = gs("incr"),
                       allincr = gs("allincr"), addincr = gs("addincr"),
                       allstudies = gs("allstudies"),
                       level = gs("level"),
                       level.comb = gs("level.comb"),
                       comb.fixed = gs("comb.fixed"),
                       comb.random = method == "Inverse" &
                         (gs("comb.random") | !is.null(tau.preset)),
                       ##
                       prediction = FALSE,
                       level.predict = gs("level.predict"),
                       ##
                       reference.group = "",
                       baseline.reference = TRUE,
                       all.treatments = NULL,
                       seq = NULL,
                       ##
                       tau.preset = NULL,
                       ##
                       tol.multiarm = 0.0005,
                       details.chkmultiarm = FALSE,
                       ##
                       sep.trts = ":",
                       nchar.trts = 666,
                       ##
                       backtransf = gs("backtransf"),
                       ##
                       title = "",
                       keepdata = gs("keepdata"),
                       warn = TRUE) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkchar <- meta:::chkchar
  chklevel <- meta:::chklevel
  chklogical <- meta:::chklogical
  chknull <- meta:::chknull
  chknumeric <- meta:::chknumeric
  setchar <- meta:::setchar
  pvalQ <- meta:::pvalQ
  ##
  modtext <- paste0("must be equal to 'Inverse' (classic network meta-analysis), ",
                    "'MH' (Mantel-Haenszel, the default) or ",
                    "'NCH' (common-effects non-central hypergeometric).")
  method <- meta:::setchar(method, c("Inverse", "MH", "NCH"), modtext)
  ##
  chklogical(allincr)
  chklogical(addincr)
  chklogical(allstudies)
  chklogical(cc.pooled)
  ##
  chklevel(level)
  chklevel(level.comb)
  chklevel(level.predict)
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
  ##
  if (method != "Inverse" & !comb.fixed) {
    warning("Argument 'comb.fixed' set to TRUE for Mantel-Haenszel ",
            "method and non-central hypergeometric distribution.",
            call. = FALSE)
    comb.fixed <- TRUE
  }
  ##
  if (method != "Inverse" & comb.random) {
    warning("Argument 'comb.random' set to FALSE for Mantel-Haenszel ",
            "method and non-central hypergeometric distribution.",
            call. = FALSE)
    comb.random <- FALSE
  }
  ##
  if (method != "Inverse" & prediction) {
    warning("Argument 'prediction' set to FALSE for Mantel-Haenszel ",
            "method and non-central hypergeometric distribution.",
            call. = FALSE)
    prediction <- FALSE
  }
  ##
  chklogical(baseline.reference)
  ##
  if (!is.null(all.treatments))
    chklogical(all.treatments)
  ##
  if (!is.null(tau.preset))
    chknumeric(tau.preset, min = 0, single = TRUE)
  ##
  chknumeric(tol.multiarm, min = 0, single = TRUE)
  chklogical(details.chkmultiarm)
  ##
  missing.sep.trts <- missing(sep.trts)
  chkchar(sep.trts)
  chknumeric(nchar.trts, min = 1, single = TRUE)
  ##
  chklogical(backtransf)
  ##
  chkchar(title)
  chklogical(keepdata)
  chklogical(warn)
  ##
  ## Check value for reference group
  ##
  if (is.null(all.treatments))
    if (reference.group == "")
      all.treatments <- TRUE
    else
      all.treatments <- FALSE
  ##
  chklogical(baseline.reference)


  ##
  ##
  ## (2) Read data
  ##
  ##
  nulldata <- is.null(data)
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch 'event1', 'event2', 'n1', 'n2', 'treat1', 'treat2', and
  ## 'studlab' from data:
  ##
  event1 <- eval(mf[[match("event1", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  if (inherits(event1, "pairwise")) {
    is.pairwise <- TRUE
    ##
    if (missing(sm) & method == "Inverse")
      sm <- attr(event1, "sm")
    else if (method != "Inverse") {
      if (!missing(sm) && tolower(sm) != "or")
        warning("Argument 'sm' set to 'OR'.", call. = FALSE)
      sm <- "OR"
    }
    ##
    n1 <- event1$n1
    event2 <- event1$event2
    n2 <- event1$n2
    treat1 <- event1$treat1
    treat2 <- event1$treat2
    studlab <- event1$studlab
    ##
    pairdata <- event1
    data <- event1
    ##
    event1 <- event1$event1
  }
  else {
    is.pairwise <- FALSE
    ##
    if (missing(sm) & method == "Inverse") {
      if (!nulldata && !is.null(attr(data, "sm")))
        sm <- attr(data, "sm")
      else
        sm <- ""
    }
    else if (method != "Inverse") {
      if (!missing(sm) && tolower(sm) != "or")
        warning("Argument 'sm' set to 'OR'.", call. = FALSE)
      sm <- "OR"
    }
    ##
    n1 <- eval(mf[[match("n1", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
    ##
    event2 <- eval(mf[[match("event2", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
    ##
    n2 <- eval(mf[[match("n2", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
    ##
    treat1 <- eval(mf[[match("treat1", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
    ##
    treat2 <- eval(mf[[match("treat2", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
    ##
    studlab <- eval(mf[[match("studlab", names(mf))]],
                    data, enclos = sys.frame(sys.parent()))
  }
  ##
  k.Comp <- length(event1)
  ##
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  if (is.factor(studlab))
    studlab <- as.character(studlab)
  ##
  ## Keep original order of studies
  ##
  .order <- seq_along(studlab)
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  missing.subset <- is.null(subset)
  
  
  ##
  ##
  ## (2b) Store complete dataset in list object data
  ##      (if argument keepdata is TRUE)
  ##
  ##
  if (keepdata) {
    if (nulldata & !is.pairwise)
      data <- data.frame(.event1 = event1)
    else if (nulldata & is.pairwise) {
      data <- pairdata
      data$.order <- .order
      data$.event1 <- event1
    }
    else
      data$.event1 <- event1
    ##
    data$.n1 <- n1
    data$.event2 <- event2
    data$.n2 <- n2
    data$.treat1 <- treat1
    data$.treat2 <- treat2
    data$.studlab <- studlab
    data$.order <- .order
    ##
    if (!missing.subset) {
      if (length(subset) == dim(data)[1])
        data$.subset <- subset
      else {
        data$.subset <- FALSE
        data$.subset[subset] <- TRUE
      }
    }
  }
  
  
  ##
  ##
  ## (3) Use subset for analysis
  ##
  ##
  if (!missing.subset) {
    if ((is.logical(subset) & (sum(subset) > k.Comp)) ||
        (length(subset) > k.Comp))
      stop("Length of subset is larger than number of studies.",
           call. = FALSE)
    ##
    studlab <- studlab[subset]
    treat1 <- treat1[subset]
    treat2 <- treat2[subset]
    ##
    event1 <- event1[subset]
    event2 <- event2[subset]
    n1 <- n1[subset]
    n2 <- n2[subset]
    ##
    .order <- .order[subset]
  }
  ##
  ## Check for correct treatment order within comparison
  ##
  wo <- treat1 > treat2
  ##
  if (any(wo)) {
    if (warn)
      warning("Note, treatments within a comparison have been ",
              "re-sorted in increasing order.",
              call. = FALSE)
    ##
    tevent1 <- event1
    event1[wo] <- event2[wo]
    event2[wo] <- tevent1[wo]
    ##
    tn1 <- n1
    n1[wo] <- n2[wo]
    n2[wo] <- tn1[wo]
    ##
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
  }
  ##
  labels <- sort(unique(c(treat1, treat2)))
  ##
  if (compmatch(labels, sep.trts)) {
    if (!missing.sep.trts)
      warning("Separator '", sep.trts, "' used in at least ",
              "one treatment label. Try to use predefined separators: ",
              "':', '-', '_', '/', '+', '.', '|', '*'.",
              call. = FALSE)
    ##
    if (!compmatch(labels, ":"))
      sep.trts <- ":"
    else if (!compmatch(labels, "-"))
      sep.trts <- "-"
    else if (!compmatch(labels, "_"))
      sep.trts <- "_"
    else if (!compmatch(labels, "/"))
      sep.trts <- "/"
    else if (!compmatch(labels, "+"))
      sep.trts <- "+"
    else if (!compmatch(labels, "."))
      sep.trts <- "-"
    else if (!compmatch(labels, "|"))
      sep.trts <- "|"
    else if (!compmatch(labels, "*"))
      sep.trts <- "*"
    else
      stop("All predefined separators (':', '-', '_', '/', '+', '.', '|', '*') ",
           "are used in at least one treatment label.\n   ",
           "Please specify a different character that should be used ",
           " as separator (argument 'sep.trts').",
           call. = FALSE)
  }
  ##
  if (reference.group != "")
    reference.group <- setref(reference.group, labels)
  ##
  rm(labels)


  ##
  ##
  ## (4) Additional checks
  ##
  ##
  if (any(treat1 == treat2))
    stop("Treatments must be different (arguments 'treat1' and 'treat2').",
         call. = FALSE)
  ##
  if (length(studlab) != 0)
    studlab <- as.character(studlab)
  else {
    if (warn)
      warning("No information given for argument 'studlab'. ",
              "Assuming that comparisons are from independent studies.",
              call. = FALSE)
    studlab <- seq(along = event1)
  }
  ##
  ## Check for correct number of comparisons
  ##
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)
      abs(x - round(x)) < tol
  ##
  tabnarms <- table(studlab)
  sel.narms <- !is.wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  ##
  if (sum(sel.narms) == 1)
    stop(paste0("Study '", names(tabnarms)[sel.narms],
                "' has a wrong number of comparisons.",
                "\n  Please provide data for all treatment comparisons",
                " (two-arm: 1; three-arm: 3; four-arm: 6, ...)."),
         call. = FALSE)
  if (sum(sel.narms) > 1)
    stop(paste0("The following studies have a wrong number of comparisons: ",
                paste(paste0("'", names(tabnarms)[sel.narms], "'"),
                      collapse = " - "),
                "\n  Please provide data for all treatment comparisons",
                " (two-arm: 1; three-arm: 3; four-arm: 6, ...)."),
         call. = FALSE)
  ##
  ## Check number of subgraphs
  ##
  n.subnets <- netconnection(treat1, treat2, studlab)$n.subnets
  ##
  if (n.subnets > 1)
    stop(paste0("Network consists of ", n.subnets, " separate sub-networks.\n",
                "  Use R function 'netconnection' to identify sub-networks."),
         call. = FALSE)
  ##
  ## Catch 'incr' from data:
  ##
  if (!missing(incr))
    incr <- eval(mf[[match("incr", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknumeric(incr, min = 0)
  ##
  if (method != "Inverse" & length(incr) > 1) {
    warning("Argument 'incr' must be a single value for ",
            "Mantel-Haenszel and common-effects non-central ",
            "hypergeometric method. Set to zero (default).",
            call. = FALSE)
    incr <- 0
  }
  ##
  if (all(incr == 0) & method != "Inverse" & cc.pooled == TRUE)
    cc.pooled <- FALSE
  ##
  lengthunique <- function(x) length(unique(x))


  ##
  ##
  ## (5) Create analysis dataset
  ##
  ##
  dat.long <- rbind(data.frame(studlab, treat = treat1, event = event1, n = n1,
                               .order,
                               stringsAsFactors = FALSE),
                    data.frame(studlab, treat = treat2, event = event2, n = n2,
                               .order,
                               stringsAsFactors = FALSE))
  ##
  dat.wide <- data.frame(studlab = studlab, treat1 = treat1, treat2 = treat2,
                         event1 = event1, n1 = n1, event2 = event2, n2 = n2,
                         .order,
                         stringsAsFactors = FALSE)
  dat.wide <- dat.wide[order(dat.wide$studlab,
                             dat.wide$treat1, dat.wide$treat2), ]
  ##
  get.designs <- function(x) {
    ##
    res <- nma.krahn(netmeta(pairwise(studlab = x$studlab,
                                      treat = x$treat,
                                      event = x$event,
                                      n = x$n,
                                      sm = "RD", warn = FALSE),
                             warn = FALSE)
                     )$studies[, c("studlab", "design")]
    ##
    if (is.null(res)) {
      net1 <- netmeta(pairwise(studlab = x$studlab,
                               treat = x$treat,
                               event = x$event,
                               n = x$n,
                               sm = "RD", warn = FALSE),
                      warn = FALSE)
      ##
      res <- data.frame(studlab = net1$studlab, design = net1$designs)
    }
    ##
    res$design <- as.character(res$design)
    res <- unique(res)
    ##
    res
  }
  ##
  set.designs <- function(x, tdat.design, all.designs) {
    x$design <- NULL
    x <- merge(x, tdat.design, by = "studlab", all = TRUE)
    x$design <- as.character(x$design)
    ##
    if (any(is.na(x$design))) {
      sel.studlab <- unique(x$studlab[is.na(x$design)])
      for (i in sel.studlab)
        x$design[x$studlab == i] <-
          all.designs$design[all.designs$studlab == i]
    }
    x
  }
  ##
  data$.drop <- rep(FALSE, nrow(data))
  ##
  ## Add variable 'non.event'
  ##
  dat.long$non.event <- dat.long$n - dat.long$event
  ##
  dat.wide$non.event1 <- dat.wide$n1 - dat.wide$event1
  dat.wide$non.event2 <- dat.wide$n2 - dat.wide$event2
  ##
  ## Remove duplicated rows from dataset (due to multi-arm studies)
  ##
  dupl <- duplicated(dat.long[, c("studlab", "treat", "event", "n")])
  ##
  if (any(dupl))
    dat.long <- dat.long[!dupl, ]
  ##
  rm(dupl)
  ##
  all.designs <- get.designs(dat.long)
  ##
  ## Add variable 'design' with study design
  ##
  tdat.design <- get.designs(dat.long)
  ##
  dat.long <- merge(dat.long, tdat.design, by = "studlab")
  dat.wide <- merge(dat.wide, tdat.design, by = "studlab")
  ##
  dat.long <- dat.long[order(dat.long$.order), ]
  dat.wide <- dat.wide[order(dat.wide$.order), ]
  ##
  names(tdat.design) <- c("studlab", ".design")
  ##
  data <- merge(data, tdat.design, by = "studlab", all.x = TRUE)
  data <- data[order(data$.order), ]
  ##
  rm(tdat.design)
  
  
  ##
  ##
  ## (6) Stage I: setting up the data (Efthimiou et al., 2018)
  ##
  ##
  ## Step i. Remove all-zero studies (only MH and NCH methods)
  ##
  if (method != "Inverse") {
    n.events <- with(dat.long, tapply(event, studlab, sum))
    ##
    if (any(n.events == 0)) {
      allzero <- n.events == 0
      ##
      if (warn)
        if (sum(allzero) == 1)
          warning("Study '", names(n.events)[allzero],
                  "' without any events excluded from network meta-analysis.",
                  call. = FALSE)
        else
          warning("Studies without any events excluded ",
                  "from network meta-analysis: ",
                  paste(paste0("'", names(n.events)[allzero], "'"),
                        collapse = " - "),
                  call. = FALSE)
      ##
      dat.long <- dat.long[dat.long$studlab %in% names(n.events)[!allzero], ,
                           drop = FALSE]
      dat.wide <- dat.wide[dat.wide$studlab %in% names(n.events)[!allzero], ,
                           drop = FALSE]
      ##
      data$.drop <- data$.drop | data$studlab %in% names(n.events)[allzero]
      ##
      rm(allzero)
    }
    ##
    rm(n.events)
  }
  ##
  ## Add variable 'study' with study numbers
  ##
  dat.long$study <- as.numeric(factor(dat.long$studlab,
                                      levels = unique(dat.long$studlab)))
  dat.wide$study <- as.numeric(factor(dat.wide$studlab,
                                      levels = unique(dat.wide$studlab)))
  ##
  ## Sort by study and treatment
  ##
  dat.long <- dat.long[order(dat.long$studlab, dat.long$treat), ]
  dat.wide <- dat.wide[order(dat.wide$studlab,
                             dat.wide$treat1, dat.wide$treat2), ]
  ##
  dat.long$incr <- 0
  dat.wide$incr <- 0
  data$.incr <- 0
  ##
  ## Step iii. Drop treatment arms without events from individual
  ##           designs (argument 'cc.pooled' is FALSE) or add
  ##           increment if argument 'cc.pooled' is TRUE and argument
  ##           'incr' is larger than zero
  ##
  if (method != "Inverse") {
    ##
    d.events <- with(dat.long, tapply(event, list(design, treat), sum))
    ##
    if (any(d.events == 0, na.rm = TRUE)) {
      zero <- d.events == 0
      zerocells <- as.data.frame(which(zero, arr.ind = TRUE),
                                 stringsAsFactors = FALSE)
      ##
      zerocells$design <- rownames(zero)[zerocells$row]
      zerocells$treat <- colnames(zero)[zerocells$col]
      ##
      zero.long <- rep(0, nrow(dat.long))
      zero.wide <- rep(0, nrow(dat.wide))
      zero.data <- rep(0, nrow(data))
      ##
      for (i in seq_along(zerocells$design)) {
        zero.long <- zero.long + (dat.long$design == zerocells$design[i])
        ##
        zero.wide <- zero.wide + (dat.wide$design == zerocells$design[i] &
                                  (dat.wide$treat1 == zerocells$treat[i] |
                                   dat.wide$treat2 == zerocells$treat[i]))
        ##
        zero.data <- zero.data + (data$.design == zerocells$design[i] &
                                  (data$.treat1 == zerocells$treat[i] |
                                   data$.treat2 == zerocells$treat[i]))
      }
      ##
      zero.long <- zero.long > 0
      zero.wide <- zero.wide > 0
      zero.data <- zero.data > 0
      ##
      if (!cc.pooled) {
        if (warn)
          if (sum(zero, na.rm = TRUE) == 1)
            warning("Treatment arm '", zerocells$treat,
                    "' without events in design '",
                    zerocells$design, "' excluded from network meta-analysis.",
                    call. = FALSE)
          else
            warning("Treatment arms without events in a design excluded ",
                    "from network meta-analysis:\n    ",
                    paste0("'",
                           paste0(paste0(zerocells$treat, " in "),
                                  zerocells$design),
                           "'",
                           collapse = " - "),
                    call. = FALSE)
        ##
        dat.long <- dat.long[!zero.long, , drop = FALSE]
        dat.wide <- dat.wide[!zero.wide, , drop = FALSE]
        data$.drop <- data$.drop | zero.data
        ##
        rm(zero, zerocells, zero.long, zero.wide, zero.data)
      }
      else {
        dat.long$incr[zero.long] <- incr
        ##
        dat.long$event[zero.long] <- dat.long$event[zero.long] + incr
        dat.long$non.event[zero.long] <- dat.long$non.event[zero.long] + incr
        dat.long$n[zero.long] <- dat.long$n[zero.long] + 2 * incr
        ##
        dat.wide$incr[zero.wide] <- incr
        ##
        data$.incr[zero.data] <- incr
        ##
        rm(zero.long)
      }
    }
    ##
    rm(d.events)
    ##
    ## (Re)Add variable 'design' with study design (as treatment arms
    ## may have been dropped)
    ##
    dat.wide$design <- NULL
    data$.design <- NULL
    ##
    tdat.design <- get.designs(dat.long)
    dat.long <- set.designs(dat.long, tdat.design, all.designs)
    all.designs <- get.designs(dat.long)
    ##
    dat.wide <- merge(dat.wide, tdat.design, by = "studlab")
    ##
    dat.long <- dat.long[order(dat.long$.order), ]
    dat.wide <- dat.wide[order(dat.wide$.order), ]
    ##
    names(tdat.design) <- c("studlab", ".design")
    ##
    data <- merge(data, tdat.design, by = "studlab", all.x = TRUE)
    data <- data[order(data$.order), ]
    ##
    rm(tdat.design)
    ##
    dat.long$design <- as.character(dat.long$design)
  }
  ##
  ## Step iv. Remove designs with single treatment arm from dataset
  ##
  if (method != "Inverse") {
    d.single <- with(dat.long, tapply(treat, design, lengthunique))
    ##
    if (any(d.single == 1, na.rm = TRUE)) {
      single <- !is.na(d.single) & d.single == 1
      design.single <- names(d.single)[single]
      ##
      if (warn)
        if (sum(single) == 1)
          warning("Design '", design.single,
                  "' with single treatment arm excluded ",
                  "from network meta-analysis.",
                  call. = FALSE)
        else
          warning("Designs with single treatment arm excluded ",
                  "from network meta-analysis: ",
                  paste(paste0("'", design.single, "'"),
                        collapse = " - "),
                  call. = FALSE)
      ##
      dat.long <- dat.long[!(dat.long$design %in% design.single), , drop = FALSE]
      dat.wide <- dat.wide[!(dat.wide$design %in% design.single), , drop = FALSE]
      ##
      data$.drop <- data$.drop | data$.design %in% design.single
      ##
      rm(single, design.single)
    }
    ##
    rm(d.single)
    ##
    ## (Re)Add variable 'design' with study design (as treatment arms
    ## may have been dropped)
    ##
    dat.wide$design <- NULL
    data$.design <- NULL
    ##
    tdat.design <- get.designs(dat.long)
    dat.long <- set.designs(dat.long, tdat.design, all.designs)
    ##
    dat.wide <- merge(dat.wide, tdat.design, by = "studlab")
    ##
    names(tdat.design) <- c("studlab", ".design")
    ##
    data <- merge(data, tdat.design, by = "studlab", all.x = TRUE)
    data <- data[order(data$.order), ]
    ##
    rm(tdat.design)
    ##
    dat.long <- dat.long[order(dat.long$design, dat.long$studlab,
                               dat.long$treat), ]
    dat.wide <- dat.wide[order(dat.wide$design, dat.wide$studlab,
                               dat.wide$treat1, dat.wide$treat2), ]
  }
  ##
  dat.long <- dat.long[, c("studlab", "treat",
                           "event", "non.event", "n", "incr",
                           "design", "study", ".order")]
  dat.wide <- dat.wide[, c("studlab", "treat1", "treat2",
                           "event1", "event2", "non.event1", "non.event2",
                           "n1", "n2", "incr",
                           "design", "study", ".order")]
  ##
  dat.long$studlab <- as.character(dat.long$studlab)
  dat.long$design <- as.character(dat.long$design)
  dat.long$treat <- as.character(dat.long$treat)
  ##
  dat.wide$studlab <- as.character(dat.wide$studlab)
  dat.wide$design <- as.character(dat.wide$design)
  dat.wide$treat1 <- as.character(dat.wide$treat1)
  dat.wide$treat2 <- as.character(dat.wide$treat2)
  ##
  o <- order(dat.wide$.order)
  studlab <- dat.wide$studlab[o]
  treat1 <- dat.wide$treat1[o]
  treat2 <- dat.wide$treat2[o]
  ##
  event1 <- dat.wide$event1[o]
  event2 <- dat.wide$event2[o]
  n1 <- dat.wide$n1[o]
  n2 <- dat.wide$n2[o]
  ##
  trts <- sort(unique(c(treat1, treat2)))
  n.treat <- length(trts)
  n.studies <- length(unique(studlab))
  m <- length(studlab)
  ##
  treat1.pos <- match(treat1, trts)
  treat2.pos <- match(treat2, trts)
  ##
  ## Calculate adjacency matrix
  ##
  B <- createB(treat1.pos, treat2.pos, ncol = n.treat)
  ##
  M <- t(B) %*% B    # unweighted Laplacian matrix
  D <- diag(diag(M)) # diagonal matrix
  A <- D - M         # adjacency matrix (n x n)
  ##
  rownames(A) <- colnames(A) <- trts
  ##
  rm(B, M, D)
  ##
  ## Empty matrices for results
  ##
  Z <- matrix(NA, nrow = n.treat, ncol = n.treat)
  rownames(Z) <- colnames(Z) <- trts
  ##
  TE.fixed <- seTE.fixed <- Z
  TE.direct.fixed <- seTE.direct.fixed <- Z
  ##
  rm(Z)


  ##
  ##
  ## (7) Conduct classic network meta-analysis using inverse variance
  ##     method
  ##
  ##
  if (method == "Inverse") {
    dat.long <- dat.long[order(dat.long$.order), ]
    ##
    p.iv <- pairwise(studlab = dat.long$studlab,
                     treat = dat.long$treat,
                     event = dat.long$event,
                     n = dat.long$n,
                     data = dat.long,
                     sm = sm,
                     incr = incr,
                     allincr = allincr, addincr = addincr,
                     allstudies = allstudies)
    ##
    net.iv <- netmeta(p.iv,
                      level = level, level.comb = level.comb,
                      comb.fixed = comb.fixed, comb.random = comb.random,
                      prediction = prediction, level.predict = level.predict,
                      reference.group = reference.group,
                      baseline.reference = baseline.reference,
                      all.treatments = all.treatments,
                      seq = seq,
                      tau.preset = tau.preset,
                      tol.multiarm = tol.multiarm,
                      details.chkmultiarm = details.chkmultiarm,
                      sep.trts = sep.trts,
                      nchar.trts = nchar.trts,
                      backtransf = backtransf,
                      title = title,
                      keepdata = keepdata,
                      warn = warn)
    return(net.iv)
  }


  ##
  ##
  ## (8) Stage 2: Direct meta-analyses per design (MH and NCH methods)
  ##
  ##
  n.d <- tapply(dat.long$treat,   dat.long$design, lengthunique)
  k.d <- tapply(dat.long$studlab, dat.long$design, lengthunique)
  ##
  designs <- names(k.d)
  d <- length(designs)
  seq.d <- seq_len(d)
  ##
  dat.design <- vector("list", d)
  names(dat.design) <- designs
  ##
  for (i in seq.d)
    dat.design[[i]] <- dat.long[dat.long$design == designs[i], , drop = FALSE]
  ##
  if (method == "MH") {
    ##
    ## MH method
    ##
    treat.per.design <- vector("list", d)
    studies.in.design <- vector("list", d)
    ##
    c.xy <- vector("list", d)
    d.xy <- vector("list", d)
    C.xy <- vector("list", d)
    L.xy <- vector("list", d)
    ##
    U.xyy <- vector("list", d)
    U.xyz <- vector("list", d)
    ##
    t.pl <- vector("list", d)
    ##
    ps1 <- ps2 <- ps3 <- ps4 <- vector("list", d)
    ##
    L.bar.xy <- vector("list", d)
    ##
    U.plus.xx <- vector("list", d)
    ##
    U.new <- vector("list", d)
    ##
    U.plus.xy <- vector("list", d)
    ##
    Var.Lbar <- vector("list", d)
    ##
    CoVar.Lbar <- vector("list", d)
    ##
    length.y <- sum(n.d - 1)
    ##
    y <- c(rep(0, length.y))
    ##
    V1 <- matrix(rep(0, length.y * length.y), nrow = length.y)
    V2 <- matrix(rep(0, length.y * length.y), nrow = length.y)
    ##
    X <- matrix(0, nrow = length.y, ncol = n.treat - 1)
    ##
    counter1 <- counter2 <- counter3 <- 0
    ##
    list1 <- matrix(0, nrow = length.y, ncol = 2)
    ##
    N.j <- rep(0, d)
    ##
    if (d > 1)
      for (i in 2:d)
        N.j[i] <- N.j[i - 1] + n.d[i - 1] - 1
    ##
    ## Loop for designs
    ##
    for (i in seq.d) {
      treat.per.design[[i]] <- unique(dat.design[[i]]$treat)
      studies.in.design[[i]] <- unique(dat.design[[i]]$study)
      ##
      ## List of studies in each design
      ##
      ## Recode the studies in each design
      ##
      for (j in seq_along(dat.design[[i]]$study))
        for (k in seq_len(k.d[i]))
          if (dat.design[[i]]$study[j] == studies.in.design[[i]][k])
            dat.design[[i]]$id.s[j] <- k
      ##
      ## Recode the treatments in each design
      ##
      for (j in seq_along(dat.design[[i]]$study))
        for (k in seq_len(n.d[i]))
          if (dat.design[[i]]$treat[j] == treat.per.design[[i]][k])
            dat.design[[i]]$id.t[j] <- k
      ##
      ## Calculate total patients per studies
      ##
      for (j in seq_along(dat.design[[i]]$studlab))
        dat.design[[i]]$n.study[j] <-
          with(dat.design[[i]], sum(n[which(study == study[j])]))
      ##
      ## Define c.xy and d.xy
      ##
      c.xy[[i]] <- array(rep(0, k.d[i] * n.d[i] * n.d[i]),
                         dim = c(n.d[i], n.d[i], k.d[i]))
      ##
      d.xy[[i]] <- array(rep(0, k.d[i] * n.d[i] * n.d[i]),
                         dim = c(n.d[i], n.d[i], k.d[i]))
      ##
      for (st in seq_len(k.d[i]))
        for (t1 in seq_len(n.d[i]))
          for (t2 in seq_len(n.d[i])) {
            c.xy[[i]][t1, t2, st] <-
              with(dat.design[[i]],
                   sum(    event[which(id.s == st & id.t == t1)]) *
                   sum(non.event[which(id.s == st & id.t == t2)]) /
                   n.study[id.s == st][1])
            ##
            d.xy[[i]][t1, t2, st] <-
              with(dat.design[[i]],
              (sum(    event[which(id.s == st & id.t == t1)]) +
               sum(non.event[which(id.s == st & id.t == t2)])) /
              n.study[id.s == st][1])
          }
      ##
      ## Define C.xy
      ##
      C.xy[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      ##
      for (j in seq_len(n.d[i]))
        for (k in seq_len(n.d[i]))
          C.xy[[i]][j, k] <- sum(c.xy[[i]][j, k, ])
      ##
      ## Define L.xy
      ##
      L.xy[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      for (j in seq_len(n.d[i]))
        for (k in seq_len(n.d[i]))
          L.xy[[i]][j, k] <- log(C.xy[[i]][j, k] / C.xy[[i]][k, j])
      ##
      ## Calculate the variance of L.xy (dimension n.d[i] x n.d[i])
      ##
      U.xyy[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      for (j in seq_len(n.d[i]))
        for (k in seq_len(n.d[i]))
          U.xyy[[i]][j, k] <-
            sum(c.xy[[i]][j, k, ] * d.xy[[i]][j, k, ]) /
            (2 * C.xy[[i]][j, k]^2) +
            sum(c.xy[[i]][j, k, ] * d.xy[[i]][k, j, ] + c.xy[[i]][k, j, ] *
                d.xy[[i]][j, k, ]) / (2 * C.xy[[i]][j, k] * C.xy[[i]][k, j]) +
            sum(c.xy[[i]][k, j, ] * d.xy[[i]][k, j, ]) / (2 * C.xy[[i]][k, j]^2)
      ##
      ## Calculate the covariance matrix U.xyz
      ##
      t.pl[[i]] = array(rep(0, k.d[i] * n.d[i] * n.d[i]),
                        dim = c(n.d[i], n.d[i], k.d[i]))
      ##
      for (j in seq_len(k.d[i]))
        t.pl[[i]][ , , j] <- with(dat.design[[i]], n.study[id.s == j][1])
      ##
      U.xyz[[i]] <- array(rep(0, n.d[i] * n.d[i] * n.d[i]),
                          dim = c(n.d[i], n.d[i], n.d[i]))
      ##
      ## Per study ...
      ##
      for (t1 in seq_len(n.d[i])) {
        for (t2 in seq_len(n.d[i])) {
          for (t3 in seq_len(n.d[i])) {
            for (st in seq_len(k.d[i])) {
              sel1 <- with(dat.design[[i]], which(id.s == st & id.t == t1))
              sel2 <- with(dat.design[[i]], which(id.s == st & id.t == t2))
              sel3 <- with(dat.design[[i]], which(id.s == st & id.t == t3))
              ##
              ps1[[i]][st] <-
                with(dat.design[[i]],
                     event[sel1] /
                     (t.pl[[i]][1, 1, st])^2 *
                     non.event[sel2] * non.event[sel3] *
                     (t1 != t2) * (t1 != t3) * (t2 != t3))
              ##
              ps2[[i]][st] <-
                with(dat.design[[i]],
                     n[sel1] / (t.pl[[i]][1, 1, st])^2 *
                     non.event[sel2] * event[sel3] *
                     (t1 != t2) * (t1 != t3) * (t2 != t3))
              ##
              ps3[[i]][st] <-
                with(dat.design[[i]],
                     n[sel1] / (t.pl[[i]][1, 1, st])^2 *
                     event[sel2] * non.event[sel3] *
                     (t1 != t2) * (t1 != t3) * (t2 != t3))
              ##
              ps4[[i]][st] <-
                with(dat.design[[i]],
                     non.event[sel1] / (t.pl[[i]][1, 1, st])^2 *
                     event[sel2] * event[sel3] *
                     (t1 != t2) * (t1 != t3) * (t2 != t3))
            }
            ##
            U.xyz[[i]][t1, t2, t3] <-
              sum(ps1[[i]][]) / (3 * C.xy[[i]][t1, t2] * C.xy[[i]][t1, t3]) +
              sum(ps2[[i]][]) / (3 * C.xy[[i]][t1, t2] * C.xy[[i]][t3, t1]) +
              sum(ps3[[i]][]) / (3 * C.xy[[i]][t2, t1] * C.xy[[i]][t1, t3]) +
              sum(ps4[[i]][]) / (3 * C.xy[[i]][t2, t1] * C.xy[[i]][t3, t1])
          }
        }
      }
      ##
      ## Calculate L.bar.xy
      ##
      L.bar.xy[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          L.bar.xy[[i]][t1, t2] <-
            (sum(L.xy[[i]][t1, ]) - sum(L.xy[[i]][t2, ])) / n.d[i]
      ##
      ## Calculate U.plus.xx
      ##
      for (t1 in seq_len(n.d[i]))
        U.xyy[[i]][t1, t1] <- 0
      ##
      U.plus.xx[[i]] <- rep(0, n.d[i])
      for (t1 in seq_len(n.d[i]))
        U.plus.xx[[i]][t1] <-
          sum(U.xyy[[i]][t1, seq_len(n.d[i])]) + sum(U.xyz[[i]][t1, , ])
      ##
      ## Calculate U.plus.xy
      ##
      U.new[[i]] <- array(rep(0, n.d[i] * n.d[i] * n.d[i]),
                          dim = c(n.d[i], n.d[i], n.d[i]))
      ##
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          for (t3 in seq_len(n.d[i]))
            U.new[[i]][t1, t2, t3] <-
              (t1 != t2) * (t1 != t3) * (t2 != t3) *
              U.xyz[[i]][t1, t2, t3] +
              (t1 != t2) * (t2 == t3) * U.xyy[[i]][t1, t2]
      ##
      U.plus.xy[[i]] = matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      ##
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          U.plus.xy[[i]][t1, t2] <-
            sum(U.new[[i]][seq_len(n.d[i]), t1, t2]) -
            sum(U.new[[i]][t1, t2, seq_len(n.d[i])]) -
            sum(U.new[[i]][t2, t1, seq_len(n.d[i])]) +
            U.new[[i]][t1, t2, t2]
      ##
      ## Variance of L.bar.xy
      ##
      Var.Lbar[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          Var.Lbar[[i]][t1, t2] <-
            (U.plus.xx[[i]][t1] - 2 * U.plus.xy[[i]][t1, t2] +
             U.plus.xx[[i]][t2]) / n.d[i]^2
      ##
      ## Covariance of L.bar.xy. Only a subset of covariances are
      ## calculate here, i.e. cov(L.bar_(1, t1), L.bar_(1, t2))
      ##
      CoVar.Lbar[[i]] <- matrix(rep(0, n.d[i]^2), nrow = n.d[i])
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          CoVar.Lbar[[i]][t1, t2] <-
            (U.plus.xx[[i]][1] - U.plus.xy[[i]][1, t2] -
             U.plus.xy[[i]][t1, 1] +
             U.plus.xy[[i]][t1, t2]) / n.d[i]^2
      ##
      ## Define matrix X
      ##
      for (j in seq_along(dat.design[[i]]$studlab))
        for (k in seq_along(trts))
          if (dat.design[[i]]$treat[j] == trts[k])
            dat.design[[i]]$id.t2[j] <- k
      ##
      ##
      ## Create y, the vector of treatment effects from each study. Only
      ## effects vs the first treatment are needed. y is coded as OR
      ## 1vsX
      ##
      for (j in seq_len(n.d[i] - 1))
        y[j + counter1] <- L.bar.xy[[i]][1, j + 1]
      ##
      counter1 <- counter1 + n.d[i] - 1
      ##
      ## Create V, the matrix of covariances of treatment effects from
      ## each study. Only covariances of effects vs the first treatment
      ## are needed.
      ##
      for (j in seq_len(n.d[i] - 1))
        for (k in seq_len(n.d[i] - 1))
          V1[j + counter2, k + counter2] <- (k == j) * Var.Lbar[[i]][1, j + 1]
      ##
      counter2 <- counter2 + n.d[i] - 1
      ##
      for (j in 2:(n.d[i]))
        for (k in 2:(n.d[i]))
          V2[j + counter3 - 1, k + counter3 - 1] <-
            (k != j) * CoVar.Lbar[[i]][j, k]
      ##
      counter3 <- counter3 + n.d[i] - 1
      ##
      for (j in seq_len(n.d[i] - 1)) {
        list1[N.j[i] + j, 1] <- dat.design[[i]]$id.t2[[1]]
        list1[N.j[i] + j, 2] <- dat.design[[i]]$id.t2[[j + 1]]
      }
      ##
      ## Drop unnecessary variables
      ##
      dat.design[[i]]$id.s <- NULL
      dat.design[[i]]$id.t <- NULL
      dat.design[[i]]$id.t2 <- NULL
      dat.design[[i]]$n.study <- NULL
    }
    ##
    V <- V1 + V2
    ##
    basic.contrasts <- 2:n.treat
    ##
    for (i in 1:length.y)
      for (k in 1:(n.treat - 1)) {
        if (list1[i, 1] == basic.contrasts[k])
          X[i, k] = -1
        if (list1[i, 2] == basic.contrasts[k])
          X[i, k] = 1
      }
    ##
    W <- solve(V)
    ##
    TE.basic <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
  }
  else if (method == "NCH") {
    ##
    ## NCH method
    ##
    dat.long$ttt <- 0
    ##
    for (k in 1:n.treat)
      for (i in seq_along(dat.long$studlab))
        if (dat.long$treat[i] == trts[k])
          dat.long$ttt[i] <- k
    ##
    dat.long$treat <- as.numeric(dat.long$ttt)
    dat.long$ttt <- NULL
    ##
    dat.long$treat <- as.numeric(dat.long$treat)
    ##
    dat.long <- dat.long[order(dat.long$studlab, dat.long$treat), ] # necessary ???
    ##
    for (i in unique(dat.long$studlab)) {
      counter <- 1
      for (j in seq_along(dat.long$studlab))
        if (dat.long$studlab[j] == i) {
          dat.long$count[j] <- counter
          counter <- counter + 1
        }
    }
    ##
    max.arms <- max(dat.long$count)
    ##
    d1 <- reshape(dat.long, idvar = "studlab",
                  timevar = "count", direction = "wide")
    ##
    d1 <- d1[order(d1$studlab), ]
    ##
    d1$N1 <- 0 # d1$n.all = 0
    d1$narms <- 0
    ##
    for (i in seq_len(max.arms)) {
      ev <- paste0("event.", i)
      tot <- paste0("n.", i)
      d1$N1 <- rowSums(cbind(d1$N1, d1[, colnames(d1) == ev]),
                       na.rm = TRUE)
      ## d1$n.all <- rowSums(cbind(d1$n.all, d1[, colnames(d1) == tot]), na.rm = TRUE)
    }
    ##
    for (i in seq_along(d1$studlab))
      for (k in 1:max.arms) {
        ev <- paste0("event.", k)
        if (!is.na(d1[, colnames(d1) == ev][i]))
          d1$narms[i] <- k
      }
    ##
    for (i in seq_along(colnames(d1)))
      for (k in seq_along(colnames(d1))) {
        if (colnames(d1)[k] == paste0("treat.", i))
          colnames(d1)[k] = paste0("t", i)
        ##
        if (colnames(d1)[k] == paste0("event.", i))
          colnames(d1)[k] = paste0("r", i)
        ##
        if (colnames(d1)[k] == paste0("n.", i))
          colnames(d1)[k] = paste0("n", i)
      }
    ##
    ## Likelihood function
    ##
    myLik1 <- function(mypar) {
      x <- mypar
      ##
      myLogLik1 <- myLogLik2 <- myLogLik3 <- rep(0, length(d1$studlab))
      ##
      for (i in seq_along(d1$studlab)) {
        for (k in 2:d1$narms[i]) {
          myLogLik1[i] <-
            myLogLik1[i] +
            d1[, colnames(d1) == paste0("r", k)][i] *
            (x[d1[, colnames(d1) == paste0("t", k)][i]] -
             x[d1$t1[i]] * (d1$t1[i] != 1))
          ##
          myLogLik2[i] <-
            myLogLik2[i] +
            d1[, colnames(d1) == paste0("n", k)][i] *
            exp(x[d1[, colnames(d1) == paste0("t", k)][i]] -
                x[d1$t1[i]] * (d1$t1[i] != 1))
          ##
          myLogLik3[i] <- -d1$N1[i] * log(d1$n1[i] + myLogLik2[i])
        }
      }
      ##
      myLogLik <- sum(myLogLik1 + myLogLik3)
      ##
      myLogLik
    }
    ##
    opt <- optim(rep(0, n.treat), myLik1, method = "L-BFGS-B",
                 lower = -Inf, upper = Inf,
                 control = list(fnscale = -1, maxit = 10000),
                 hessian = TRUE)
    ##
    W <- solve(-opt$hessian[2:n.treat, 2:n.treat])
    ##
    TE.basic <- -(opt$par)[2:n.treat]
  }
  ##
  ## Define H matrix
  ##
  H <- matrix(0,
              nrow = n.treat * ((n.treat - 1) / 2),
              ncol = n.treat - 1)
  ##
  diag(H) <- 1
  ##
  if (n.treat > 2) {
    t1 <- c()
    t2 <- c()
    for (i in 2:(n.treat - 1))
      for (j in (i + 1):(n.treat)) {
        t1 <- rbind(t1, i)
        t2 <- rbind(t2, j)
      }
    ##
    h1 <- matrix(c(t1, t2), nrow = nrow(t1))
    ##
    for (i in 1:((n.treat - 1) * (n.treat - 2) / 2))
      for (j in 1:(n.treat - 1))
        H[i + n.treat - 1, j] <- -(h1[i, 1] == j + 1) + (h1[i, 2] == j + 1)
  }
  ##
  ## Fixed effects matrix
  ##
  d.hat <- H %*% TE.basic
  ##
  TE.fixed[lower.tri(TE.fixed, diag = FALSE)] <-  d.hat
  TE.fixed <- t(TE.fixed)
  TE.fixed[lower.tri(TE.fixed, diag = FALSE)] <- -d.hat
  diag(TE.fixed) <- 0
  ##
  ## Matrix with standard errors
  ##
  if (method == "MH")
    cov.d.hat <- H %*% solve(t(X) %*% W %*% X) %*% t(H)
  else
    cov.d.hat <- H %*% W %*% t(H)
  ##
  seTE.fixed[lower.tri(seTE.fixed, diag = FALSE)] <- sqrt(diag(cov.d.hat))
  seTE.fixed <- t(seTE.fixed)
  seTE.fixed[lower.tri(seTE.fixed, diag = FALSE)] <- sqrt(diag(cov.d.hat))
  diag(seTE.fixed) <- 0
  ##
  ## Confidence intervals
  ##
  ci.f <- ci(TE.fixed, seTE.fixed, level = level.comb)
  ##
  ## Inconsistency global
  ##
  if (method == "MH") {
    Q <- as.vector(t(y - X %*% TE.basic) %*% solve(V) %*%
                    (y - X %*% TE.basic))
    df.Q <- sum(n.d - 1) - (n.treat - 1)
    pval.Q <- pvalQ(Q, df.Q)
  }
  else {
    Q <- NA
    df.Q <- NA
    pval.Q <- NA
  }
  
  
  ##
  ##
  ## (9) Inconsistency evaluation: direct MH estimates
  ##
  ##
  B.matrix <- createB(treat1.pos, treat2.pos)
  B.matrix <- unique(B.matrix)
  colnames(B.matrix) <- trts
  ##
  for (i in seq_len(nrow(B.matrix))) {
    sel.treat1 <- colnames(B.matrix)[B.matrix[i, ] ==  1]
    sel.treat2 <- colnames(B.matrix)[B.matrix[i, ] == -1]
    selstud <- treat1 == sel.treat1 & treat2 == sel.treat2
    ##
    m.i <- metabin(event1, n1, event2, n2, subset = selstud,
                   sm = "OR",
                   incr = incr, allincr = allincr, addincr = addincr,
                   allstudies = allstudies, MH.exact = !cc.pooled)
    ##
    TE.i   <- m.i$TE.fixed
    seTE.i <- m.i$seTE.fixed
    ##
    TE.direct.fixed[sel.treat1, sel.treat2]   <- TE.i
    seTE.direct.fixed[sel.treat1, sel.treat2] <- seTE.i
    ##
    TE.direct.fixed[sel.treat2, sel.treat1]   <- -TE.i
    seTE.direct.fixed[sel.treat2, sel.treat1] <- seTE.i
  }
  ##
  rm(sel.treat1, sel.treat2, selstud, m.i, TE.i, seTE.i)
  ##
  ci.d <- meta::ci(TE.direct.fixed, seTE.direct.fixed, level = level.comb)
  
  
  labels <- sort(unique(c(treat1, treat2)))
  ##
  if (!is.null(seq))
    seq <- setseq(seq, labels)
  else {
    seq <- labels
    if (is.numeric(seq))
      seq <- as.character(seq)
  }
  
  
  res <- list(studlab = studlab,
              treat1 = treat1,
              treat2 = treat2,
              ##
              TE = data$TE[!data$.drop],
              seTE = data$seTE[!data$.drop],
              seTE.adj = rep(NA, sum(!data$.drop)),
              ##
              event1 = event1,
              event2 = event2,
              n1 = n1,
              n2 = n2,
              ##
              k = n.studies,
              m = m,
              n = length(trts),
              d = d,
              ##
              trts = trts,
              k.trts = rowSums(A),
              n.trts = NA,
              events.trts = NA,
              ##
              studies = NA,
              narms = NA,
              ##
              designs = designs,
              ##
              TE.fixed = TE.fixed,
              seTE.fixed = seTE.fixed,
              lower.fixed = ci.f$lower,
              upper.fixed = ci.f$upper,
              zval.fixed = ci.f$z,
              pval.fixed = ci.f$p,
              ##
              TE.random = NA,
              seTE.random = NA,
              lower.random = NA,
              upper.random = NA,
              zval.random = NA,
              pval.random = NA,
              ##
              prop.direct.fixed = NA,
              prop.direct.random = NA,
              ##
              TE.direct.fixed = TE.direct.fixed,
              seTE.direct.fixed = seTE.direct.fixed,
              lower.direct.fixed = ci.d$lower,
              upper.direct.fixed = ci.d$upper,
              zval.direct.fixed = ci.d$z,
              pval.direct.fixed = ci.d$p,
              ##
              TE.indirect.fixed = NA,
              seTE.indirect.fixed = NA,
              lower.indirect.fixed = NA,
              upper.indirect.fixed = NA,
              zval.indirect.fixed = NA,
              pval.indirect.fixed = NA,
              ##
              Q = Q,
              df.Q = df.Q,
              pval.Q = pval.Q,
              I2 = NA,
              tau = NA,
              ##
              Q.heterogeneity = NA,
              df.Q.heterogeneity = NA,
              pval.Q.heterogeneity = NA,
              Q.inconsistency = NA,
              df.Q.inconsistency = NA,
              pval.Q.inconsistency = NA,
              ##
              Q.decomp = NA,
              ##
              A.matrix = A,
              H.matrix = H,
              ##
              n.matrix = NA,
              events.matrix = NA,
              ##
              P.fixed = NA,
              P.random = NA,
              ##
              Cov.fixed = NA,
              Cov.random = NA,
              ##
              sm = sm,
              method = method,
              ##
              incr = incr,
              allincr = allincr,
              addincr = addincr,
              allstudies = allstudies,
              cc.pooled = cc.pooled,
              ##
              level = level,
              level.comb = level.comb,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              ##
              prediction = prediction,
              level.predict = level.predict,
              ##
              reference.group = reference.group,
              baseline.reference = baseline.reference,
              all.treatments = all.treatments,
              seq = seq,
              ##
              tau.preset = tau.preset,
              ##
              tol.multiarm = tol.multiarm,
              details.chkmultiarm = details.chkmultiarm,
              ##
              sep.trts = sep.trts,
              nchar.trts = nchar.trts,
              ##
              backtransf = backtransf,
              ##
              title = title,
              ##
              data = data,
              data.wide = dat.wide,
              data.long = dat.long,
              data.design = dat.design,
              ##
              y = if (method == "MH") y else NULL,
              V = if (method == "MH") V else NULL,
              X = if (method == "MH") X else NULL,
              ##
              warn = warn,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- c("netmetabin", "netmeta")
  
  
  ##
  ## Study overview
  ##
  p0 <- prepare(rep(1, nrow(res$data.wide)),
                rep(1, nrow(res$data.wide)),
                res$data.wide$treat1,
                res$data.wide$treat2,
                res$data.wide$studlab)
  tdata <- data.frame(studies = p0$studlab, narms = p0$narms,
                      stringsAsFactors = FALSE)
  tdata <- unique(tdata[order(tdata$studies, tdata$narms), ])
  res$studies <- tdata$studies
  res$narms <- tdata$narms
  ##
  res$data <- merge(res$data,
                    data.frame(.studlab = res$studies,
                               .narms = res$narms),
                    by = ".studlab", all.x = TRUE)
  
  
  res$data <- res$data[order(res$data$.order), ]
  res$data$.order <- NULL
  rownames(res$data) <- seq_len(nrow(res$data))
  ##
  res$data.wide <- res$data.wide[order(res$data.wide$.order), ]
  res$data.wide$.order <- NULL
  rownames(res$data.wide) <- seq_len(nrow(res$data.wide))
  ##
  res$data.long <- res$data.long[order(res$data.long$.order), ]
  res$data.long$.order <- NULL
  rownames(res$data.long) <- seq_len(nrow(res$data.long))
  
  
  res$events.matrix <- netmatrix(res, event1 + event2, func = "sum")
  ##
  dat.e <- bySummary(c(event1, event2), c(treat1, treat2), long = FALSE)
  rownames(dat.e) <- dat.e$indices
  res$events.trts <- dat.e[trts, "sum"]
  names(res$events.trts) <- trts
  ##
  res$n.matrix <- netmatrix(res, n1 + n2, func = "sum")
  ##
  dat.n <- bySummary(c(n1, n2), c(treat1, treat2), long = FALSE)
  rownames(dat.n) <- dat.n$indices
  res$n.trts <- dat.n[trts, "sum"]
  names(res$n.trts) <- trts
  
  
  res
}
