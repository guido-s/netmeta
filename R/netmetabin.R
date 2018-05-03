netmetabin <- function(event1, n1, event2, n2,
                       treat1, treat2, studlab,
                       data = NULL, subset = NULL,
                       sm,
                       method = "MH",
                       incr = gs("incr"),
                       allincr = gs("allincr"), addincr = gs("addincr"),
                       allstudies = gs("allstudies"),
                       MH.exact = TRUE,
                       level = 0.95, level.comb = 0.95,
                       comb.fixed = gs("comb.fixed"),
                       comb.random = method == "Inverse" &
                         (gs("comb.random") | !is.null(tau.preset)),
                       ##
                       prediction = FALSE,
                       level.predict = 0.95,
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
  chklogical(MH.exact)
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
            "method and non-central hypergeometric distribution.")
    comb.fixed <- TRUE
  }
  ##
  if (method != "Inverse" & comb.random) {
    warning("Argument 'comb.random' set to FALSE for Mantel-Haenszel ",
            "method and non-central hypergeometric distribution.")
    comb.random <- FALSE
  }
  ##
  if (method != "Inverse" & prediction) {
    warning("Argument 'prediction' set to FALSE for Mantel-Haenszel ",
            "method and non-central hypergeometric distribution.")
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
  else
    data$...order <- seq_len(nrow(data))
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
        warning("Argument 'sm' set to 'OR'.")
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
        warning("Argument 'sm' set to 'OR'.")
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
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  missing.subset <- is.null(subset)
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
      data$...order <- seq_len(nrow(data))
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
  }
  ##
  labels <- sort(unique(c(treat1, treat2)))
  ##
  if (compmatch(labels, sep.trts)) {
    if (!missing.sep.trts)
      warning("Separator '", sep.trts, "' used in at least one treatment label. Try to use predefined separators: ':', '-', '_', '/', '+', '.', '|', '*'.")
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
      stop("All predefined separators (':', '-', '_', '/', '+', '.', '|', '*') are used in at least one treatment label.",
           "\n   Please specify a different character that should be used as separator (argument 'sep.trts').",
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
              "Assuming that comparisons are from independent studies.")
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
    stop(paste("Study '", names(tabnarms)[sel.narms],
               "' has a wrong number of comparisons.",
               "\n  Please provide data for all treatment comparisons",
               " (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
               sep = ""),
         call. = FALSE)
  if (sum(sel.narms) > 1)
    stop(paste("The following studies have a wrong number of comparisons: ",
               paste(paste("'", names(tabnarms)[sel.narms], "'", sep = ""),
                     collapse = ", "),
               "\n  Please provide data for all treatment comparisons",
               " (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
               sep = ""),
         call. = FALSE)
  ##
  ## Check number of subgraphs
  ##
  n.subnets <- netconnection(treat1, treat2, studlab)$n.subnets
  ##
  if (n.subnets > 1)
    stop(paste("Network consists of ", n.subnets, " separate sub-networks.\n  ",
               "Use R function 'netconnection' to identify sub-networks.",
               sep = ""),
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
            "hypergeometric method. Set to zero.")
    incr <- 0
  }
  ##
  lengthunique <- function(x) length(unique(x))
  
  
  ##
  ##
  ## (5) Create analysis dataset
  ##
  ##
  dat <- rbind(data.frame(studlab, treat = treat1, event = event1, n = n1,
                          stringsAsFactors = FALSE),
               data.frame(studlab, treat = treat2, event = event2, n = n2,
                          stringsAsFactors = FALSE))
  dat$...order <- seq_len(nrow(dat))
  ##
  d1 <- d2 <- data.frame(studlab = studlab, treat1 = treat1, treat2 = treat2,
                         event1 = event1, n1 = n1, event2 = event2, n2 = n2,
                         stringsAsFactors = FALSE)
  d1$.first <- TRUE
  d2$.first <- FALSE
  dat1 <- rbind(d1, d2)
  dat1$...order <- seq_len(nrow(dat1))
  rm(d1, d2)
  ##
  data$.drop <- rep(FALSE, nrow(data))
  ##
  ## Add variable 'non.event'
  ##
  dat$non.event <- dat$n - dat$event
  ##
  dat1$non.event1 <- dat1$n1 - dat1$event1
  dat1$non.event2 <- dat1$n2 - dat1$event2
  ##
  ## Pool events and totals in studies with multiple arms of the same
  ## treat
  ##
  dupl <- duplicated(dat[, c("studlab", "treat")])
  ##
  if (any(dupl)) {
    dat.D <- dat[ dupl, ]
    dat   <- dat[!dupl, ]
    ##
    for (i in seq_len(nrow(dat.D))) {
      sel <- (dat$studlab == dat.D$studlab[i] & dat$treat == dat.D$treat[i])
      ##
      dat$event[sel] <- dat$event[sel] + dat.D$event[i]
      dat$n[sel] <- dat$n[sel] + dat.D$n[i]
      ##
      sel1 <- dat1$studlab == dat.D$studlab[i] & dat1$treat1 == dat.D$treat[i]
      sel2 <- dat1$studlab == dat.D$studlab[i] & dat1$treat2 == dat.D$treat[i]
      ##
      dat1$event1[sel1] <- dat1$event1[sel1] + dat.D$event[i]
      dat1$event2[sel2] <- dat1$event2[sel2] + dat.D$event[i]
      dat1$n1[sel1] <- dat1$n1[sel1] + dat.D$n[i]
      dat1$n2[sel2] <- dat1$n2[sel2] + dat.D$n[i]
    }
    ##
    rm(dat.D, sel, sel1, sel2)
  }
  ##
  rm(dupl)
  ##
  ## Remove single arms studies from dataset
  ##
  sel.first <- !duplicated(dat$studlab)
  sel.other <-  duplicated(dat$studlab)
  ##
  first <- unique(dat$studlab[sel.first])
  other <- unique(dat$studlab[sel.other])
  ##
  if (length(first) != length(other)) {
    single <- !(first %in% other)
    ##
    if (warn)
      if (sum(single) == 1)
        warning("Single-arm study '", first[single],
                "' excluded from network meta-analysis.")
      else
        warning("Single-arm studies excluded from network meta-analysis: ",
                paste(paste0("'", first[single], "'"),
                      collapse = ", "))
    ##
    dat  <-  dat[!(dat$studlab  %in% first[single]), ]
    dat1 <- dat1[!(dat1$studlab %in% first[single]), ]
    ##
    data$.drop <- data$.drop | data$studlab %in% first[single]
    ##
    rm(single)
  }
  ##
  rm(sel.first, sel.other, first, other)
  ##
  ## Add variable 'design' with study design
  ##
  dat.design <- nma.krahn(netmeta(pairwise(studlab = dat$studlab,
                                           treat = dat$treat,
                                           event = dat$event,
                                           n = dat$n,
                                           sm = "RD"))
                          )$studies[, c("studlab", "design")]
  ##
  dat  <- merge(dat,  unique(dat.design), by = "studlab")
  dat1 <- merge(dat1, unique(dat.design), by = "studlab")
  ##
  dat  <-  dat[order(dat$...order), ]
  dat1 <- dat1[order(dat1$...order), ]
  ##
  names(dat.design) <- c("studlab", ".design")
  data <- merge(data, unique(dat.design), by = "studlab")
  data <- data[order(data$...order), ]
  ##
  rm(dat.design)
  
  
  ##
  ##
  ## (6) Stage I: setting up the data (Efthimiou et al., 2018)
  ##
  ##
  ## Step i. Remove all-zero studies
  ##
  n.events <- with(dat, tapply(event, studlab, sum))
  ##
  if (any(n.events == 0) & !allstudies) {
    allzero <- n.events == 0
    ##
    if (warn)
      if (sum(allzero) == 1)
        warning("Study '", names(n.events)[allzero],
                "' without any events excluded from network meta-analysis.")
      else
        warning("Studies without any events excluded ",
                "from network meta-analysis: ",
                paste(paste0("'", names(n.events)[allzero], "'"),
                      collapse = ", "))
    ##
    dat  <-  dat[dat$studlab %in% names(n.events)[!allzero], , drop = FALSE]
    dat1 <- dat1[dat$studlab %in% names(n.events)[!allzero], , drop = FALSE]
    ##
    data$.drop <- data$.drop | data$studlab %in% names(n.events)[allzero]
    ##
    rm(allzero)
  }
  ##
  rm(n.events)
  ##
  ## Add variable 'study' with study numbers
  ##
  dat$study  <- as.numeric(factor(dat$studlab, levels = unique(dat$studlab)))
  dat1$study <- as.numeric(factor(dat1$studlab, levels = unique(dat1$studlab)))
  ##
  ## Sort by study and treatment
  ##
  o <- order(dat$study, dat$treat)
  dat  <-  dat[o, ]
  dat1 <- dat1[o, ]
  ##
  ## Step iii. Drop treatment arms without events from individual
  ##           designs (argument 'MH.exact' is TRUE) or add increment
  ##           if argument 'MH.exact' is FALSE and argument 'incr' is
  ##           larger than zero
  ##
  if (method == "MH") {
    ##
    d.events <- with(dat, tapply(event, list(design, treat), sum))
    ##
    if (any(d.events == 0, na.rm = TRUE) | addincr) {
      zero <- d.events == 0
      zerocells <- as.data.frame(which(zero, arr.ind = TRUE),
                                 stringsAsFactors = FALSE)
      ##
      zerocells$design <- rownames(zero)[zerocells$row]
      zerocells$treat <- colnames(zero)[zerocells$col]
      ##
      if (!allincr & !addincr) {
        sel  <- rep(0, nrow(dat))
        sel1 <- rep(0, nrow(dat1))
        seld <- rep(0, nrow(data))
        ##
        for (i in seq_along(zerocells$design)) {
          sel <- sel + (dat$design == zerocells$design[i] &
                        dat$treat == zerocells$treat[i])
          ##
          sel1 <- sel1 + (dat1$design == zerocells$design[i] &
                          (dat1$treat1 == zerocells$treat[i] |
                           dat1$treat2 == zerocells$treat[i]))
          ##
          seld <- seld + (data$.design == zerocells$design[i] &
                          (data$.treat1 == zerocells$treat[i] |
                           data$.treat2 == zerocells$treat[i]))
        }
        ##
        sel  <- sel > 0
        sel1 <- sel1 > 0
        seld <- seld > 0
      }
      else {
        sel  <- rep(TRUE, length(dat$event))
        sel1 <- rep(TRUE, length(dat1$event1))
        seld <- rep(TRUE, length(data$.event1))
      }
      ##
      if (MH.exact) {
        if (warn)
          if (sum(zero, na.rm = TRUE) == 1)
            warning("Treatment arm '", zerocells$treat,
                    "' without events in design '",
                    zerocells$design, "' excluded from network meta-analysis.")
          else
            warning("Treatment arms without events in a design excluded ",
                    "from network meta-analysis:\n    ",
                    paste0("'",
                           paste0(paste0(zerocells$treat, " in "),
                                  zerocells$design),
                           "'",
                           collapse = ", "))
        ##
        dat  <-  dat[!sel, , drop = FALSE]
        dat1 <- dat1[!sel1, , drop = FALSE]
        ##
        data$.drop <- data$.drop | seld
        ##
        dat.incr <- dat
        ##
        rm(zero, zerocells, sel)
      }
      else {
        dat.incr <- dat
        dat.incr$event[sel] <- dat.incr$event[sel] + incr
        dat.incr$non.event[sel] <- dat.incr$non.event[sel] + incr
        dat.incr$n[sel] <- dat.incr$n[sel] + 2 * incr
        ##
        dat1$incr <- 0
        dat1$incr[sel] <- incr
        ##
        rm(sel)
      }
    }
    else
      dat.incr <- dat
    ##
    rm(d.events)
  }
  else
    dat.incr <- dat
  ##
  dat.incr$design <- as.character(dat.incr$design)
  ##
  ## Step iv. Remove designs with single treatment arm from dataset
  ##
  d.single <- with(dat, tapply(treat, design, lengthunique))
  ##
  if (any(d.single == 1, na.rm = TRUE)) {
    single <- !is.na(d.single) & d.single == 1
    design.single <- names(d.single)[single]
    ##
    if (warn)
      if (sum(single) == 1)
        warning("Design '", design.single,
                "' with single treatment arm excluded ",
                "from network meta-analysis.")
      else
        warning("Designs with single treatment arm excluded ",
                "from network meta-analysis: ",
                paste(paste0("'", design.single, "'"),
                      collapse = ", "))
    ##
    dat <- dat[!(dat$design %in% design.single), , drop = FALSE]
    dat1 <- dat1[!(dat1$design %in% design.single), , drop = FALSE]
    ##
    data$.drop <- data$.drop | data$.design %in% design.single
    ##
    rm(single, design.single)
  }
  ##
  rm(d.single)
  ##
  dat <-  dat[order(dat$design, dat$studlab, dat$treat), ]
  dat1 <- dat1[order(dat1$design, dat1$studlab, dat1$treat1, dat1$treat2), ]
  dat1 <- dat1[dat1$.first, ]
  ##
  dat  <-  dat[, c("studlab", "treat", "event", "non.event", "n",
                   "design", "study", "...order")]
  dat1 <- dat1[, c("studlab", "treat1", "treat2",
                   "event1", "event2", "non.event1", "non.event2",
                   "n1", "n2",
                   "design", "study", "...order")]
  ##
  dat$studlab <- as.character(dat$studlab)
  dat$design <- as.character(dat$design)
  dat$treat <- as.character(dat$treat)
  ##
  dat1$studlab <- as.character(dat1$studlab)
  dat1$design <- as.character(dat1$design)
  dat1$treat1 <- as.character(dat1$treat1)
  dat1$treat2 <- as.character(dat1$treat2)
  ##
  o <- order(dat1$...order)
  studlab <- dat1$studlab[o]
  treat1 <- dat1$treat1[o]
  treat2 <- dat1$treat2[o]
  ##
  event1 <- dat1$event1[o]
  event2 <- dat1$event2[o]
  n1 <- dat1$n1[o]
  n2 <- dat1$n2[o]
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
    p.iv <- pairwise(studlab = dat$studlab,
                     treat = dat$treat,
                     event = dat$event,
                     n = dat$n,
                     sm = sm,
                     incr = incr,
                     allincr = allincr, addincr = addincr,
                     allstudies = allstudies)
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
  n.d <- tapply(dat$treat,   dat$design, lengthunique)
  k.d <- tapply(dat$studlab, dat$design, lengthunique)
  ##
  designs <- names(k.d)
  d <- length(designs)
  seq.d <- seq_len(d)
  ##
  ddat <- vector("list", d)
  names(ddat) <- designs
  ##
  for (i in seq.d) 
    ddat[[i]] <- dat.incr[dat.incr$design == designs[i], , drop = FALSE]
  ##
  if (method == "MH") {
    ##
    ## MH method
    ##
    treat.per.design <- vector("list", d)
    ##
    for (i in seq.d)
      treat.per.design[[i]] <- unique(ddat[[i]]$treat)
    ##
    ## List of studies in each design
    ##
    studies.in.design <- vector("list", d)
    ##
    for (i in seq.d)
      studies.in.design[[i]] <- unique(ddat[[i]]$study)
    ##
    ## Recode the studies in each design
    ##
    for (i in seq.d)
      for (j in seq_along(ddat[[i]]$study))
        for (k in seq_len(k.d[i]))
          if (ddat[[i]]$study[j] == studies.in.design[[i]][k])
            ddat[[i]]$id.s[j] <- k
    ##
    ## Recode the treatments in each design
    ##
    for (i in seq.d)
      for (j in seq_along(ddat[[i]]$study))
        for (k in seq_len(n.d[i]))
          if (ddat[[i]]$treat[j] == treat.per.design[[i]][k])
            ddat[[i]]$id.t[j] <- k
    ##
    ## Calculate total patients per studies
    ##
    for (i in seq.d)
      for (j in seq_along(ddat[[i]]$studlab))
        ddat[[i]]$n.study[j] <-
          sum(ddat[[i]]$n[which(ddat[[i]]$study == ddat[[i]]$study[j])])
    
    
    ##
    ## Define c.xy and d.xy
    ##
    c.xy <- vector("list", d)
    d.xy <- vector("list", d)
    ##
    for (i in seq.d)
      c.xy[[i]] <- array(rep(0, k.d[i] * n.d[i] * n.d[i]),
                         dim = c(n.d[i], n.d[i], k.d[i]))
    ##
    for (i in seq.d)
      d.xy[[i]] <- array(rep(0, k.d[i] * n.d[i] * n.d[i]),
                         dim = c(n.d[i], n.d[i], k.d[i]))
    ##
    for (i in seq.d)
      for (st in seq_len(k.d[i]))
        for (t1 in seq_len(n.d[i]))
          for (t2 in seq_len(n.d[i])) {
            c.xy[[i]][t1, t2, st] <-
              sum(ddat[[i]]$event[which(ddat[[i]]$id.s == st &
                                        ddat[[i]]$id.t == t1)]) *
              sum(ddat[[i]]$non.event[which(ddat[[i]]$id.s == st &
                                            ddat[[i]]$id.t == t2)]) /
              ddat[[i]]$n.study[ddat[[i]]$id.s == st][1]
            ##
            d.xy[[i]][t1, t2, st] <-
              (sum(ddat[[i]]$event[which(ddat[[i]]$id.s == st &
                                         ddat[[i]]$id.t == t1)]) + 
               sum(ddat[[i]]$non.event[which(ddat[[i]]$id.s == st &
                                             ddat[[i]]$id.t == t2)])) /
              ddat[[i]]$n.study[ddat[[i]]$id.s == st][1]
          }
    ##
    ## Define C.xy
    ##
    C.xy <- vector("list", d)
    ##
    for (i in seq.d) {
      C.xy[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      for (j in seq_len(n.d[i]))
        for (k in seq_len(n.d[i]))
          C.xy[[i]][j, k] <- sum(c.xy[[i]][j, k, ])
    }
    ##
    ## Define L.xy
    ##
    L.xy <- vector("list", d)
    ##
    for (i in seq.d) {
      L.xy[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      for (j in seq_len(n.d[i]))
        for (k in seq_len(n.d[i]))
          L.xy[[i]][j, k] <- log(C.xy[[i]][j, k] / C.xy[[i]][k, j])
    }
    ##
    ## Calculate the variance of L.xy (dimension n.d[i] x n.d[i])
    ##
    U.xyy <- vector("list", d)
    for (i in seq.d) {
      U.xyy[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      for (j in seq_len(n.d[i]))
        for (k in seq_len(n.d[i]))
          U.xyy[[i]][j, k] <-
            sum(c.xy[[i]][j, k, ] * d.xy[[i]][j, k, ]) /
            (2 * C.xy[[i]][j, k]^2) +
            sum(c.xy[[i]][j, k, ] * d.xy[[i]][k, j, ] + c.xy[[i]][k, j, ] *
                d.xy[[i]][j, k, ]) / (2 * C.xy[[i]][j, k] * C.xy[[i]][k, j]) +
            sum(c.xy[[i]][k, j, ] * d.xy[[i]][k, j, ]) / (2 * C.xy[[i]][k, j]^2)
    }
    ##
    ## Calculate the covariance matrix U.xyz
    ##
    t.pl <- vector("list", d)
    ##
    for (i in seq.d)
      t.pl[[i]] = array(rep(0, k.d[i] * n.d[i] * n.d[i]),
                        dim = c(n.d[i], n.d[i], k.d[i]))
    ##
    for (i in seq.d)
      for (j in seq_len(k.d[i]))
        t.pl[[i]][,, j] <- ddat[[i]]$n.study[ddat[[i]]$id.s == j][1]
    ##    
    U.xyz <- vector("list", d)
    ##
    for (i in seq.d)
      U.xyz[[i]] <- array(rep(0, n.d[i] * n.d[i] * n.d[i]),
                          dim = c(n.d[i], n.d[i], n.d[i]))
    ##
    ## Per study ...
    ##
    ps1 <- ps2 <- ps3 <- ps4 <- vector("list", d)
    ##
    for (i in seq.d) {
      ps1[[i]] <- ps2[[i]] <-
        ps3[[i]] <- ps4[[i]] <- rep(0, k.d[i])
      ##
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          for (t3 in seq_len(n.d[i]))
            for (st in seq_len(k.d[i])) {
              sel1 <- which(ddat[[i]]$id.s == st & ddat[[i]]$id.t == t1)
              sel2 <- which(ddat[[i]]$id.s == st & ddat[[i]]$id.t == t2)
              sel3 <- which(ddat[[i]]$id.s == st & ddat[[i]]$id.t == t3)
              ##
              ps1[[i]][st] <-
                ddat[[i]]$event[sel1] /
                (t.pl[[i]][1, 1, st])^2 *
                                        ddat[[i]]$non.event[sel2] * 
                                        ddat[[i]]$non.event[sel3] * 
                                        (t1 != t2) * (t1 != t3) * (t2 != t3)
              ##
              ps2[[i]][st] <- ddat[[i]]$n[sel1] /
                (t.pl[[i]][1, 1, st])^2 * 
                                        ddat[[i]]$non.event[sel2] * 
                                        ddat[[i]]$event[sel3] *
                                        (t1 != t2) * (t1 != t3) * (t2 != t3)
              ##
              ps3[[i]][st] <- ddat[[i]]$n[sel1] /
                (t.pl[[i]][1, 1, st])^2 * 
                                        ddat[[i]]$event[sel2] * 
                                        ddat[[i]]$non.event[sel3] * 
                                        (t1 != t2) * (t1 != t3) * (t2 != t3)
              ##
              ps4[[i]][st] <- ddat[[i]]$non.event[sel1] /
                (t.pl[[i]][1, 1, st])^2 * 
                                        ddat[[i]]$event[sel2] * 
                                        ddat[[i]]$event[sel3] * 
                                        (t1 != t2) * (t1 != t3) * (t2 != t3)
            }
            ##
            U.xyz[[i]][t1, t2, t3] <- sum(ps1[[i]][]) / (3 * C.xy[[i]][t1, t2] * C.xy[[i]][t1, t3]) + 
              sum(ps2[[i]][]) / (3 * C.xy[[i]][t1, t2] * C.xy[[i]][t3, t1]) + 
              sum(ps3[[i]][]) / (3 * C.xy[[i]][t2, t1] * C.xy[[i]][t1, t3]) + 
              sum(ps4[[i]][]) / (3 * C.xy[[i]][t2, t1] * C.xy[[i]][t3, t1])
          }
    ##
    ## Calculate L.bar.xy
    ##
    L.bar.xy <- vector("list", d)
    ##
    for (i in seq.d) {
      L.bar.xy[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          L.bar.xy[[i]][t1, t2] <-
            (sum(L.xy[[i]][t1, ]) - sum(L.xy[[i]][t2, ])) / n.d[i]
    }
    ##
    ## Calculate U.plus.xx
    ##
    for (i in seq.d)
      for (t1 in seq_len(n.d[i]))
        U.xyy[[i]][t1, t1] <- 0
    ##
    U.plus.xx <- vector("list", d)
    ##
    for (i in seq.d) {
      U.plus.xx[[i]] <- rep(0, n.d[i])
      for (t1 in seq_len(n.d[i]))
        U.plus.xx[[i]][t1] <-
          sum(U.xyy[[i]][t1, seq_len(n.d[i])]) + sum(U.xyz[[i]][t1,, ])
    }
    ##
    ## Calculate U.plus.xy
    ##
    U.new <- vector("list", d)
    ##
    for (i in seq.d)
      U.new[[i]] <- array(rep(0, n.d[i] * n.d[i] * n.d[i]),
                          dim = c(n.d[i], n.d[i], n.d[i]))
    ##
    for (i in seq.d)
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          for (t3 in seq_len(n.d[i]))
            U.new[[i]][t1, t2, t3] <-
              (t1 != t2) * (t1 != t3) * (t2 != t3) *
              U.xyz[[i]][t1, t2, t3] + 
              (t1 != t2) * (t2 == t3) * U.xyy[[i]][t1, t2]
    ##
    U.plus.xy <- vector("list", d)
    for (i in seq.d)
      U.plus.xy[[i]] = matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
    ##
    for (i in seq.d)
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
    Var.Lbar <- vector("list", d)
    ##
    for (i in seq.d) {
      Var.Lbar[[i]] <- matrix(rep(0, n.d[i] * n.d[i]), nrow = n.d[i])
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          Var.Lbar[[i]][t1, t2] <-
            (U.plus.xx[[i]][t1] - 2 * U.plus.xy[[i]][t1, t2] +
             U.plus.xx[[i]][t2]) / n.d[i]^2
    }
    ##
    ## Covariance of L.bar.xy. Only a subset of covariances are
    ## calculate here, i.e. cov(L.bar_(1, t1), L.bar_(1, t2))
    ##
    CoVar.Lbar <- vector("list", d)
    ##
    for (i in seq.d) {
      CoVar.Lbar[[i]] <- matrix(rep(0, n.d[i]^2), nrow = n.d[i])
      for (t1 in seq_len(n.d[i]))
        for (t2 in seq_len(n.d[i]))
          CoVar.Lbar[[i]][t1, t2] <-
            (U.plus.xx[[i]][1] - U.plus.xy[[i]][1, t2] -
             U.plus.xy[[i]][t1, 1] +
             U.plus.xy[[i]][t1, t2]) / n.d[i]^2
    }    
    ##
    ## Create y, the vector of treatment effects from each study. Only
    ## effects vs the first treatment are needed. y is coded as OR
    ## 1vsX
    ##
    Dim.y <- sum(n.d - 1)
    y <- c(rep(0, Dim.y))
    ##
    counter <- 0
    ##
    for (i in seq.d) {
      for (j in seq_len(n.d[i] - 1))
        y[j + counter] <- L.bar.xy[[i]][1, j + 1]
      ##
      counter <- counter + n.d[i] - 1 
    }    
    ##
    ## Create V, the matrix of covariances of treatment effects from
    ## each study. Only covariances of effects vs the first treatment
    ## are needed.
    ##
    V1 <- matrix(rep(0, Dim.y * Dim.y), nrow = Dim.y)
    V2 <- matrix(rep(0, Dim.y * Dim.y), nrow = Dim.y)
    ##
    counter = 0
    ##
    for (i in seq.d) {
      for (j in seq_len(n.d[i] - 1))
        for (k in seq_len(n.d[i] - 1))
          V1[j + counter, k + counter] <- (k == j) * Var.Lbar[[i]][1, j + 1]
      ##
      counter <- counter + n.d[i] - 1
    }
    ##
    counter <- 0
    ##
    for (i in seq.d) {
      for (j in 2:(n.d[i]))
        for (k in 2:(n.d[i]))
          V2[j + counter - 1, k + counter - 1] <-
            (k != j) * CoVar.Lbar[[i]][j, k]
      ##
      counter <- counter + n.d[i] - 1
    }
    ##
    V <- V1 + V2
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
    ## Define matrix X
    ##
    X <- matrix(0, nrow = Dim.y, ncol = n.treat - 1)
    ##
    for (i in seq.d)
      for (j in seq_along(ddat[[i]]$studlab))
        for (k in seq_along(trts))
          if (ddat[[i]]$treat[j] == trts[k])
            ddat[[i]]$id.t2[j] <- k
    ##
    list1 <- matrix(0, nrow = Dim.y, ncol = 2)
    ##
    N.j <- rep(0, d)
    ##
    if (d > 1)
      for (i in 2:d)
        N.j[i] <- N.j[i - 1] + n.d[i - 1] - 1
    ##
    for (i in seq.d)
      for (j in seq_len(n.d[i] - 1)) {
        list1[N.j[i] + j, 1] <- ddat[[i]]$id.t2[[1]]
        list1[N.j[i] + j, 2] <- ddat[[i]]$id.t2[[j + 1]]
      }
    ##    
    basic.contrasts <- c(2:n.treat)
    ##
    for (i in 1:Dim.y)
      for (k in 1:(n.treat - 1)) {
        if (list1[i, 1] == basic.contrasts[k])
          X[i, k] = -1
        if (list1[i, 2] == basic.contrasts[k])
          X[i, k] = 1
      }
    ##
    ## Estimate NMA 
    ##
    W <- solve(V)
    ## Basic parameters
    TE.basic <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
    ##
    d.hat <- H %*% TE.basic
    cov.d.hat <- H %*% solve(t(X) %*% W %*% X) %*% t(H)    
    ##
    ## Inconsistency global
    ##
    Q <- as.vector(t(y - X %*% TE.basic) %*% solve(V) %*% (y - X %*% TE.basic))
    df.Q <- sum(n.d - 1) - n.treat + 1
    pval.Q <- pvalQ(Q, df.Q)
    ##
    ## Drop unnecessary variables
    ##
    for (i in seq.d) {
      ddat[[i]]$id.s <- NULL
      ddat[[i]]$id.t <- NULL
      ddat[[i]]$id.t2 <- NULL
      ddat[[i]]$n.study <- NULL
    }
  }
  else if (method == "NCH") {
    ##
    ## NCH method
    ##    
    dat$ttt <- 0
    ##
    for (k in 1:n.treat)
      for (i in seq_along(dat$studlab))
        if (dat$treat[i] == trts[k])
          dat$ttt[i] <- k
    ##
    dat$treat <- as.numeric(dat$ttt)
    dat$ttt <- NULL
    ##
    dat$treat <- as.numeric(dat$treat)
    ##
    dat <- dat[order(dat$studlab, dat$treat), ] # necessary ???
    ##
    for (j in unique(dat$studlab)) {
      k <- 1
      for (i in seq_along(dat$studlab))
        if (dat$studlab[i] == j) {
          dat$count[i] <- k
          k <- k + 1
        }
    }
    ##
    max.arms <- max(dat$count)
    ##
    d1 <- reshape(dat, idvar = "studlab",
                  timevar = "count", direction = "wide")
    ##
    d1 <- d1[order(d1$studlab), ]
    ##
    d1$N1 <- 0 # d1$n.all = 0
    d1$narms <- 0
    ##
    for (i in seq_len(max.arms)) {
      ev <- paste("event.", i, sep = "")
      tot <- paste("n.", i, sep = "")
      d1$N1 <- rowSums(cbind(d1$N1, d1[, colnames(d1) == ev]),
                       na.rm = TRUE)
      ## d1$n.all <- rowSums(cbind(d1$n.all, d1[, colnames(d1) == tot]), na.rm = TRUE)
    }
    ##
    for (i in seq_along(d1$studlab))
      for (k in 1:max.arms) {
        ev <- paste("event.", k, sep = "")
        if (!is.na(d1[, colnames(d1) == ev][i]))
          d1$narms[i] <- k
      }
    ##
    for (i in seq_along(colnames(d1)))
      for (k in seq_along(colnames(d1))) {
        if (colnames(d1)[k] == paste("treat.", i, sep = ""))
          colnames(d1)[k] = paste("t", i, sep = "")
        ##
        if (colnames(d1)[k] == paste("event.", i, sep = ""))
          colnames(d1)[k] = paste("r", i, sep = "")
        ##
        if (colnames(d1)[k] == paste("n.", i, sep = ""))
          colnames(d1)[k] = paste("n", i, sep = "")
      }
    ##
    ## Likelihood function
    ##
    myLik1 <- function(mypar) {
      x <- mypar
      myLogLik1 <- rep(0, length(d1$studlab))
      myLogLik2 <- rep(0, length(d1$studlab))
      myLogLik3 <- rep(0, length(d1$studlab))  
      ##
      for (i in seq_along(d1$studlab))
        for (k in 2:d1$narms[i])
          myLogLik1[i] <-
            myLogLik1[i] + 
            d1[, colnames(d1) == paste("r", k, sep = "")][i] *
            (x[d1[, colnames(d1) == paste("t", k, sep = "")][i]] -
             x[d1$t1[i]] * (d1$t1[i] != 1))
      ##
      for (i in seq_along(d1$studlab)) {
        for (k in 2:d1$narms[i])
          myLogLik2[i] <-
            myLogLik2[i] +
            d1[, colnames(d1) == paste("n", k, sep = "")][i] *
            exp(x[d1[, colnames(d1) == paste("t", k, sep = "")][i]] -
                x[d1$t1[i]] * (d1$t1[i] != 1))
        ##
        myLogLik3[i] <- -d1$N1[i] * log(d1$n1[i] + myLogLik2[i])
      }
      ##
      myLogLik <- sum(myLogLik1 + myLogLik3)
      ##
      ## res1 <- list(myLogLik, myLogLik1, myLogLik2, myLogLik3)
      ##
      myLogLik
    }
    ##    
    init.val <- rep(0, n.treat)
    ##
    results <- optim(init.val, myLik1, method = "L-BFGS-B",
                     lower = -Inf, upper = Inf, 
                     control = list(fnscale = -1, maxit = 10000),
                     hessian = TRUE)
    ##    
    W <- solve(-results$hessian[2:n.treat, 2:n.treat])
    TE.basic <- -(results$par)[2:n.treat]
    ##
    ## Se <- sqrt(diag(W))  
    ## round(exp(TE.basic), digits = n.treat)
    ## round(exp(TE.basic + 1.96 * Se), digits = n.treat)
    ## round(exp(TE.basic - 1.96 * Se), digits = n.treat)
    ##
    ## define H matrix
    ##
    H <- matrix(rep(0, n.treat * ((n.treat - 1) / 2) * (n.treat - 1)),
                nrow = n.treat * ((n.treat - 1) / 2))
    ##
    for (i in 1:n.treat - 1)
      H[i, i] <- 1
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
    ## All relative effects
    ##
    d.hat <- H %*% TE.basic
    cov.d.hat <- H %*% W %*% t(H)
    ##
    Q <- NA
    df.Q <- NA
    pval.Q <- NA
    y <- NA
    V <- NA
    X <- NA
  }
  ##  
  TE.fixed[lower.tri(TE.fixed, diag = FALSE)] <- d.hat
  TE.fixed <- t(TE.fixed)
  TE.fixed[lower.tri(TE.fixed, diag = FALSE)] <- -d.hat
  diag(TE.fixed) <- 0
  ##
  seTE.fixed[lower.tri(seTE.fixed, diag = FALSE)] <- sqrt(diag(cov.d.hat))
  seTE.fixed <- t(seTE.fixed)
  seTE.fixed[lower.tri(seTE.fixed, diag = FALSE)] <- sqrt(diag(cov.d.hat))
  diag(seTE.fixed) <- 0
  ##
  ci.f <- ci(TE.fixed, seTE.fixed, level = level.comb)
  
  
  ##
  ##
  ## (9) Inconsistency evaluation: direct MH estimates
  ##
  ##
  C.matrix <- createB(treat1.pos, treat2.pos)
  C.matrix <- unique(C.matrix)
  colnames(C.matrix) <- trts
  ##
  for (i in seq_len(nrow(C.matrix))) {
    sel.treat1 <- colnames(C.matrix)[C.matrix[i, ] ==  1]
    sel.treat2 <- colnames(C.matrix)[C.matrix[i, ] == -1]
    selstud <- treat1 == sel.treat1 & treat2 == sel.treat2
    ##
    m.i <- metabin(event1, n1, event2, n2,
                   subset = selstud,
                   sm = "OR", MH.exact = TRUE)
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
              MH.exact = MH.exact,
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
              data.wide = dat1,
              data.long = dat,
              data.design = ddat,
              ##
              y = y, V = V, X = X,
              ##
              warn = warn,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- c("netmetabin", "netmeta")
  

  res$data <- res$data[order(res$data$...order), ]
  res$data$...order <- NULL
  rownames(res$data) <- seq_len(nrow(res$data))
  ##
  res$data.wide <- res$data.wide[order(res$data.wide$...order), ]
  res$data.wide$...order <- NULL
  rownames(res$data.wide) <- seq_len(nrow(res$data.wide))
  ##
  res$data.long <- res$data.long[order(res$data.long$...order), ]
  res$data.long$...order <- NULL
  rownames(res$data.long) <- seq_len(nrow(res$data.long))


  tab <- table(res$studlab)
  ##
  res$studies <- names(tab)
  res$narms <- as.vector(tab)
  
  
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
