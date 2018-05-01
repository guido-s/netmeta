netmeta <- function(TE, seTE,
                    treat1, treat2, studlab,
                    data = NULL, subset = NULL,
                    sm,
                    level = 0.95, level.comb = 0.95,
                    comb.fixed = gs("comb.fixed"),
                    comb.random = gs("comb.random") | !is.null(tau.preset),
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
                    n1 = NULL,
                    n2 = NULL,
                    event1 = NULL,
                    event2 = NULL,
                    ##
                    backtransf = gs("backtransf"),
                    ##
                    title = "",
                    keepdata = gs("keepdata"),
                    warn = TRUE
                    ) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkchar <- meta:::chkchar
  chklevel <- meta:::chklevel
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  ##
  chklevel(level)
  chklevel(level.comb)
  chklevel(level.predict)
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
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
  ## Catch TE, treat1, treat2, seTE, studlab from data:
  ##
  TE <- eval(mf[[match("TE", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  ##
  if (inherits(TE, "pairwise")) {
    is.pairwise <- TRUE
    ##
    sm <- attr(TE, "sm")
    ##
    seTE <- TE$seTE
    treat1 <- TE$treat1
    treat2 <- TE$treat2
    studlab <- TE$studlab
    ##
    if (!is.null(TE$n1))
      n1 <- TE$n1
    if (!is.null(TE$n2))
      n2 <- TE$n2
    if (!is.null(TE$event1))
      event1 <- TE$event1
    if (!is.null(TE$event2))
      event2 <- TE$event2
    ##
    pairdata <- TE
    data <- TE
    ##
    TE <- TE$TE
  }
  else {
    is.pairwise <- FALSE
    if (missing(sm))
      if (!is.null(data) && !is.null(attr(data, "sm")))
        sm <- attr(data, "sm")
      else
        sm <- ""
    ##
    seTE <- eval(mf[[match("seTE", names(mf))]],
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
    ##
    n1 <- eval(mf[[match("n1", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
    ##
    n2 <- eval(mf[[match("n2", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
    ##
    event1 <- eval(mf[[match("event1", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
    ##
    event2 <- eval(mf[[match("event2", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
  }
  ##
  k.Comp <- length(TE)
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
  if (!is.null(event1) & !is.null(event2))
    available.events <- TRUE
  else
    available.events <- FALSE
  ##
  if (!is.null(n1) & !is.null(n2))
    available.n <- TRUE
  else
    available.n <- FALSE
  
  
  ##
  ##
  ## (2b) Store complete dataset in list object data
  ##      (if argument keepdata is TRUE)
  ##
  ##
  if (keepdata) {
    if (nulldata & !is.pairwise)
      data <- data.frame(.studlab = studlab, stringsAsFactors = FALSE)
    else if (nulldata & is.pairwise) {
      data <- pairdata
      data$.studlab <- studlab
    }
    else
      data$.studlab <- studlab
    ##
    data$.treat1 <- treat1
    data$.treat2 <- treat2
    ##
    data$.TE <- TE
    data$.seTE <- seTE
    ##
    data$.event1 <- event1
    data$.n1 <- n1
    data$.event2 <- event2
    data$.n2 <- n2
    ##
    ## Check for correct treatment order within comparison
    ##
    wo <- data$.treat1 > data$.treat2
    ##
    if (any(wo)) {
      data$.TE[wo] <- -data$.TE[wo]
      ttreat1 <- data$.treat1
      data$.treat1[wo] <- data$.treat2[wo]
      data$.treat2[wo] <- ttreat1[wo]
      ##
      if (!is.null(data$.n1) & !is.null(data$.n2)) {
        tn1 <- data$.n1
        data$.n1[wo] <- data$.n2[wo]
        data$.n2[wo] <- tn1[wo]
      }
      ##
      if (!is.null(data$.event1) & !is.null(data$.event2)) {
        tevent1 <- data$.event1
        data$.event1[wo] <- data$.event2[wo]
        data$.event2[wo] <- tevent1[wo]
      }
    }
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
    TE <- TE[subset]
    seTE <- seTE[subset]
    treat1 <- treat1[subset]
    treat2 <- treat2[subset]
    studlab <- studlab[subset]
    ##
    if (!is.null(n1))
      n1 <- n1[subset]
    if (!is.null(n2))
      n2 <- n2[subset]
    if (!is.null(event1))
      event1 <- event1[subset]
    if (!is.null(event2))
      event2 <- event2[subset]
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
  if (!is.null(seq))
    seq <- setseq(seq, labels)
  else {
    seq <- labels
    if (is.numeric(seq))
      seq <- as.character(seq)
  }
  ##
  if (reference.group != "")
    reference.group <- setref(reference.group, labels)
  
  
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
      warning("No information given for argument 'studlab'. Assuming that comparisons are from independent studies.")
    studlab <- seq(along = TE)
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
  ## Check NAs and zero standard errors
  ##
  excl <- is.na(TE) | is.na(seTE) | seTE <= 0
  ##
  if (any(excl)) {
    data$.excl <- excl
    ##
    dat.NAs <- data.frame(studlab = studlab[excl],
                          treat1 = treat1[excl],
                          treat2 = treat2[excl],
                          TE = format(round(TE[excl], 4)),
                          seTE = format(round(seTE[excl], 4)),
                          stringsAsFactors = FALSE
                          )
    warning("Comparison",
            if (sum(excl) > 1) "s",
            " with missing TE / seTE or zero seTE not considered in network meta-analysis.",
            call. = FALSE)
    cat(paste("Comparison",
              if (sum(excl) > 1) "s",
              " not considered in network meta-analysis:\n", sep = ""))
    prmatrix(dat.NAs, quote = FALSE, right = TRUE,
             rowlab = rep("", sum(excl)))
    cat("\n")
    ##
    studlab <- studlab[!(excl)]
    treat1  <- treat1[!(excl)]
    treat2  <- treat2[!(excl)]
    TE      <- TE[!(excl)]
    seTE    <- seTE[!(excl)]
    ##
    if (!is.null(n1))
      n1 <- n1[!excl]
    if (!is.null(n2))
      n2 <- n2[!excl]
    if (!is.null(event1))
      event1 <- event1[!excl]
    if (!is.null(event2))
      event2 <- event2[!excl]
    ##
    seq <- seq[seq %in% unique(c(treat1, treat2))]
    labels <- labels[labels %in% unique(c(treat1, treat2))]
  }
  ##
  ## Check for correct number of comparisons (after removing
  ## comparisons with missing data)
  ##
  tabnarms <- table(studlab)
  sel.narms <- !is.wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  ##
  if (sum(sel.narms) == 1)
    stop(paste("After removing comparisons with missing treatment effects",
               " or standard errors,\n  study '",
               names(tabnarms)[sel.narms],
               "' has a wrong number of comparisons.",
               " Please check data and\n  consider to remove study",
               " from network meta-analysis.",
               sep = ""),
         call. = FALSE)
  if (sum(sel.narms) > 1)
    stop(paste("After removing comparisons with missing treatment effects",
               " or standard errors,\n  the following studies have",
               " a wrong number of comparisons: ",
               paste(paste("'", names(tabnarms)[sel.narms], "'", sep = ""),
                     collapse = ", "),
               "\n  Please check data and consider to remove studies",
               " from network meta-analysis.",
               sep = ""),
         call. = FALSE)
  ##
  ## Check number of subgraphs
  ##
  n.subnets <- netconnection(treat1, treat2, studlab)$n.subnets
  ##
  if (n.subnets > 1)
    stop(paste("After removing comparisons with missing treatment effects",
               " or standard errors,\n  network consists of ",
               n.subnets, " separate sub-networks.\n  ",
               "Please check data and consider to remove studies",
               " from network meta-analysis.",
               sep = ""),
         call. = FALSE)
  ##
  ## Check for correct treatment order within comparison
  ##
  wo <- treat1 > treat2
  ##
  if (any(wo)) {
    if (warn)
      warning("Note, treatments within a comparison have been re-sorted in increasing order.", call. = FALSE)
    TE[wo] <- -TE[wo]
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
    ##
    if (available.n) {
      tn1 <- n1
      n1[wo] <- n2[wo]
      n2[wo] <- tn1[wo]
    }
    ##
    if (available.events) {
      tevent1 <- event1
      event1[wo] <- event2[wo]
      event2[wo] <- tevent1[wo]
    }       
  }
  
  
  ##
  ##
  ## (5) Generate analysis dataset
  ##
  ##
  ##
  ## Generate ordered data set, with added numbers of arms per study
  ##
  p0 <- prepare(TE, seTE, treat1, treat2, studlab)
  ##
  ## Check consistency of treatment effects and standard errors in
  ## multi-arm studies
  ##
  chkmultiarm(p0$treat1, p0$treat2, p0$TE, p0$seTE, p0$studlab,
              tol = tol.multiarm, details = details.chkmultiarm)
  ##
  ## Study overview
  ##
  tdata <- data.frame(studies = p0$studlab, narms = p0$narms,
                      stringsAsFactors = FALSE)
  tdata <- unique(tdata[order(tdata$studies, tdata$narms), ])
  studies <- tdata$studies
  narms <- tdata$narms
  
  
  ##
  ##
  ## (6) Conduct network meta-analysis
  ##
  ##
  ## Fixed effect model
  ##
  res.f <- nma.ruecker(p0$TE, sqrt(1 / p0$weights),
                       p0$treat1, p0$treat2,
                       p0$treat1.pos, p0$treat2.pos,
                       p0$narms, p0$studlab,
                       sm,
                       level, level.comb,
                       p0$seTE, sep.trts = sep.trts)
  ##
  ## Random effects model
  ##
  if (is.null(tau.preset))
    tau <- res.f$tau
  else
    tau <- tau.preset
  ##
  p1 <- prepare(TE, seTE, treat1, treat2, studlab, tau)
  ##
  res.r <- nma.ruecker(p1$TE, sqrt(1 / p1$weights),
                       p1$treat1, p1$treat2,
                       p1$treat1.pos, p1$treat2.pos,
                       p1$narms, p1$studlab, 
                       sm,
                       level, level.comb,
                       p1$seTE, tau, sep.trts = sep.trts)
  ##
  TE.random <- res.r$TE.pooled
  seTE.random <- res.r$seTE.pooled
  df.Q <- res.f$df
  ##
  ## Prediction intervals
  ##
  if (df.Q == 0)
    prediction <- FALSE
  ##
  if (df.Q >= 2) {
    seTE.predict <- sqrt(seTE.random^2 + tau^2)
    ci.p <- ci(TE.random, seTE.predict, level.predict, df.Q - 1)
    p.lower <- ci.p$lower
    p.upper <- ci.p$upper
    diag(p.lower) <- 0
    diag(p.upper) <- 0
  }
  else {
    seTE.predict <- p.lower <- p.upper <- seTE.random
    seTE.predict[!is.na(seTE.predict)] <- NA
    p.lower[!is.na(p.lower)] <- NA
    p.upper[!is.na(p.upper)] <- NA
  }
  
  
  ##
  ##
  ## (7) Generate R object
  ##
  ##
  trts <- rownames(res.f$A.matrix)
  ##
  o <- order(p0$order)
  ##
  res <- list(studlab = res.f$studlab[o],
              treat1 = res.f$treat1[o],
              treat2 = res.f$treat2[o],
              ##
              TE = res.f$TE[o],
              seTE = res.f$seTE.orig[o],
              seTE.adj = res.f$seTE[o],
              ##
              n1 = n1,
              n2 = n2,
              event1 = event1,
              event2 = event2,
              ##
              k = res.f$k,
              m = res.f$m,
              n = res.f$n,
              d = NA,
              ##
              trts = trts,
              k.trts = rowSums(res.f$A.matrix),
              n.trts = if (available.n) NA else NULL,
              events.trts = if (available.events) NA else NULL,
              ##
              studies = studies,
              narms = narms,
              ##
              designs = NA,
              ##
              TE.nma.fixed = res.f$TE.nma[o],
              seTE.nma.fixed = res.f$seTE.nma[o],
              lower.nma.fixed = res.f$lower.nma[o],
              upper.nma.fixed = res.f$upper.nma[o],
              ##
              leverage.fixed = res.f$leverage[o],
              w.fixed = res.f$w.pooled[o],
              Q.fixed = res.f$Q.pooled[o],
              ##
              TE.fixed = res.f$TE.pooled,
              seTE.fixed = res.f$seTE.pooled,
              lower.fixed = res.f$lower.pooled,
              upper.fixed = res.f$upper.pooled,
              zval.fixed = res.f$zval.pooled,
              pval.fixed = res.f$pval.pooled,
              ##
              TE.nma.random = res.r$TE.nma[o],
              seTE.nma.random = res.r$seTE.nma[o],
              lower.nma.random = res.r$lower.nma[o],
              upper.nma.random = res.r$upper.nma[o],
              ##
              w.random = res.r$w.pooled[o],
              ##
              TE.random = TE.random,
              seTE.random = seTE.random,
              lower.random = res.r$lower.pooled,
              upper.random = res.r$upper.pooled,
              zval.random = res.r$zval.pooled,
              pval.random = res.r$pval.pooled,
              ##
              seTE.predict = seTE.predict,
              lower.predict = p.lower,
              upper.predict = p.upper,
              ##
              prop.direct.fixed = NA,
              prop.direct.random = NA,
              ##
              TE.direct.fixed = res.f$TE.direct,
              seTE.direct.fixed = res.f$seTE.direct,
              lower.direct.fixed = res.f$lower.direct,
              upper.direct.fixed = res.f$upper.direct,
              zval.direct.fixed = res.f$zval.direct,
              pval.direct.fixed = res.f$pval.direct,
              ##
              TE.direct.random = res.r$TE.direct,
              seTE.direct.random = res.r$seTE.direct,
              lower.direct.random = res.r$lower.direct,
              upper.direct.random = res.r$upper.direct,
              zval.direct.random = res.r$zval.direct,
              pval.direct.random = res.r$pval.direct,
              ##
              TE.indirect.fixed = NA,
              seTE.indirect.fixed = NA,
              lower.indirect.fixed = NA,
              upper.indirect.fixed = NA,
              zval.indirect.fixed = NA,
              pval.indirect.fixed = NA,
              ##
              TE.indirect.random = NA,
              seTE.indirect.random = NA,
              lower.indirect.random = NA,
              upper.indirect.random = NA,
              zval.indirect.random = NA,
              pval.indirect.random = NA,
              ##
              Q = res.f$Q,
              df.Q = df.Q,
              pval.Q = res.f$pval.Q,
              I2 = res.f$I2,
              tau = tau,
              Q.heterogeneity = NA,
              df.Q.heterogeneity = NA,
              pval.Q.heterogeneity = NA,
              Q.inconsistency = NA,
              df.Q.inconsistency = NA,
              pval.Q.inconsistency = NA,
              ##
              Q.decomp = res.f$Q.decomp,
              ##
              A.matrix = res.f$A.matrix,
              B.matrix = res.f$B.matrix[o, ],
              L.matrix = res.f$L.matrix,
              Lplus.matrix = res.f$Lplus.matrix,
              Q.matrix = res.f$Q.matrix,
              ##
              G.matrix = res.f$G.matrix[o, o],
              H.matrix = res.f$H.matrix[o, o],
              ##
              n.matrix = if (available.n) NA else NULL,
              events.matrix = if (available.events) NA else NULL,
              ##
              P.fixed = NA,
              P.random = NA,
              ##
              Cov.fixed = res.f$Cov,
              Cov.random = res.r$Cov,
              ##
              treat1.pos = res.f$treat1.pos[o],
              treat2.pos = res.f$treat2.pos[o],
              ##
              sm = sm,
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
              warn = warn,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- "netmeta"
  ##
  ## Add results for indirect treatment estimates
  ##
  n <- res$n
  ##
  res$prop.direct.fixed  <- netmeasures(res, random = FALSE)$proportion
  ## Print warning(s) in call of netmeasures() once
  oldopts <- options(warn = -1)
  res$prop.direct.random <- netmeasures(res, random = TRUE,
                                        tau.preset = res$tau,
                                        warn = FALSE)$proportion
  options(oldopts)
  if (is.logical(res$prop.direct.fixed))
    res$prop.direct.fixed <- as.numeric(res$prop.direct.fixed)
  if (is.logical(res$prop.direct.random))
    res$prop.direct.random <- as.numeric(res$prop.direct.random)
  ##
  P.fixed <- P.random <- matrix(NA, n, n)
  colnames(P.fixed) <- rownames(P.fixed) <-
    colnames(P.random) <- rownames(P.random) <- trts
  ##
  if (n == 2) {
    ##
    ## For two treatments only direct evidence is available
    ##
    res$prop.direct.fixed <- 1
    res$prop.direct.random <- 1
    names(res$prop.direct.fixed) <-
      names(res$prop.direct.random) <- paste(labels, collapse = sep.trts)
    ##
    sel <- row(P.fixed) != col(P.fixed)
    P.fixed[sel] <- 1
    P.random[sel] <- 1
  }
  else {
    k <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        k <- k + 1
        P.fixed[i, j] <- P.fixed[j, i] <- res$prop.direct.fixed[k]
        P.random[i, j] <- P.random[j, i] <- res$prop.direct.random[k]
      }
    }
  }
  ##
  ## Set direct evidence estimates to 0 if only indirect evidence is available
  ## (otherwise indirect estimates would be NA as direct estimates are NA)
  ##
  TE.direct.fixed <- res$TE.direct.fixed
  TE.direct.random <- res$TE.direct.random
  ##
  TE.direct.fixed[abs(P.fixed) < .Machine$double.eps^0.5] <- 0
  TE.direct.random[abs(P.random) < .Machine$double.eps^0.5] <- 0
  ##
  ## Indirect estimate is NA if only direct evidence is available
  ##
  res$P.fixed <- P.fixed
  res$P.random <- P.random
  ##
  P.fixed[abs(P.fixed - 1) < .Machine$double.eps^0.5] <- NA
  P.random[abs(P.random - 1) < .Machine$double.eps^0.5] <- NA
  ##
  ## Fixed effect model
  ##
  ci.if <- ci((res$TE.fixed - P.fixed * TE.direct.fixed) / (1 - P.fixed),
              sqrt(res$seTE.fixed^2 / (1 - P.fixed)),
              level = level)
  ##
  res$TE.indirect.fixed   <- ci.if$TE
  res$seTE.indirect.fixed <- ci.if$seTE
  ##
  res$lower.indirect.fixed <- ci.if$lower
  res$upper.indirect.fixed <- ci.if$upper
  ##
  res$zval.indirect.fixed <- ci.if$z
  res$pval.indirect.fixed <- ci.if$p
  ##
  ## Random effects model
  ##
  ci.ir <- ci((res$TE.random - P.random * TE.direct.random) / (1 - P.random),
              sqrt(res$seTE.random^2 / (1 - P.random)),
              level = level)
  ##
  res$TE.indirect.random   <- ci.ir$TE
  res$seTE.indirect.random <- ci.ir$seTE
  ##
  res$lower.indirect.random <- ci.ir$lower
  res$upper.indirect.random <- ci.ir$upper
  ##
  res$zval.indirect.random <- ci.ir$z
  res$pval.indirect.random <- ci.ir$p
  ##
  ## Number of designs
  ##
  krahn <- nma.krahn(res)
  res$d <- krahn$d
  if (is.null(res$d))
    res$d <- 1
  ##
  res$designs <- as.character(krahn$design$design)
  
  
  ##
  ## Calculate heterogeneity and inconsistency statistics
  ##
  if (res$d > 1) {
    dd <- decomp.design(res, warn = FALSE)
    res$Q.heterogeneity <- dd$Q.decomp$Q[2]
    res$Q.inconsistency <- dd$Q.decomp$Q[3]
    ##
    res$df.Q.heterogeneity <- dd$Q.decomp$df[2]
    res$df.Q.inconsistency <- dd$Q.decomp$df[3]
    ##
    res$pval.Q.heterogeneity <- dd$Q.decomp$pval[2]
    res$pval.Q.inconsistency <- dd$Q.decomp$pval[3]
  }
  
  
  if (keepdata) {
    if (is.null(krahn))
      ddat <- data.frame(.studlab = data$.studlab,
                         .design = paste(data$.treat1, data$.treat2,
                                        sep = sep.trts),
                         stringsAsFactors = FALSE)
    else {
      ddat <- unique(krahn$studies[, c("studlab", "design")])
      names(ddat) <- paste0(".", names(ddat))
    }
    ##
    res$data <- merge(data, ddat,
                      by = ".studlab",
                      suffixes = c(".orig", ""),
                      stringsAsFactors = FALSE)
  }
  ##
  if (available.events) {
    res$events.matrix <- netmatrix(res, event1 + event2, func = "sum")
    ##
    dat.e <- bySummary(c(event1, event2), c(treat1, treat2), long = FALSE)
    rownames(dat.e) <- dat.e$indices
    res$events.trts <- dat.e[trts, "sum"]
    names(res$events.trts) <- trts
  }
  ##
  if (available.n) {
    res$n.matrix <- netmatrix(res, n1 + n2, func = "sum")
    ##
    dat.n <- bySummary(c(n1, n2), c(treat1, treat2), long = FALSE)
    rownames(dat.n) <- dat.n$indices
    res$n.trts <- dat.n[trts, "sum"]
    names(res$n.trts) <- trts
  }
  
  
  res
}
