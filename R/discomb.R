discomb <- function(TE, seTE,
                    treat1, treat2,
                    studlab, data = NULL, subset = NULL,
                    ##
                    inactive = NULL,
                    sep.comps = "+", 
                    C.matrix,
                    ##
                    sm,
                    level = gs("level"),
                    level.comb = gs("level.comb"),
                    comb.fixed = gs("comb.fixed"),
                    comb.random = gs("comb.random") | !is.null(tau.preset),
                    ##
                    reference.group = "",
                    baseline.reference = TRUE,
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
                    warn = TRUE) {
  
  
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
  chkchar(sep.comps, nchar = 1)
  ##
  chklevel(level)
  chklevel(level.comb)
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  ##
  chklogical(baseline.reference)
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
  chklogical(warn)
  
  
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
  k.Comp <- length(TE)
  ##
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  if (is.factor(studlab))
    studlab <- as.character(studlab)
  
  
  ##
  ##
  ## (3) Use subset for analysis
  ##
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  if (!is.null(subset)) {
    if ((is.logical(subset) & (sum(subset) > k.Comp)) ||
        (length(subset) > k.Comp))
      stop("Length of subset is larger than number of studies.")
    ##
    TE <- TE[subset]
    seTE <- seTE[subset]
    treat1 <- treat1[subset]
    treat2 <- treat2[subset]
    studlab <- studlab[subset]
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
    stop("Treatments must be different (arguments 'treat1' and 'treat2').")
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
               sep = ""))
  if (sum(sel.narms) > 1)
    stop(paste("The following studies have a wrong number of comparisons: ",
               paste(paste("'", names(tabnarms)[sel.narms], "'", sep = ""),
                     collapse = ", "),
               "\n  Please provide data for all treatment comparisons",
               " (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
               sep = ""))
  ##
  ## Check NAs and zero standard errors
  ##
  excl <- is.na(TE) | is.na(seTE) | seTE <= 0
  ##
  if (any(excl)) {
    dat.NAs <- data.frame(studlab = studlab[excl],
                          treat1 = treat1[excl],
                          treat2 = treat2[excl],
                          TE = format(round(TE[excl], 4)),
                          seTE = format(round(seTE[excl], 4))
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
               sep = ""))
  if (sum(sel.narms) > 1)
    stop(paste("After removing comparisons with missing treatment effects",
               " or standard errors,\n  the following studies have",
               " a wrong number of comparisons: ",
               paste(paste("'", names(tabnarms)[sel.narms], "'", sep = ""),
                     collapse = ", "),
               "\n  Please check data and consider to remove studies",
               " from network meta-analysis.",
               sep = ""))
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
  }
  
  
  ##
  ##
  ## (5) Create C.matrix
  ##
  ##
  netc <- netconnection(treat1, treat2, studlab)
  ##
  if (missing(C.matrix)) {
    C.matrix <- createC(netc, sep.comps, inactive)
    C.matrix <- C.matrix[labels, , drop = FALSE]
  }
  else {
    ##
    if (!(is.matrix(C.matrix) | is.data.frame(C.matrix))) 
      stop("Argument 'C.matrix' must be a matrix or data frame.", 
           call. = FALSE)
    ##
    wrong.labels <- FALSE
    if (is.null(rownames(C.matrix)))
      wrong.labels <- TRUE
    else {
      if (length(unique(labels)) == length(unique(tolower(labels)))) 
        idx <- charmatch(tolower(rownames(C.matrix)), 
                         tolower(labels), nomatch = NA)
      else idx <- charmatch(rownames(C.matrix), labels, nomatch = NA)
      if (any(is.na(idx)) || any(idx == 0)) 
        wrong.labels <- TRUE
    }
    if (wrong.labels) 
      stop(paste("Row names of argument 'C.matrix' must be a ", 
                 "permutation of treatment names:\n  ",
                 paste(paste("'", labels, "'", sep = ""), collapse = " - "),
                 sep = ""),
           call. = FALSE)
    ##
    C.matrix <- C.matrix[labels, , drop = FALSE]
  }
  if (is.data.frame(C.matrix))
    C.matrix <- as.matrix(C.matrix)
  ##
  c <- ncol(C.matrix) # number of components
  
  
  p0 <- prepare(TE, seTE, treat1, treat2, studlab)
  ##
  o <- order(p0$order)
  ##
  B.matrix <- createB(p0$treat1.pos[o], p0$treat2.pos[o])
  ##
  colnames(B.matrix) <- labels
  rownames(B.matrix) <- studlab
  
  
  ##
  ## Design matrix based on treatment components
  ##
  X <- B.matrix %*% C.matrix
  ##
  colnames(X) <- colnames(C.matrix)
  rownames(X) <- studlab
  
  
  tdata <- data.frame(studies = p0$studlab[o], narms = p0$narms[o])
  tdata <- unique(tdata[order(tdata$studies, tdata$narms), ])
  ##
  studies <- tdata$studies
  narms <- tdata$narms
  n.a <- sum(narms)  
  ##
  comps <- colnames(C.matrix) # treatment components
  ##
  n <- length(labels)
  m <- length(TE)
  k <- length(unique(studlab))
  
  
  ##
  ## Fixed effects models
  ##
  df.Q.additive <- n.a - k - qr(X)$rank
  ##
  if (netc$n.subnets == 1) {
    net <- netmeta(TE, seTE, treat1, treat2, studlab)
    ##
    Q <- net$Q
    df.Q <- net$df.Q
    pval.Q <- net$pval.Q
    ##
    df.Q.diff <- n - 1 - qr(X)$rank
  }
  else {
    Q <- df.Q <- pval.Q <- NA
    df.Q.diff <- NA
  }
  
  
  res.f <- nma.additive(p0$TE[o], p0$weights[o], p0$studlab[o],
                        p0$treat1[o], p0$treat2[o], level.comb,
                        X, C.matrix, B.matrix,
                        Q, df.Q.additive, df.Q.diff,
                        n, sep.trts)
  
  
  ##
  ## Calculate heterogeneity statistics (additive model)
  ##
  Q.additive <- res.f$Q.additive
  ##
  if (!is.null(tau.preset))
    tau <- tau.preset
  else
    tau <- res.f$tau
  ##
  I2 <- res.f$I2
  
  
  ##
  ## Random effects models
  ##
  p1 <- prepare(TE, seTE, treat1, treat2, studlab, tau)
  ##
  res.r <- nma.additive(p1$TE[o], p1$weights[o], p1$studlab[o],
                        p1$treat1[o], p1$treat2[o], level.comb,
                        X, C.matrix, B.matrix,
                        Q, df.Q.additive, df.Q.diff,
                        n, sep.trts)
  
  
  NAs <- rep(NA, length(res.f$comparisons$TE))
  
  
  res <- list(studlab = p0$studlab[o],
              treat1 = p0$treat1[o],
              treat2 = p0$treat2[o],
              ##
              TE = p0$TE[o],
              seTE = seTE[o],
              seTE.adj = sqrt(1 / p0$weights[o]),
              ##
              event1 = NA,
              event2 = NA,
              n1 = NA,
              n2 = NA,
              ##
              k = k,
              m = m,
              n = n,
              d = NA,
              c = c,
              ##
              trts = labels,
              k.trts = NA,
              n.trts = NA,
              events.trts = NA,
              ##
              studies = NA,
              narms = NA,
              ##
              designs = NA,
              ##
              comps = names(res.f$components$TE),
              k.comps = NA,
              n.comps = NA,
              events.comps = NA,
              ##
              TE.nma.fixed = NAs,
              seTE.nma.fixed = NAs,
              lower.nma.fixed = NAs,
              upper.nma.fixed = NAs,
              zval.nma.fixed = NAs,
              pval.nma.fixed = NAs,
              ##
              TE.cnma.fixed = res.f$comparisons$TE,
              seTE.cnma.fixed = res.f$comparisons$seTE,
              lower.cnma.fixed = res.f$comparisons$lower,
              upper.cnma.fixed = res.f$comparisons$upper,
              zval.cnma.fixed = res.f$comparisons$z,
              pval.cnma.fixed = res.f$comparisons$p,
              ##
              TE.fixed = res.f$all.comparisons$TE,
              seTE.fixed = res.f$all.comparisons$seTE,
              lower.fixed = res.f$all.comparisons$lower,
              upper.fixed = res.f$all.comparisons$upper,
              zval.fixed = res.f$all.comparisons$z,
              pval.fixed = res.f$all.comparisons$p,
              ##
              TE.nma.random = NAs,
              seTE.nma.random = NAs,
              lower.nma.random = NAs,
              upper.nma.random = NAs,
              zval.nma.random = NAs,
              pval.nma.random = NAs,
              ##
              TE.cnma.random = res.r$comparisons$TE,
              seTE.cnma.random = res.r$comparisons$seTE,
              lower.cnma.random = res.r$comparisons$lower,
              upper.cnma.random = res.r$comparisons$upper,
              zval.cnma.random = res.r$comparisons$z,
              pval.cnma.random = res.r$comparisons$p,
              ##
              TE.random = res.r$all.comparisons$TE,
              seTE.random = res.r$all.comparisons$seTE,
              lower.random = res.r$all.comparisons$lower,
              upper.random = res.r$all.comparisons$upper,
              zval.random = res.r$all.comparisons$z,
              pval.random = res.r$all.comparisons$p,
              ##
              Comp.fixed = unname(res.f$components$TE),
              seComp.fixed = unname(res.f$components$seTE),
              lower.Comp.fixed = unname(res.f$components$lower),
              upper.Comp.fixed = unname(res.f$components$upper),
              zval.Comp.fixed = unname(res.f$components$z),
              pval.Comp.fixed = unname(res.f$components$p),
              ##
              Comp.random = unname(res.r$components$TE),
              seComp.random = unname(res.r$components$seTE),
              lower.Comp.random = unname(res.r$components$lower),
              upper.Comp.random = unname(res.r$components$upper),
              zval.Comp.random = unname(res.r$components$z),
              pval.Comp.random = unname(res.r$components$p),
              ##
              Comb.fixed = unname(res.f$combinations$TE),
              seComb.fixed = unname(res.f$combinations$seTE),
              lower.Comb.fixed = unname(res.f$combinations$lower),
              upper.Comb.fixed = unname(res.f$combinations$upper),
              zval.Comb.fixed = unname(res.f$combinations$z),
              pval.Comb.fixed = unname(res.f$combinations$p),
              ##
              Comb.random = unname(res.r$combinations$TE),
              seComb.random = unname(res.r$combinations$seTE),
              lower.Comb.random = unname(res.r$combinations$lower),
              upper.Comb.random = unname(res.r$combinations$upper),
              zval.Comb.random = unname(res.r$combinations$z),
              pval.Comb.random = unname(res.r$combinations$p),
              ##
              Q.additive = Q.additive, 
              df.Q.additive = df.Q.additive, 
              pval.Q.additive = res.f$pval.Q.additive,
              tau = tau,
              I2 = I2,
              ##
              Q.standard = Q,
              df.Q.standard = df.Q,
              pval.Q.standard = pval.Q, 
              ##
              Q.diff = res.f$Q.diff,
              df.Q.diff = df.Q.diff,
              pval.Q.diff = res.f$pval.Q.diff, 
              ##
              B.matrix = B.matrix,
              C.matrix = C.matrix,
              ##
              n.matrix = NA,
              events.matrix = NA,
              ##
              sm = sm,
              method = "Inverse",
              level = level,
              level.comb = level.comb,
              comb.fixed = comb.fixed,
              comb.random = comb.random, 
              ##
              reference.group = reference.group,
              baseline.reference = baseline.reference,
              all.treatments = NULL,
              seq = seq,
              ##
              tau.preset = tau.preset,
              ##
              sep.trts = sep.trts,
              nchar.trts = nchar.trts,
              ##
              inactive = inactive,
              sep.comps = sep.comps,
              ##
              backtransf = backtransf, 
              ##
              title = title,
              ##
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- c("discomb", "netcomb")
  
  
  res
}
