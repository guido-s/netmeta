pairwise <- function(treat,
                     event, n, mean, sd, TE, seTE, time,
                     data = NULL, studlab,
                     incr = 0.5, allincr = FALSE, addincr = FALSE,
                     allstudies = FALSE,
                     keep = NULL, func = "mode", ties.method = "random",
                     ...) {
  
  
  null.data <- is.null(data)
  if (null.data)
    data <- sys.frame(sys.parent())
  ##
  ## Catch studlab, treat, event, n, mean, sd, time from data:
  ##
  mf <- match.call()
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  treat <- eval(mf[[match("treat", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  event <- eval(mf[[match("event", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  n <- eval(mf[[match("n", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  mean <- eval(mf[[match("mean", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  sd <- eval(mf[[match("sd", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  TE <- eval(mf[[match("TE", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  seTE <- eval(mf[[match("seTE", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  time <- eval(mf[[match("time", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  
  
  args <- list(...)
  nam.args <- names(args)
  
  
  if (is.null(treat))
    stop("Argument 'treat' mandatory.")
  ##
  if (is.list(treat))
    chklist(treat)
  ##
  if (!is.null(event))
    if (is.list(event))
      chklist(event)
    else
      meta:::chknumeric(event)
  ##
  if (!is.null(n))
    if (is.list(n))
      chklist(n)
    else
      meta:::chknumeric(n)
  ##
  if (!is.null(mean))
    if (is.list(mean))
      chklist(mean)
    else
      meta:::chknumeric(mean)
  ##
  if (!is.null(sd))
    if (is.list(sd))
      chklist(sd)
    else
      meta:::chknumeric(sd)
  ##
  if (!is.null(TE))
    if (is.list(TE))
      chklist(TE)
    else
      meta:::chknumeric(TE)
  ##
  if (!is.null(seTE))
    if (is.list(seTE))
      chklist(seTE)
    else
      meta:::chknumeric(seTE)
  ##
  if (!is.null(time))
    if (is.list(time))
      chklist(time)
    else
      meta:::chknumeric(time)
  ##
  meta:::chknumeric(incr, min = 0, single = TRUE)
  meta:::chklogical(allincr)
  meta:::chklogical(addincr)
  meta:::chklogical(allstudies)
  
  
  if (!is.null(event) & !is.null(n) &
      is.null(mean) & is.null(sd) &
      is.null(TE) & is.null(seTE) &
      is.null(time))
    type <- "binary"
  else if (is.null(event) & !is.null(n) &
           !is.null(mean) & !is.null(sd) &
           is.null(TE) & is.null(seTE) &
           is.null(time))
    type <- "continuous"
  else if (!is.null(event) & is.null(n) &
           is.null(mean) & is.null(sd) &
           is.null(TE) & is.null(seTE) &
           !is.null(time))
    type <- "count"
  else if (is.null(event) & is.null(n) &
           is.null(mean) & is.null(sd) &
           !is.null(TE) & !is.null(seTE) &
           is.null(time))
    type <- "generic"
  else
    stop("Type of outcome unclear. Please provide the necessary information:\n  - event, n (binary outcome)\n  - n, mean, sd (continuous outcome)\n  - TE, seTE (generic outcome)\n  - event, time (incidence rates).")
  
  
  
  
  
  ##
  ## Determine whether data is in wide or long arm-based format
  ##
  if (type == "binary")
    wide.armbased <- is.list(event) & is.list(n)
  else if (type == "continuous")
    wide.armbased <- is.list(n) & is.list(mean) & is.list(sd)
  else if (type == "count")
    wide.armbased <- is.list(event) & is.list(time)
  else if (type == "generic")
    wide.armbased <- is.list(TE) & is.list(seTE)
  
  
  
  
  
  ##
  ## Transform long arm-based format to list format
  ##
  if (!wide.armbased) {
    if (is.null(studlab))
      stop("Argument 'studlab' mandatory if argument 'event' is a vector.")
    ##
    studlab <- as.character(studlab)
    ##
    treat.list <- list()
    event.list <- list()
    n.list     <- list()
    mean.list  <- list()
    sd.list    <- list()
    TE.list    <- list()
    seTE.list  <- list()
    time.list  <- list()
    ##
    if (!is.null(keep))
      if (null.data) {
        warning("Argument 'keep' ignored as argument 'data' is missing.")
        keep <- NULL
      }
      else if (any(!(keep %in% names(data))))
        stop("Argument 'keep' must contain variable names from dataset (argument 'data').")
    ##
    if (!is.null(keep) & !is.null(func)) {
      for (i in seq(along = func))
        func[i] <- meta:::setchar(func[i],
                                  c("mode", "min", "max", "mean", "median", "sum"),
                                  name = "func")
      if (length(func) == 1)
        func <- rep(func, length(keep))
    }
    ##
    ## Add variable(s) to exported dataset (defined by argument 'keep')
    ##
    if (!is.null(keep)) {
      j <- 0
        for (i in keep) {
          j <- j + 1
          ##
          dat.i <- bySummary(data[, i],
                             studlab, ties.method = ties.method)
          dat.i <- dat.i[, c("index1", func[j])]
          ##
          dat.i <- dat.i[!duplicated(dat.i), , drop = FALSE]
          ##
          if (j == 1) {
            adddata <- data.frame(studlab = dat.i$index1)
            adddata[i] <- dat.i[, 2]
          }
          else
            adddata[i] <- dat.i[, 2]
        }
    }
    else
      adddata <- NULL
    ##
    treat <- as.character(treat)
    ##
    ttab <- table(studlab, treat)
    n.arms <- apply(ttab, 1, sum)
    ##
    if (type == "binary") {
      ##
      ## Generate lists
      ##
      tdat <- data.frame(studlab, treat, event, n,
                         stringsAsFactors = FALSE)
      tdat <- tdat[order(tdat$studlab, tdat$treat), ]
      ##
      studlab <- names(n.arms)
      tres <- data.frame(studlab = studlab, stringsAsFactors = FALSE)
      ##
      for (i in 1:max(n.arms)) {
        tdat.i <- tdat[!duplicated(tdat$studlab), ]
        tres.i <- merge(tres, tdat.i, by = "studlab", all.x = TRUE)
        ##
        treat.list[[i]] <- tres.i$treat
        event.list[[i]] <- tres.i$event
        n.list[[i]]     <- tres.i$n
        ##
        tdat <- tdat[duplicated(tdat$studlab), ]
      }
      ##
      treat <- treat.list
      event <- event.list
      n     <- n.list
    }
    ##
    else if (type == "continuous") {
      ##
      ## Generate lists
      ##
      tdat <- data.frame(studlab, treat, n, mean, sd,
                         stringsAsFactors = FALSE)
      tdat <- tdat[order(tdat$studlab, tdat$treat), ]
      ##
      studlab <- names(n.arms)
      tres <- data.frame(studlab = studlab, stringsAsFactors = FALSE)
      ##
      for (i in 1:max(n.arms)) {
        tdat.i <- tdat[!duplicated(tdat$studlab), ]
        tres.i <- merge(tres, tdat.i, by = "studlab", all.x = TRUE)
        ##
        treat.list[[i]] <- tres.i$treat
        n.list[[i]]     <- tres.i$n
        mean.list[[i]]  <- tres.i$mean
        sd.list[[i]]    <- tres.i$sd
        ##
        tdat <- tdat[duplicated(tdat$studlab), ]
      }
      ##
      treat <- treat.list
      n     <- n.list
      mean  <- mean.list
      sd    <- sd.list
    }
    ##
    else if (type == "count") {
      ##
      ## Generate lists
      ##
      tdat <- data.frame(studlab, treat, event, time,
                         stringsAsFactors = FALSE)
      tdat <- tdat[order(tdat$studlab, tdat$treat), ]
      ##
      studlab <- names(n.arms)
      tres <- data.frame(studlab = studlab, stringsAsFactors = FALSE)
      ##
      for (i in 1:max(n.arms)) {
        tdat.i <- tdat[!duplicated(tdat$studlab), ]
        tres.i <- merge(tres, tdat.i, by = "studlab", all.x = TRUE)
        ##
        treat.list[[i]] <- tres.i$treat
        event.list[[i]] <- tres.i$event
        time.list[[i]]  <- tres.i$time
        ##
        tdat <- tdat[duplicated(tdat$studlab), ]
      }
      ##
      treat <- treat.list
      event <- event.list
      time  <- time.list
    }
    ##
    else if (type == "generic") {
      ##
      ## Generate lists
      ##
      tdat <- data.frame(studlab, treat, TE, seTE,
                         stringsAsFactors = FALSE)
      tdat <- tdat[order(tdat$studlab, tdat$treat), ]
      ##
      studlab <- names(n.arms)
      tres <- data.frame(studlab = studlab, stringsAsFactors = FALSE)
      ##
      for (i in 1:max(n.arms)) {
        tdat.i <- tdat[!duplicated(tdat$studlab), ]
        tres.i <- merge(tres, tdat.i, by = "studlab", all.x = TRUE)
        ##
        treat.list[[i]] <- tres.i$treat
        TE.list[[i]]    <- tres.i$TE
        seTE.list[[i]]  <- tres.i$seTE
        ##
        tdat <- tdat[duplicated(tdat$studlab), ]
      }
      ##
      treat <- treat.list
      TE    <- TE.list
      seTE  <- seTE.list
    }
  }
  
  
  ##
  ## Check and set study labels
  ##
  if (is.null(studlab))
    studlab <- seq(along = treat[[1]])
  ##
  if (length(studlab) != length(unique(studlab)))
    stop("Study labels must all be distinct.")
  ##
  levs <- studlab
  
  
  narms <- length(treat)
  nstud <- length(studlab)
  
  
  ##
  ## Auxiliary function to calculate number of arms with zero events
  ## per study
  ##
  sumzero <- function(x) sum(x[!is.na(x)] == 0)
  
  
  if (type == "binary") {
    if (length(event) != narms)
      stop("Different length of lists 'treat' and 'event'.")
    if (length(n) != narms)
      stop("Different length of lists 'treat' and 'n'.")
    ##
    ## Determine increment for individual studies
    ##
    n.zeros <- apply(matrix(unlist(event), ncol = length(event)), 1, sumzero)
    n.all   <- apply(matrix(unlist(n), ncol = length(event)) -
                     matrix(unlist(event), ncol = length(event)),
                     1, sumzero)
    ##
    incr.study <- rep(0, length(n.zeros))
    ##
    if ("sm" %in% nam.args)
      sm <- args$sm
    else
      sm <- gs("smbin")
    ##
    sm <- meta:::setchar(sm, c("OR", "RD", "RR", "ASD"))
    ##
    sparse <- switch(sm,
                     OR = (n.zeros > 0) | (n.all > 0),
                     RD = (n.zeros > 0) | (n.all > 0),
                     RR = (n.zeros > 0) | (n.all > 0),
                     ASD = rep(FALSE, length(n.zeros)))
    ##
    if (!allincr & !addincr)
      incr.study[sparse] <- incr
    else if (addincr)
      incr.study[] <- incr
    else {
      if (any(n.zeros > 0))
        incr.study[] <- incr
      else
        incr.study[] <- 0
    }
    ##
    for (i in 1:(narms - 1)) {
      ##
      if (i == 1 & (length(treat[[i]]) != length(event[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'event'.")
      if (i == 1 & (length(event[[i]]) != length(n[[i]])))
        stop("Different length of element ", i, " of lists 'event' and 'n'.")
      ##
      for (j in (i + 1):narms) {
        ##
        if (length(treat[[j]]) != length(event[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'event'.")
        if (length(event[[j]]) != length(n[[j]]))
          stop("Different length of element ", j, " of lists 'event' and 'n'.")
        ##
        dat <- data.frame(TE = NA, seTE = NA,
                          studlab = studlab,
                          treat1 = treat[[i]],
                          treat2 = treat[[j]],
                          event1 = event[[i]], n1 = n[[i]],
                          event2 = event[[j]], n2 = n[[j]],
                          incr = incr.study,
                          allstudies = allstudies,
                          stringsAsFactors = FALSE)
        ##
        if (wide.armbased) {
          dat <- cbind(dat, data, stringsAsFactors = FALSE)
          dupl <- duplicated(names(dat))
          if (any(dupl))
            names(dat)[dupl] <- paste(names(dat)[dupl], "orig", sep = ".")
        }
        else if (!is.null(adddata))
          dat <- merge(dat, adddata, by = "studlab",
                       suffixes = c("",".orig"))
        ##
        dat <- dat[!(is.na(dat$event1) & is.na(dat$n1)), ]
        dat <- dat[!(is.na(dat$event2) & is.na(dat$n2)), ]
        ##
        if (nrow(dat) > 0) {
          m1 <- metabin(dat$event1, dat$n1,
                        dat$event2, dat$n2,
                        incr = dat$incr, addincr = TRUE,
                        allstudies = allstudies,
                        ...)
          dat$TE   <- m1$TE
          dat$seTE <- m1$seTE
          ##
          dat.NAs <- dat[is.na(dat$TE) | is.na(dat$seTE) | dat$seTE <= 0, ]
          ##
          if (i == 1 & j == 2) {
            res <- dat
            res.NAs <- dat.NAs
          }
          else {
            res <- rbind(res, dat)
            res.NAs <- rbind(res.NAs, dat.NAs)
          }
        }
        else
          if (i == 1 & j == 2)
            stop("No studies available for comparison of first and second treatment.")
      }
    }
  }
  
  
  if (type == "continuous") {
    if (length(n) != narms)
      stop("Different length of lists 'treat' and 'n'.")
    if (length(mean) != narms)
      stop("Different length of lists 'treat' and 'mean'.")
    if (length(sd) != narms)
      stop("Different length of lists 'treat' and 'sd'.")
    ##
    for (i in seq_len(narms)) {
      ##
      if (length(treat[[i]]) != length(n[[i]]))
        stop("Different length of element ", i, " of lists 'treat' and 'n'.",
             call. = FALSE)
      if (length(treat[[i]]) != length(mean[[i]]))
        stop("Different length of element ", i, " of lists 'treat' and 'mean'.",
             call. = FALSE)
      if (length(treat[[i]]) != length(sd[[i]]))
        stop("Different length of element ", i, " of lists 'treat' and 'sd'.",
             call. = FALSE)
      if (length(treat[[i]]) != nstud)
        stop("Different length of study labels and element ", i, " of list 'treat'.",
             call. = FALSE)
    }
    ##
    ## For standardized mean difference, calculate pooled standard
    ## deviation for multi-arm studies
    ##
    if ("sm" %in% nam.args && (tolower(args$sm) == "smd" & narms > 2)) {
      pooled.sd <- function(sd, n) {
        sel <- !is.na(sd) & !is.na(n)
        ##
        if (any(sel))
          res <- sqrt(sum((n[sel] - 1) * sd[sel]^2) / sum(n[sel] - 1))
        else
          res <- NA
        ##
        res
      }
      ##
      N <- matrix(unlist(n), ncol = narms, nrow = nstud, byrow = FALSE)
      M <- matrix(unlist(mean), ncol = narms, nrow = nstud, byrow = FALSE)
      S <- matrix(unlist(sd), ncol = narms, nrow = nstud, byrow = FALSE)
      ##
      sel.n <- apply(!is.na(N) & N > 0, 1, sum) > 2
      sel.mean <- apply(!is.na(M), 1, sum) > 2
      sel.sd <- apply(!is.na(S) & S > 0, 1, sum) > 2
      sel <- sel.n & sel.mean & sel.sd
      ##
      if (any(sel)) {
        N <- N[sel, , drop = FALSE]
        S <- S[sel, , drop = FALSE]
        sd.p <- rep_len(NA, nrow(N))
        ##
        for (i in seq_len(nrow(N)))
          sd.p[i] <- pooled.sd(S[i, ], N[i, ])
      }
      ##
      for (i in seq_len(narms))
        sd[[i]][sel] <- ifelse(is.na(sd[[i]][sel]), NA, sd.p)
    }
    ##
    for (i in 1:(narms - 1)) {
      for (j in (i + 1):narms) {
        dat <- data.frame(TE = NA, seTE = NA,
                          studlab = studlab,
                          treat1 = treat[[i]],
                          treat2 = treat[[j]],
                          n1 = n[[i]], mean1 = mean[[i]], sd1 = sd[[i]],
                          n2 = n[[j]], mean2 = mean[[j]], sd2 = sd[[j]],
                          stringsAsFactors = FALSE)
        ##
        if (wide.armbased) {
          dat <- cbind(dat, data, stringsAsFactors = FALSE)
          dupl <- duplicated(names(dat))
          if (any(dupl))
            names(dat)[dupl] <- paste(names(dat)[dupl], "orig", sep = ".")
        }
        else if (!is.null(adddata))
          dat <- merge(dat, adddata, by = "studlab",
                       suffixes = c("",".orig"))
        ##
        dat <- dat[!(is.na(dat$n1) & is.na(dat$mean1) & is.na(dat$sd1)), ]
        dat <- dat[!(is.na(dat$n2) & is.na(dat$mean2) & is.na(dat$sd2)), ]
        ##
        if (nrow(dat) > 0) {
          m1 <- metacont(dat$n1, dat$mean1, dat$sd1,
                         dat$n2, dat$mean2, dat$sd2,
                         ...)
          dat$TE   <- m1$TE
          dat$seTE <- m1$seTE
          ##
          dat.NAs <- dat[is.na(dat$TE) | is.na(dat$seTE) | dat$seTE <= 0, ]
          ##
          if (i == 1 & j == 2) {
            res <- dat
            res.NAs <- dat.NAs
          }
          else {
            res <- rbind(res, dat)
            res.NAs <- rbind(res.NAs, dat.NAs)
          }
        }
        else
          if (i == 1 & j == 2)
            stop("No studies available for comparison of first and second treatment.")
      }
    }
  }
  
  
  if (type == "generic") {
    if (length(TE) != narms)
      stop("Different length of lists 'treat' and 'TE'.")
    if (length(seTE) != narms)
      stop("Different length of lists 'treat' and 'seTE'.")
    ##
    for (i in 1:(narms - 1)) {
      ##
      if (i == 1 & (length(treat[[i]]) != length(TE[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'TE'.")
      if (i == 1 & (length(treat[[i]]) != length(seTE[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'seTE'.")
      ##
      for (j in (i + 1):narms) {
        ##
        if (length(treat[[j]]) != length(TE[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'TE'.")
        if (length(treat[[j]]) != length(seTE[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'seTE'.")
        ##
        dat <- data.frame(TE = NA, seTE = NA,
                          studlab = studlab,
                          treat1 = treat[[i]],
                          treat2 = treat[[j]],
                          TE1 = TE[[i]], seTE1 = seTE[[i]],
                          TE2 = TE[[j]], seTE2 = seTE[[j]],
                          stringsAsFactors = FALSE)
        ##
        if (wide.armbased) {
          dat <- cbind(dat, data, stringsAsFactors = FALSE)
          dupl <- duplicated(names(dat))
          if (any(dupl))
            names(dat)[dupl] <- paste(names(dat)[dupl], "orig", sep = ".")
        }
        else if (!is.null(adddata))
          dat <- merge(dat, adddata, by = "studlab",
                       suffixes = c("",".orig"))
        ##
        dat <- dat[!(is.na(dat$TE1) & is.na(dat$seTE1)), ]
        dat <- dat[!(is.na(dat$TE2) & is.na(dat$seTE2)), ]
        ##
        if (nrow(dat) > 0) {
          m1 <- metagen(dat$TE1 - dat$TE2,
                        sqrt(dat$seTE1^2 + dat$seTE2^2), ...)
          dat$TE <- m1$TE
          dat$seTE <- m1$seTE
          ##
          dat.NAs <- dat[is.na(dat$TE) | is.na(dat$seTE) | dat$seTE <= 0, ]
          ##
          if (i == 1 & j == 2) {
            res <- dat
            res.NAs <- dat.NAs
          }
          else {
            res <- rbind(res, dat)
            res.NAs <- rbind(res.NAs, dat.NAs)
          }
        }
        else
          if (i == 1 & j == 2)
            stop("No studies available for comparison of first and second treatment.")
      }
    }
  }
  
  
  if (type == "count") {
    if (length(event) != narms)
      stop("Different length of lists 'treat' and 'event'.")
    if (length(time) != narms)
      stop("Different length of lists 'treat' and 'time'.")
    ##
    ## Determine increment for individual studies
    ##
    n.zeros <- apply(matrix(unlist(event), ncol = length(event)), 1, sumzero)
    ##
    incr.study <- rep(0, length(n.zeros))
    ##
    sparse <- n.zeros > 0
    ##
    if (!allincr & !addincr)
      incr.study[sparse] <- incr
    else if (addincr)
      incr.study[] <- incr
    else {
      if (any(n.zeros > 0))
        incr.study[] <- incr
      else
        incr.study[] <- 0
    }
    ##
    for (i in 1:(narms - 1)) {
      ##
      if (i == 1 & (length(treat[[i]]) != length(event[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'event'.")
      if (i == 1 & (length(treat[[i]]) != length(time[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'time'.")
      ##
      for (j in (i + 1):narms) {
        ##
        if (length(treat[[j]]) != length(event[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'event'.")
        if (length(treat[[j]]) != length(time[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'time'.")
        ##
        dat <- data.frame(TE = NA, seTE = NA,
                          studlab = studlab,
                          treat1 = treat[[i]],
                          treat2 = treat[[j]],
                          event1 = event[[i]], time1 = time[[i]],
                          event2 = event[[j]], time2 = time[[j]],
                          incr = incr.study,
                          stringsAsFactors = FALSE)
        ##
        if (wide.armbased) {
          dat <- cbind(dat, data, stringsAsFactors = FALSE)
          dupl <- duplicated(names(dat))
          if (any(dupl))
            names(dat)[dupl] <- paste(names(dat)[dupl], "orig", sep = ".")
        }
        else if (!is.null(adddata))
          dat <- merge(dat, adddata, by = "studlab",
                       suffixes = c("",".orig"))
        ##
        dat <- dat[!(is.na(dat$event1) & is.na(dat$time1)), ]
        dat <- dat[!(is.na(dat$event2) & is.na(dat$time2)), ]
        ##
        if (nrow(dat) > 0) {
          m1 <- metainc(dat$event1, dat$time1,
                        dat$event2, dat$time2,
                        incr = dat$incr, addincr = TRUE,
                        allstudies = allstudies,
                        ...)
          dat$TE <- m1$TE
          dat$seTE <- m1$seTE
          ##
          dat.NAs <- dat[is.na(dat$TE) | is.na(dat$seTE) | dat$seTE <= 0, ]
          ##
          if (i == 1 & j == 2) {
            res <- dat
            res.NAs <- dat.NAs
          }
          else {
            res <- rbind(res, dat)
            res.NAs <- rbind(res.NAs, dat.NAs)
          }
        }
        else
          if (i == 1 & j == 2)
            stop("No studies available for comparison of first and second treatment.")
      }
    }
  }
  
  
  ##
  ## Additional checks
  ##
  ##
  ## a) Duplicate treatments ?
  ##
  sel.treat <- as.character(res$treat1) == as.character(res$treat2)
  ##
  if (any(sel.treat))
    stop(paste("Identical treatments for the following studies:\n  ",
               paste(paste("'", unique(sort(res$studlab[sel.treat])),
                           "'", sep = ""),
                     collapse = " - "), sep = ""))
  ##
  ## b) Studies missing ?
  ##
  sel.study <- !(studlab %in% unique(as.character(res$studlab)))
  ##
  if (any(sel.study))
    warning(paste("The following studies are not considered in the analysis\n  ",
                  "(due to single study arm or missing values):\n  ",
                  paste(paste("'", studlab[sel.study], "'", sep = ""),
                        collapse = " - "), sep = ""))
  ##
  ## c) Missing treatment estimates or standard errors?
  ##
  if (nrow(res.NAs) > 0) {
    warning("Comparison",
            if (nrow(res.NAs) > 1) "s",
            " with missing TE / seTE or zero seTE",
            " will not be considered in network meta-analysis.",
            call. = FALSE)
    cat(paste("Comparison",
              if (nrow(res.NAs) > 1) "s",
              " will not be considered in network meta-analysis:\n",
              sep = ""))
    ##
    prmatrix(res.NAs,
             quote = FALSE, right = TRUE, na.print = "NA",
             rowlab = rep("", nrow(res.NAs)))
  }
  
  
  attr(res, "sm") <- m1$sm
  attr(res, "method") <- m1$method
  attr(res, "version") <- packageDescription("netmeta")$Version
  
  
  res <- res[order(factor(res$studlab, levels = levs), res$treat1, res$treat2), ]
  ##
  rownames(res) <- 1:nrow(res)
  class(res) <- c(class(res), "pairwise")
  res
}
