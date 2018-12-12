pairwise <- function(treat,
                     event, n, mean, sd, TE, seTE, time,
                     data = NULL, studlab,
                     incr = 0.5, allincr = FALSE, addincr = FALSE,
                     allstudies = FALSE, warn = gs("warn"),
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
  meta:::chklogical(warn)


  if (!is.null(TE) & !is.null(seTE))
    type <- "generic"
  else if (!is.null(event) & !is.null(time) &
           is.null(mean) & is.null(sd))
    type <- "count"
  else if (!is.null(event) & !is.null(n) &
           is.null(mean) & is.null(sd))
    type <- "binary"
  else if (!is.null(n) & !is.null(mean) & !is.null(sd))
    type <- "continuous"
  else
    stop("Type of outcome unclear. Please provide the necessary ",
         "information:\n  - event, n (binary outcome)\n  - n, ",
         "mean, sd (continuous outcome)\n  - TE, seTE (generic outcome)\n",
         "  - event, time (incidence rates).")





  ##
  ## Determine whether data is in wide or long arm-based format
  ##
  if (type == "generic")
    wide.armbased <- is.list(TE) & is.list(seTE)
  if (type == "binary")
    wide.armbased <- is.list(event) & is.list(n)
  else if (type == "continuous")
    wide.armbased <- is.list(n) & is.list(mean) & is.list(sd)
  else if (type == "count")
    wide.armbased <- is.list(event) & is.list(time)





  ##
  ## Transform long arm-based format to list format
  ##
  if (!wide.armbased) {
    if (is.null(studlab))
      stop("Argument 'studlab' mandatory if argument 'event' is a vector.")
    ##
    studlab <- as.character(studlab)
    ##
    treat <- as.character(treat)
    ##
    ttab <- table(studlab, treat)
    n.arms <- apply(ttab, 1, sum)
    max.arms <- max(n.arms)
    ##
    treat.list <- vector("list", max.arms)
    event.list <- vector("list", max.arms)
    n.list     <- vector("list", max.arms)
    mean.list  <- vector("list", max.arms)
    sd.list    <- vector("list", max.arms)
    TE.list    <- vector("list", max.arms)
    seTE.list  <- vector("list", max.arms)
    time.list  <- vector("list", max.arms)
    ##
    if (!null.data)
      adddata <- vector("list", max.arms)
    ##
    if (type == "binary") {
      ##
      ## Generate lists
      ##
      tdat <- data.frame(studlab, treat, event, n,
                         .order = seq_along(studlab),
                         stringsAsFactors = FALSE)
      ##
      if (!null.data) {
        tdat <- cbind(tdat, data)
        dupl <- duplicated(names(tdat))
        if (any(dupl))
          names(tdat)[dupl] <- paste(names(tdat)[dupl], "orig", sep = ".")
      }
      ##
      studlab <- names(n.arms)
      dat.studlab <- data.frame(studlab, stringsAsFactors = FALSE)
      ##
      for (i in 1:max.arms) {
        sel.i <- !duplicated(tdat$studlab)
        tdat.i <- merge(dat.studlab, tdat[sel.i, ],
                        by = "studlab", all.x = TRUE)
        ##
        treat.list[[i]] <- tdat.i$treat
        event.list[[i]] <- tdat.i$event
        n.list[[i]]     <- tdat.i$n
        ##
        tdat.i$event <- NULL
        tdat.i$n     <- NULL
        ##
        if (!null.data)
          adddata[[i]] <- tdat.i
        ##
        tdat <- tdat[!sel.i, ]
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
                         .order = seq_along(studlab),
                         stringsAsFactors = FALSE)
      ##
      if (!null.data) {
        tdat <- cbind(tdat, data)
        dupl <- duplicated(names(tdat))
        if (any(dupl))
          names(tdat)[dupl] <- paste(names(tdat)[dupl], "orig", sep = ".")
      }
      ##
      studlab <- names(n.arms)
      dat.studlab <- data.frame(studlab, stringsAsFactors = FALSE)
      ##
      for (i in 1:max.arms) {
        sel.i <- !duplicated(tdat$studlab)
        tdat.i <- merge(dat.studlab, tdat[sel.i, ],
                        by = "studlab", all.x = TRUE)
        ##
        treat.list[[i]] <- tdat.i$treat
        n.list[[i]]     <- tdat.i$n
        mean.list[[i]]  <- tdat.i$mean
        sd.list[[i]]    <- tdat.i$sd
        ##
        tdat.i$n    <- NULL
        tdat.i$mean <- NULL
        tdat.i$sd   <- NULL
        ##
        if (!null.data)
          adddata[[i]] <- tdat.i
        ##
        tdat <- tdat[!sel.i, ]
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
                         .order = seq_along(studlab),
                         stringsAsFactors = FALSE)
      ##
      if (!is.null(n))
        tdat$n <- n
      ##
      if (!null.data) {
        tdat <- cbind(tdat, data)
        dupl <- duplicated(names(tdat))
        if (any(dupl))
          names(tdat)[dupl] <- paste(names(tdat)[dupl], "orig", sep = ".")
      }
      ##
      studlab <- names(n.arms)
      dat.studlab <- data.frame(studlab, stringsAsFactors = FALSE)
      ##
      for (i in 1:max.arms) {
        sel.i <- !duplicated(tdat$studlab)
        tdat.i <- merge(dat.studlab, tdat[sel.i, ],
                        by = "studlab", all.x = TRUE)
        ##
        treat.list[[i]] <- tdat.i$treat
        event.list[[i]] <- tdat.i$event
        time.list[[i]]  <- tdat.i$time
        ##
        tdat.i$event <- NULL
        tdat.i$time  <- NULL
        ##
        if (!null.data)
          adddata[[i]] <- tdat.i
        ##
        tdat <- tdat[!sel.i, ]
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
                         .order = seq_along(studlab),
                         stringsAsFactors = FALSE)
      ##
      if (!is.null(n))
        tdat$n <- n
      ##
      if (!is.null(event))
        tdat$event <- event
      ##
      if (!null.data) {
        tdat <- cbind(tdat, data)
        dupl <- duplicated(names(tdat))
        if (any(dupl))
          names(tdat)[dupl] <- paste(names(tdat)[dupl], "orig", sep = ".")
      }
      ##
      studlab <- names(n.arms)
      dat.studlab <- data.frame(studlab, stringsAsFactors = FALSE)
      ##
      for (i in 1:max.arms) {
        sel.i <- !duplicated(tdat$studlab)
        tdat.i <- merge(dat.studlab, tdat[sel.i, ],
                        by = "studlab", all.x = TRUE)
        ##
        treat.list[[i]] <- tdat.i$treat
        TE.list[[i]]    <- tdat.i$TE
        seTE.list[[i]]  <- tdat.i$seTE
        ##
        tdat.i$TE   <- NULL
        tdat.i$seTE <- NULL
        ##
        if (!null.data)
          adddata[[i]] <- tdat.i
        ##
        tdat <- tdat[!sel.i, ]
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
  if (length(treat) != 2 && length(studlab) != length(unique(studlab)))
    stop("Study labels must all be distinct.")
  ##
  levs <- unique(studlab)


  narms <- length(treat)
  nstud <- length(studlab)


  ##
  ## Auxiliary functions
  ##
  sumzero <- function(x)
    sum(x[!is.na(x)] == 0)
  ##
  anytrue <- function(x)
    any(x == TRUE, na.rm = TRUE)


  ##
  ##
  ## Generate dataset with variables from original dataset
  ##
  ##
  if (!null.data & !wide.armbased) {
    names.adddata <- names(adddata[[1]])
    ##
    notunique <- matrix(NA,
                        ncol = length(names.adddata),
                        nrow = narms * (narms - 1) / 2)
    colnames(notunique) <- names.adddata
    ##
    n.ij <- 0
    ##
    for (i in 1:(narms - 1)) {
      for (j in (i + 1):narms) {
        n.ij <- n.ij + 1
        notunique[n.ij, ] <- apply(adddata[[i]] != adddata[[j]], 2, anytrue)
      }
    }
    ##
    notunique <- apply(notunique, 2, anytrue)
    ##
    for (i in 1:(narms - 1)) {
      for (j in (i + 1):narms) {
        dat.i <- adddata[[i]]
        dat.j <- adddata[[j]]
        ##
        if (any(!notunique))
          dat.ij <- dat.i[, names.adddata[!notunique]]
        else
          stop("Study label must be unique for single treatment arm.")
        ##
        for (nam in names.adddata[notunique]) {
          dat.ij[, paste(nam, 1, sep = "")] <- adddata[[i]][nam]
          dat.ij[, paste(nam, 2, sep = "")] <- adddata[[j]][nam]
        }
        ##
        if (i == 1 & j == 2)
          newdata <- dat.ij
        else
          newdata <- rbind(newdata, dat.ij)
      }
    }
    ##
    names.basic <- c("studlab", "treat1", "treat2")
    names.newdata <- names(newdata)
    ##
    newdata <- newdata[, c(names.basic,
                           names.newdata[!(names.newdata %in% names.basic)])]
    newdata <- newdata[!is.na(newdata$treat1) & !is.na(newdata$treat2), ]
  }





  if (type == "binary") {
    ##
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
                          studlab,
                          treat1 = treat[[i]], treat2 = treat[[j]],
                          event1 = event[[i]], n1 = n[[i]],
                          event2 = event[[j]], n2 = n[[j]],
                          .order = seq_along(studlab),
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
          ##
          dat$TE   <- m1$TE
          dat$seTE <- m1$seTE
          ##
          dat$TE[is.infinite(dat$TE)] <- NA
          dat$seTE[is.infinite(dat$seTE)] <- NA
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
                          studlab,
                          treat1 = treat[[i]],
                          treat2 = treat[[j]],
                          n1 = n[[i]], mean1 = mean[[i]], sd1 = sd[[i]],
                          n2 = n[[j]], mean2 = mean[[j]], sd2 = sd[[j]],
                          .order = seq_along(studlab),
                          stringsAsFactors = FALSE)
        ##
        if (wide.armbased) {
          dat <- cbind(dat, data, stringsAsFactors = FALSE)
          dupl <- duplicated(names(dat))
          if (any(dupl))
            names(dat)[dupl] <- paste(names(dat)[dupl], "orig", sep = ".")
        }
        ##
        dat <- dat[!(is.na(dat$n1) & is.na(dat$mean1) & is.na(dat$sd1)), ]
        dat <- dat[!(is.na(dat$n2) & is.na(dat$mean2) & is.na(dat$sd2)), ]
        ##
        if (nrow(dat) > 0) {
          m1 <- metacont(dat$n1, dat$mean1, dat$sd1,
                         dat$n2, dat$mean2, dat$sd2,
                         ...)
          ##
          dat$TE   <- m1$TE
          dat$seTE <- m1$seTE
          ##
          dat$TE[is.infinite(dat$TE)] <- NA
          dat$seTE[is.infinite(dat$seTE)] <- NA
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
            stop("No studies available for comparison of first and second treatment.",
                 call. = FALSE)
      }
    }
  }


  if (type == "generic") {
    if (length(TE) != narms)
      stop("Different length of lists 'treat' and 'TE'.",
           call. = FALSE)
    if (length(seTE) != narms)
      stop("Different length of lists 'treat' and 'seTE'.",
           call. = FALSE)
    ##
    for (i in 1:(narms - 1)) {
      ##
      if (i == 1 & (length(treat[[i]]) != length(TE[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'TE'.",
             call. = FALSE)
      if (i == 1 & (length(treat[[i]]) != length(seTE[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'seTE'.",
             call. = FALSE)
      ##
      for (j in (i + 1):narms) {
        ##
        if (length(treat[[j]]) != length(TE[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'TE'.",
               call. = FALSE)
        if (length(treat[[j]]) != length(seTE[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'seTE'.",
               call. = FALSE)
        ##
        dat <- data.frame(TE = NA, seTE = NA,
                          studlab,
                          treat1 = treat[[i]],
                          treat2 = treat[[j]],
                          TE1 = TE[[i]], seTE1 = seTE[[i]],
                          TE2 = TE[[j]], seTE2 = seTE[[j]],
                          .order = seq_along(studlab),
                          stringsAsFactors = FALSE)
        ##
        if (!is.null(event)) {
          dat$event1 <- event[[i]]
          dat$event2 <- event[[j]]
        }
        ##
        if (!is.null(n)) {
          dat$n1 <- n[[i]]
          dat$n2 <- n[[j]]
        }
        ##
        if (!is.null(mean)) {
          dat$mean1 <- mean[[i]]
          dat$mean2 <- mean[[j]]
        }
        ##
        if (!is.null(sd)) {
          dat$sd1 <- sd[[i]]
          dat$sd2 <- sd[[j]]
        }
        ##
        if (!is.null(time)) {
          dat$time1 <- time[[i]]
          dat$time2 <- time[[j]]
        }
        ##
        if (wide.armbased) {
          dat <- cbind(dat, data, stringsAsFactors = FALSE)
          dupl <- duplicated(names(dat))
          if (any(dupl))
            names(dat)[dupl] <- paste(names(dat)[dupl], "orig", sep = ".")
        }
        ##
        dat <- dat[!(is.na(dat$TE1) & is.na(dat$seTE1)), ]
        dat <- dat[!(is.na(dat$TE2) & is.na(dat$seTE2)), ]
        ##
        if (nrow(dat) > 0) {
          m1 <- metagen(dat$TE1 - dat$TE2,
                        sqrt(dat$seTE1^2 + dat$seTE2^2), ...)
          ##
          dat$TE <- m1$TE
          dat$seTE <- m1$seTE
          ##
          dat$TE[is.infinite(dat$TE)] <- NA
          dat$seTE[is.infinite(dat$seTE)] <- NA
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
            stop("No studies available for comparison of first and second treatment.",
                 call. = FALSE)
      }
    }
  }


  if (type == "count") {
    if (length(event) != narms)
      stop("Different length of lists 'treat' and 'event'.",
           call. = FALSE)
    if (length(time) != narms)
      stop("Different length of lists 'treat' and 'time'.",
           call. = FALSE)
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
        stop("Different length of element ", i, " of lists 'treat' and 'event'.",
             call. = FALSE)
      if (i == 1 & (length(treat[[i]]) != length(time[[i]])))
        stop("Different length of element ", i, " of lists 'treat' and 'time'.",
             call. = FALSE)
      ##
      for (j in (i + 1):narms) {
        ##
        if (length(treat[[j]]) != length(event[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'event'.",
               call. = FALSE)
        if (length(treat[[j]]) != length(time[[j]]))
          stop("Different length of element ", j, " of lists 'treat' and 'time'.",
               call. = FALSE)
        ##
        dat <- data.frame(TE = NA, seTE = NA,
                          studlab,
                          treat1 = treat[[i]],
                          treat2 = treat[[j]],
                          event1 = event[[i]], time1 = time[[i]],
                          event2 = event[[j]], time2 = time[[j]],
                          incr = incr.study,
                          .order = seq_along(studlab),
                          stringsAsFactors = FALSE)
        ##
        if (wide.armbased) {
          dat <- cbind(dat, data, stringsAsFactors = FALSE)
          dupl <- duplicated(names(dat))
          if (any(dupl))
            names(dat)[dupl] <- paste(names(dat)[dupl], "orig", sep = ".")
        }
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
          ##
          dat$TE <- m1$TE
          dat$seTE <- m1$seTE
          ##
          dat$TE[is.infinite(dat$TE)] <- NA
          dat$seTE[is.infinite(dat$seTE)] <- NA
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
            stop("No studies available for comparison of first and second treatment.",
                 call. = FALSE)
      }
    }
  }
  ##
  if (!null.data & !wide.armbased)
    res <- merge(res, newdata,
                 by = c("studlab", "treat1", "treat2"),
                 suffixes = c("",".orig"),
                 all.x = TRUE)


  ##
  ## Additional checks
  ##
  ##
  ## a) Duplicate treatments ?
  ##
  sel.treat <- as.character(res$treat1) == as.character(res$treat2)
  ##
  if (any(sel.treat)) {
    sel.stud <- unique(sort(res$studlab[sel.treat]))
    ##
    stop(paste0("Identical treatments for the following stud",
                if (length(sel.stud) == 1) "y: " else "ies:\n  ",
                paste0(paste0("'", sel.stud, "'"),
                       collapse = " - "),
                "\n  Please check dataset."),
         call. = FALSE)
  }
  ##
  ## b) Studies missing ?
  ##
  sel.study <- !(studlab %in% unique(as.character(res$studlab)))
  ##
  if (any(sel.study) & warn)
    warning(paste0("The following stud",
                   if (sum(sel.study) == 1) "y is " else "ies are ",
                   "excluded from the analysis\n  ",
                   "(due to a single study arm or missing values):",
                   if (sum(sel.study) == 1) " " else "\n  ",
                   paste0(paste0("'", studlab[sel.study], "'"),
                          collapse = " - ")),
            call. = FALSE)
  ##
  ## c) Missing treatment estimates or standard errors?
  ##
  if (nrow(res.NAs) > 0 & warn) {
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
    res.NAs$.order <- NULL
    res.NAs$.order1 <- NULL
    res.NAs$.order2 <- NULL
    ##
    prmatrix(res.NAs,
             quote = FALSE, right = TRUE, na.print = "NA",
             rowlab = rep("", nrow(res.NAs)))
  }


  attr(res, "sm") <- m1$sm
  attr(res, "method") <- m1$method
  attr(res, "version") <- packageDescription("netmeta")$Version


  if (!is.null(res$.order1)) {
    res <- res[order(res$.order1), ]
    res$.order1 <- NULL
    res$.order2 <- NULL
    res$.order <- NULL
    res <- unique(res)
  }
  else if (!is.null(res$.order)) {
    res <- res[order(res$.order), ]
    res$.order <- NULL
    res$.order.orig <- NULL
    res <- unique(res)
  }
  else {
    res <- res[order(factor(res$studlab, levels = levs),
                     res$treat1, res$treat2), ]
  }
  ##
  rownames(res) <- 1:nrow(res)
  class(res) <- c(class(res), "pairwise")
  res
}
