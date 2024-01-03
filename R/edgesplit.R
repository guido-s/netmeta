edgesplit <- function(x,
                      sep.trts = " vs ", verbose = TRUE, warn = TRUE,
                      is.tictoc = FALSE) {
  
  chkclass(x, "netmeta")
  is.bin <- inherits(x, "netmetabin")
  
  if (is.null(x$data))
    stop("Node-splitting method only available for network meta-analysis ",
         "objects created with argument 'keepdata' equal to TRUE.")
  #
  if (verbose)
    cat("Start computations for node-splitting method\n")
  #
  dat <- x$data
  dat <- dat[order(dat$.studlab, dat$.treat1, dat$.treat2), ]
  #
  if (!is.null(dat$.subset))
    dat <- dat[dat$.subset, , drop = FALSE]
  #
  if (!is.null(dat$.drop))
    dat <- dat[!dat$.drop, , drop = FALSE]
  #
  # Determine comparisons with direct evidence
  #
  idx.d <- which(!is.na(x$TE.direct.common), arr.ind = TRUE)
  idx.d <- idx.d[idx.d[, 1] < idx.d[, 2], , drop = FALSE]
  #
  rownames(idx.d) <- seq_len(nrow(idx.d))
  idx1 <- idx.d[, 1]
  idx2 <- idx.d[, 2]
  #
  n.comps <- nrow(idx.d)
  #
  trts <- x$trts
  #
  # Perform network meta-analyses for indirect evidence
  # (by dropping one direct comparison at a time)
  #
  TE.indirect.common <- x$TE.direct.common
  TE.indirect.common[!is.na(TE.indirect.common)] <- NA
  seTE.indirect.common <- TE.indirect.common
  #
  TE.indirect.random <- x$TE.direct.random
  TE.indirect.random[!is.na(TE.indirect.random)] <- NA
  seTE.indirect.random <- TE.indirect.random
  #
  if (is.tictoc) {
    tictoc <- rep(NA, n.comps)
    names.tictoc <- ""
  }
  #
  for (i in seq_len(n.comps)) {
    #
    idx1.i <- idx1[i]
    idx2.i <- idx2[i]
    #
    if (is.tictoc)
      tictoc::tic()
    #
    if (verbose)
      cat(paste0("- ",
                 paste(trts[idx1.i], trts[idx2.i], sep = sep.trts),
                 " (", i, "/", n.comps, ")\n"))
    #
    # Determine all pairwise comparisons of trts[idx1.i] vs trts[idx2.i]
    #
    direct.i <-
      (dat$.treat1 == trts[idx1.i] & dat$.treat2 == trts[idx2.i]) |
      (dat$.treat2 == trts[idx1.i] & dat$.treat1 == trts[idx2.i])
    #
    # Determine all studies with pairwise comparison trts[idx1.i] vs
    # trts[idx2.i]
    #
    study.direct.i <- unique(dat$.studlab[direct.i])
    #
    # Only drop pairwise comparisons of trts[idx1.i] vs trts[idx2.i],
    # i.e., keep other comparisons of a multi-arm
    # study
    #
    dat.i <- dat[!direct.i, , drop = FALSE]
    #
    # All study labels must be different as adjusted standard errors
    # are used to calculate indirect estimates
    #
    dat.i$.studlab <- seq_len(nrow(dat.i))
    #
    dat.i$.design <- NULL
    #
    if (nrow(dat.i) > 0)
      con <- netconnection(dat.i$.treat1, dat.i$.treat2, dat.i$.studlab)
    else
      con <- list(n.subnets = 0)
    #
    if (con$n.subnets == 1) {
      #
      if (is.bin)
        net.i <- netmetabin(dat.i$.event1, dat.i$.n1,
                            dat.i$.event2, dat.i$.n2,
                            dat.i$.treat1, dat.i$.treat2,
                            dat.i$.studlab,
                            data = dat.i,
                            sm = x$sm, method = x$method,
                            warn = warn)
      else {
        net.i <- netmeta(dat.i$.TE, dat.i$.seTE.adj.common,
                         dat.i$.treat1, dat.i$.treat2,
                         dat.i$.studlab,
                         data = dat.i,
                         method.tau = x$method.tau,
                         warn = warn)
        #
        net.r.i <- netmeta(dat.i$.TE, dat.i$.seTE.adj.random,
                           dat.i$.treat1, dat.i$.treat2,
                           dat.i$.studlab,
                           data = dat.i,
                           method.tau = x$method.tau,
                           warn = warn)
      }
      #
      if (trts[idx1.i] %in% rownames(net.i$TE.common) &
          trts[idx2.i] %in% colnames(net.i$TE.common)) {
        TE.indirect.common[idx1.i, idx2.i] <-
          net.i$TE.common[trts[idx1.i], trts[idx2.i]]
        TE.indirect.common[idx2.i, idx1.i] <-
          net.i$TE.common[trts[idx2.i], trts[idx1.i]]
        #
        seTE.indirect.common[idx1.i, idx2.i] <-
          seTE.indirect.common[idx2.i, idx1.i] <-
          net.i$seTE.common[trts[idx1.i], trts[idx2.i]]
      }
      #
      if (!is.bin) {
        if (trts[idx1.i] %in% rownames(net.r.i$TE.random) &
            trts[idx2.i] %in% colnames(net.r.i$TE.random)) {
          TE.indirect.random[idx1.i, idx2.i] <-
            net.r.i$TE.common[trts[idx1.i], trts[idx2.i]]
          TE.indirect.random[idx2.i, idx1.i] <-
            net.r.i$TE.common[trts[idx2.i], trts[idx1.i]]
          #
          seTE.indirect.random[idx1.i, idx2.i] <-
            seTE.indirect.random[idx2.i, idx1.i] <-
            net.r.i$seTE.common[trts[idx1.i], trts[idx2.i]]
        }
      }
    }
    #
    if (is.tictoc) {
      tictoc.i <- tictoc::toc(func.toc = NULL)
      tictoc[i] <- as.numeric(tictoc.i$toc) - as.numeric(tictoc.i$tic)
      names.tictoc[i] <- paste(trts[idx1.i], trts[idx2.i], sep = sep.trts)
      #
      if (verbose)
        cat(paste(round(tictoc[i], 3), "sec elapsed\n"))
    }
  }
  
  
  res <- list(TE.indirect.common = TE.indirect.common,
              seTE.indirect.common = seTE.indirect.common)
  #
  if (!is.bin) {
    res$TE.indirect.random <- TE.indirect.random
    res$seTE.indirect.random <- seTE.indirect.random
  }
  #
  if (is.tictoc) {
    names(tictoc) <- names.tictoc
    res$tictoc <- tictoc
  }
  #
  res
}
