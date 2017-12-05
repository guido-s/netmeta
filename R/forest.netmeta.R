forest.netmeta <- function(x,
                           pooled = ifelse(x$comb.random, "random", "fixed"),
                           reference.group = x$reference.group,
                           baseline.reference = x$baseline.reference,
                           leftcols = "studlab",
                           leftlabs = "Treatment",
                           rightcols = c("effect", "ci"),
                           rightlabs = NULL,
                           digits = gs("digits.forest"),
                           small.values = "good",
                           digits.Pscore = 2,
                           smlab = NULL,
                           sortvar = x$seq,
                           backtransf = x$backtransf,
                           lab.NA = ".",
                           add.data,
                           drop.reference.group = FALSE,
                           ...) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  meta:::chkclass(x, "netmeta")
  x <- upgradenetmeta(x)
  ##
  chklogical <- meta:::chklogical
  formatN <- meta:::formatN
  ##
  pooled <- meta:::setchar(pooled, c("fixed", "random"))
  ##
  meta:::chknumeric(digits, min = 0, single = TRUE)
  ##
  chklogical(baseline.reference)
  chklogical(drop.reference.group)
  ##
  chklogical(backtransf)
  meta:::chkchar(lab.NA)
  
  
  ##
  ##
  ## (2) Extract results for fixed effect and random effects model
  ##     and calculate P-scores
  ##
  ##
  labels <- colnames(x$TE.fixed)
  ##
  if (reference.group == "") {
    warning("First treatment used as reference as argument 'reference.group' is unspecified.")
    reference.group <- labels[1]
  }
  else
    reference.group <- setref(reference.group, labels)
  ##
  if (pooled == "fixed") {
    TE   <- x$TE.fixed
    seTE <- x$seTE.fixed
    prop.direct <- x$P.fixed
    if (is.null(smlab))
      if (baseline.reference)
        smlab <- paste("Comparison: other vs '",
                       reference.group, "'\n(Fixed Effect Model)",
                       sep = "")
      else
        smlab <- paste("Comparison: '",
                       reference.group, "' vs other\n(Fixed Effect Model)",
                       sep = "")
    ##
    Pscore <- netrank(x, small.values = small.values)$Pscore.fixed
  }
  ##
  if (pooled == "random") {
    TE   <- x$TE.random
    seTE <- x$seTE.random
    prop.direct <- x$P.random
    if (is.null(smlab))
      if (baseline.reference)
        smlab <- paste("Comparison: other vs '",
                       reference.group, "'\n(Random Effects Model)",
                       sep = "")
      else
        smlab <- paste("Comparison: '",
                       reference.group, "' vs other\n(Random Effects Model)",
                       sep = "")
    ##
    Pscore <- netrank(x, small.values = small.values)$Pscore.random
  }
  
  
  ##
  ##
  ## (3) Extract comparisons with reference group
  ##
  ##
  if (baseline.reference)
    dat <- data.frame(TE = TE[, colnames(TE) == reference.group],
                      seTE = seTE[, colnames(seTE) == reference.group],
                      trts = colnames(TE),
                      k = x$A.matrix[, colnames(TE) == reference.group],
                      prop.direct = prop.direct[, colnames(TE) == reference.group],
                      row.names = colnames(TE),
                      as.is = TRUE)
  else
    dat <- data.frame(TE = TE[rownames(TE) == reference.group, ],
                      seTE = seTE[rownames(seTE) == reference.group, ],
                      trts = rownames(TE),
                      k = x$A.matrix[rownames(TE) == reference.group, ],
                      prop.direct = prop.direct[rownames(TE) == reference.group, ],
                      row.names = colnames(TE),
                      as.is = TRUE)
  ##
  rm(TE)
  rm(seTE)
  ##
  idx1 <- charmatch(tolower(rightcols), "pscore", nomatch = NA)
  sel1 <- !is.na(idx1) & idx1 == 1
  if (any(sel1)) {
    dat$Pscore <- formatN(Pscore, digits = digits.Pscore,
                          text.NA = lab.NA)
    rightcols[sel1] <- "Pscore"
  }
  ##
  idx2 <- charmatch(tolower(leftcols), "pscore", nomatch = NA)
  sel2 <- !is.na(idx2) & idx2 == 1
  if (any(sel2)) {
    dat$Pscore <- formatN(Pscore, digits = digits.Pscore,
                          text.NA = lab.NA)
    leftcols[sel2] <- "Pscore"
  }
  ##
  if (!missing(add.data)) {
    if (!is.data.frame(add.data))
      stop("Argument 'add.data' must be a data frame.",
           call. = FALSE)
    if (nrow(add.data) != length(labels))
      stop("Dataset 'add.data' must have ", nrow(dat),
           " rows (corresponding to number of treatments)",
           call. = FALSE)
    if (any(rownames(add.data) != labels))
      stop("Dataset 'add.data' must have the following row names:\n",
           paste(paste("'", labels, "'", sep = ""), collapse = " - "),
           call. = FALSE)
    ##
    dat <- cbind(dat, add.data)
  }
  
  
  ##
  ##
  ## (4) Sort dataset according to argument sortvar
  ##
  ##
  sortvar.c <- deparse(substitute(sortvar))
  sortvar.c <- gsub("\"", "", sortvar.c)
  ##
  idx3 <- charmatch(tolower(sortvar.c), "pscore", nomatch = NA)
  sel3 <- !is.na(idx3) & idx3 == 1
  if (any(sel3))
    sortvar <- Pscore
  ##
  idx4 <- charmatch(tolower(sortvar.c), "-pscore", nomatch = NA)
  sel4 <- !is.na(idx4) & idx4 == 1
  if (any(sel4))
    sortvar <- -Pscore
  ##
  idx5 <- charmatch(tolower(sortvar.c), "te", nomatch = NA)
  sel5 <- !is.na(idx5) & idx5 == 1
  if (any(sel5))
    sortvar <- dat$TE
  ##
  idx6 <- charmatch(tolower(sortvar.c), "-te", nomatch = NA)
  sel6 <- !is.na(idx6) & idx6 == 1
  if (any(sel6))
    sortvar <- -dat$TE
  ##
  idx7 <- charmatch(tolower(sortvar.c), "sete", nomatch = NA)
  sel7 <- !is.na(idx7) & idx7 == 1
  if (any(sel7))
    sortvar <- dat$seTE
  ##
  idx8 <- charmatch(tolower(sortvar.c), "-sete", nomatch = NA)
  sel8 <- !is.na(idx8) & idx8 == 1
  if (any(sel8))
    sortvar <- -dat$seTE
  ##
  idx9 <- charmatch(tolower(sortvar.c), "k", nomatch = NA)
  sel9 <- !is.na(idx9) & idx9 == 1
  if (any(sel9))
    sortvar <- dat$k
  ##
  idx10 <- charmatch(tolower(sortvar.c), "-k", nomatch = NA)
  sel10 <- !is.na(idx10) & idx10 == 1
  if (any(sel10))
    sortvar <- -dat$k
  ##
  idx11 <- charmatch(tolower(sortvar.c), "prop.direct", nomatch = NA)
  sel11 <- !is.na(idx11) & idx11 == 1
  if (any(sel11))
    sortvar <- dat$prop.direct
  ##
  idx12 <- charmatch(tolower(sortvar.c), "-prop.direct", nomatch = NA)
  sel12 <- !is.na(idx12) & idx12 == 1
  if (any(sel12))
    sortvar <- -dat$prop.direct
  ##  
  if (!is.null(sortvar)) {
    if (is.character(sortvar))
      sort <- setseq(sortvar, labels)
    else
      sort <- order(sortvar)
    ##
    dat <- dat[sort, ]
  }
  
  
  ##
  ##
  ## (5) Generate forest plot
  ##
  ##
  if (drop.reference.group)
    dat <- subset(dat, trts != reference.group)
  ##
  dat$prop.direct <- formatN(dat$prop.direct,
                             digits = digits.Pscore, text.NA = lab.NA)
  ##
  trts <- dat$trts
  m1 <- metagen(TE, seTE, data = dat,
                sm = x$sm,
                studlab = trts, backtransf = backtransf,
                warn = FALSE)
  ##
  forest.meta(m1,
              digits = digits,
              comb.fixed = FALSE, comb.random = FALSE,
              hetstat = FALSE,
              leftcols = leftcols,
              leftlabs = leftlabs,
              rightcols = rightcols,
              rightlabs = rightlabs,
              smlab = smlab,
              lab.NA = lab.NA,
              ...)
  
  
  invisible(NULL)
}
