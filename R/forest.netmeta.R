forest.netmeta <- function(x,
                           pooled = ifelse(x$comb.random, "random", "fixed"),
                           reference.group = x$reference.group,
                           baseline.reference = x$baseline.reference,
                           leftcols = "studlab",
                           leftlabs,
                           rightcols = c("effect", "ci"),
                           rightlabs,
                           digits = gs("digits.forest"),
                           small.values = "good",
                           digits.Pscore = 2,
                           smlab = NULL,
                           sortvar = x$seq,
                           backtransf = x$backtransf,
                           lab.NA = ".",
                           add.data,
                           drop.reference.group = FALSE,
                           ##
                           col.by = "black",
                           print.byvar = FALSE,
                           ##
                           ...) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  meta:::chkclass(x, "netmeta")
  x <- upgradenetmeta(x)
  ##
  is.bin <- inherits(x, "netmetabin")
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
  chklogical(print.byvar)
  ##
  chklogical(backtransf)
  meta:::chkchar(lab.NA)
  ##
  if (missing(leftlabs))
    if (length(reference.group) > 1)
      leftlabs <- "Comparison"
    else
      leftlabs <- "Treatment"
  
  
  ##
  ##
  ## (2) Extract results for fixed effect and random effects model
  ##     and calculate P-scores
  ##
  ##
  labels <- colnames(x$TE.fixed)
  ##
  one.rg <- length(reference.group) == 1
  ##
  if (one.rg) {
    if (reference.group == "") {
      warning("First treatment used as reference as argument ",
              "'reference.group' is unspecified.",
              call. = FALSE)
      reference.group <- labels[1]
    }
    else
      reference.group <- setref(reference.group, labels)
  }
  else
    for (i in seq_along(reference.group))
      reference.group[i] <- setref(reference.group[i], labels)
  ##
  if (pooled == "fixed") {
    TE   <- x$TE.fixed
    seTE <- x$seTE.fixed
    ##
    prop.direct <- x$P.fixed
    ##
    Pscore <- netrank(x, small.values = small.values)$Pscore.fixed
    ##
    text.pooled <- "Fixed Effect Model"
    ##
    if (x$method == "MH")
      text.pooled <- "Mantel-Haenszel Method"
    else if (x$method == "NCH")
      text.pooled <- "Non-Central Hypergeometric"
  }
  ##
  if (pooled == "random") {
    TE   <- x$TE.random
    seTE <- x$seTE.random
    ##
    prop.direct <- x$P.random
    ##
    Pscore <- netrank(x, small.values = small.values)$Pscore.random
    ##
    text.pooled <- "Random Effects Model"
  }
  ##
  if (is.null(smlab)) {
    if (one.rg) {
      if (baseline.reference)
        smlab <- paste0("Comparison: other vs '",
                        reference.group, "'\n(",
                        text.pooled,
                        ")")
      else
        smlab <- paste0("Comparison: '",
                        reference.group,
                        "' vs other \n(",
                        text.pooled,
                        ")")
    }
    else
      smlab  <- text.pooled
  }
  ##
  idx1 <- charmatch(tolower(rightcols), "pscore", nomatch = NA)
  sel1 <- !is.na(idx1) & idx1 == 1
  if (any(sel1))
    rightcols[sel1] <- "Pscore"
  ##
  idx2 <- charmatch(tolower(leftcols), "pscore", nomatch = NA)
  sel2 <- !is.na(idx2) & idx2 == 1
  if (any(sel2))
    leftcols[sel2] <- "Pscore"
  ##
  sortvar.c <- deparse(substitute(sortvar))
  sortvar.c <- gsub("\"", "", sortvar.c)
  ##
  idx3 <- charmatch(tolower(sortvar.c), "pscore", nomatch = NA)
  sel3 <- !is.na(idx3) & idx3 == 1
  idx4 <- charmatch(tolower(sortvar.c), "-pscore", nomatch = NA)
  sel4 <- !is.na(idx4) & idx4 == 1
  ##
  idx5 <- charmatch(tolower(sortvar.c), "te", nomatch = NA)
  sel5 <- !is.na(idx5) & idx5 == 1
  idx6 <- charmatch(tolower(sortvar.c), "-te", nomatch = NA)
  sel6 <- !is.na(idx6) & idx6 == 1
  ##
  idx7 <- charmatch(tolower(sortvar.c), "sete", nomatch = NA)
  sel7 <- !is.na(idx7) & idx7 == 1
  idx8 <- charmatch(tolower(sortvar.c), "-sete", nomatch = NA)
  sel8 <- !is.na(idx8) & idx8 == 1
  ##
  idx9 <- charmatch(tolower(sortvar.c), "k", nomatch = NA)
  sel9 <- !is.na(idx9) & idx9 == 1
  idx10 <- charmatch(tolower(sortvar.c), "-k", nomatch = NA)
  sel10 <- !is.na(idx10) & idx10 == 1
  ##
  idx11 <- charmatch(tolower(sortvar.c), "prop.direct", nomatch = NA)
  sel11 <- !is.na(idx11) & idx11 == 1
  idx12 <- charmatch(tolower(sortvar.c), "-prop.direct", nomatch = NA)
  sel12 <- !is.na(idx12) & idx12 == 1
  
  
  ##
  ##
  ## (3) Extract comparisons with reference group
  ##
  ##
  dat <- data.frame(comparison = character(0),
                    treat = character(0),
                    TE = numeric(0), seTE = numeric(0),
                    Pscore = numeric(0),
                    k = numeric(0),
                    prop.direct = numeric(0),
                    stringsAsFactors = FALSE)
  ##
  for (i in seq_along(reference.group)) {
    rg.i <- reference.group[i]
    ##
    if (baseline.reference)
      dat.i <- data.frame(comparison = rg.i,
                          treat = colnames(TE),
                          TE = TE[, colnames(TE) == rg.i],
                          seTE = seTE[, colnames(seTE) == rg.i],
                          Pscore = Pscore,
                          k = x$A.matrix[, colnames(TE) == rg.i],
                          prop.direct =
                            if (is.bin) prop.direct
                            else prop.direct[, colnames(TE) == rg.i],
                          stringsAsFactors = FALSE)
    else
      dat.i <- data.frame(comparison = rg.i,
                          Pscore = Pscore,
                          TE = TE[rownames(TE) == rg.i, ],
                          seTE = seTE[rownames(seTE) == rg.i, ],
                          treat = rownames(TE),
                          k = x$A.matrix[rownames(TE) == rg.i, ],
                          prop.direct =
                            if (is.bin) prop.direct
                            else prop.direct[rownames(TE) == rg.i, ],
                          stringsAsFactors = FALSE)
    ##
    if (!missing(add.data)) {
      if (!is.data.frame(add.data))
        stop("Argument 'add.data' must be a data frame.",
             call. = FALSE)
      if (nrow(add.data) != length(labels))
        stop("Dataset 'add.data' must have ", nrow(dat.i),
             " rows (corresponding to number of treatments)",
             call. = FALSE)
      if (any(rownames(add.data) != labels))
        stop("Dataset 'add.data' must have the following row names:\n",
             paste(paste("'", labels, "'", sep = ""), collapse = " - "),
             call. = FALSE)
      ##
      dat.i <- cbind(dat.i, add.data)
    }
    ##
    ## Sort dataset according to argument sortvar
    ##
    if (any(sel3))
      sortvar <- Pscore
    ##
    if (any(sel4))
      sortvar <- -Pscore
    ##
    if (any(sel5))
      sortvar <- dat.i$TE
    ##
    if (any(sel6))
      sortvar <- -dat.i$TE
    ##
    if (any(sel7))
      sortvar <- dat.i$seTE
    ##
    if (any(sel8))
      sortvar <- -dat.i$seTE
    ##
    if (any(sel9))
      sortvar <- dat.i$k
    ##
    if (any(sel10))
      sortvar <- -dat.i$k
    ##
    if (any(sel11))
      sortvar <- dat.i$prop.direct
    ##
    if (any(sel12))
      sortvar <- -dat.i$prop.direct
    ##
    if (!is.null(sortvar)) {
      if (is.character(sortvar))
        sort <- setseq(sortvar, labels)
      else
        sort <- order(sortvar)
      ##
      dat.i <- dat.i[sort, ]
    }
    ##
    if (drop.reference.group)
      dat.i <- subset(dat.i, treat != rg.i)
    ##
    if (baseline.reference)
      dat.i$comparison <- paste0("Other vs '", dat.i$comparison, "'")
    else
      dat.i$comparison <- paste0("'", dat.i$comparison, "' vs other")
    ##
    dat <- rbind(dat, dat.i)
  }
  ##
  dat.out <- dat
  ##
  if ("Pscore" %in% names(dat))
    dat$Pscore <- formatN(dat$Pscore, digits = digits.Pscore,
                          text.NA = lab.NA)
  ##
  if ("prop.direct" %in% names(dat))
    dat$prop.direct <- formatN(dat$prop.direct,
                               digits = digits.Pscore, text.NA = lab.NA)
  ##
  rm(TE)
  rm(seTE)
  
  
  ##
  ##
  ## (5) Generate forest plot
  ##
  ##
  treat <- dat$treat
  ##
  if (one.rg)
    m1 <- metagen(TE, seTE, data = dat,
                  sm = x$sm,
                  studlab = treat, backtransf = backtransf,
                  warn = FALSE)
  else
    m1 <- metagen(TE, seTE, data = dat,
                  byvar = dat$comparison,
                  sm = x$sm,
                  studlab = treat, backtransf = backtransf,
                  warn = FALSE)
  ##
  forest(m1,
         digits = digits,
         comb.fixed = FALSE, comb.random = FALSE,
         hetstat = FALSE,
         leftcols = leftcols,
         leftlabs = leftlabs,
         rightcols = rightcols,
         rightlabs = if (missing(rightlabs)) NULL else rightlabs,
         smlab = smlab,
         lab.NA = lab.NA,
         ##
         col.by = col.by,
         print.byvar = print.byvar,
         ##
         ...)
  

  rownames(dat.out) <- seq_len(nrow(dat.out))
  ##
  attr(dat.out, "pooled") <- pooled
  attr(dat.out, "small.values") <- small.values
  attr(dat.out, "small.values") <- small.values
  ##
  invisible(dat.out)
}
