forest.netmeta <- function(x,
                           pooled = ifelse(x$comb.random, "random", "fixed"),
                           reference.group = x$reference.group,
                           leftcols = "studlab",
                           leftlabs = "Treatment",
                           rightcols = c("effect", "ci"),
                           rightlabs = NULL,
                           small.values = "good",
                           digits.Pscore = 2,
                           smlab = NULL,
                           sortvar = x$seq,
                           ...){
  
  
  meta:::chkclass(x, "netmeta")
  
  
  ipool <- charmatch(tolower(pooled), c("fixed", "random"), nomatch = NA)
  ##
  if (is.na(ipool)) 
        stop("Argument 'pooled' should be \"fixed\" or \"random\"")
  ##
  pooled <- c("fixed", "random")[ipool]
  
  if (pooled == "fixed"){
    TE   <- x$TE.fixed
    seTE <- x$seTE.fixed
    if (is.null(smlab))
      smlab <- "Fixed Effect Model"
  }
  ##
  if (pooled == "random"){
    TE   <- x$TE.random
    seTE <- x$seTE.random
    if (is.null(smlab))
      smlab <- "Random Effects Model"
  }
  
  
  Pscore <- netrank(x, small.values = small.values)$Pscore
  
  
  sortvar.c <- deparse(substitute(sortvar))
  sortvar.c <- gsub("\"", "", sortvar.c)
  ##
  idx1 <- charmatch(tolower(sortvar.c), "pscore", nomatch = NA)
  sel1 <- !is.na(idx1) & idx1 == 1
  if (any(sel1))
    sortvar <- Pscore
  ##
  idx2 <- charmatch(tolower(sortvar.c), "-pscore", nomatch = NA)
  sel2 <- !is.na(idx2) & idx2 == 1
  if (any(sel2))
    sortvar <- -Pscore
  
  
  labels <- colnames(TE)
  ##
  if (!is.null(sortvar)){
    if (is.character(sortvar)){
      seq <- setseq(sortvar, labels)
      TE <- TE[seq, seq]
      seTE <- seTE[seq, seq]
      Pscore <- Pscore[seq]
    }
    else{
      o <- order(sortvar)
      TE <- TE[o, o]
      seTE <- seTE[o, o]
      Pscore <- Pscore[o]
    }
  }
  
  
  if (reference.group == ""){
    warning("First treatment used as reference as argument 'reference.group' is unspecified.")
    reference.group <- labels[1]
  }
  else
    reference.group <- setref(reference.group, labels)
  ##
  TE.b <- TE[, colnames(TE) == reference.group]
  seTE.b <- seTE[, colnames(seTE) == reference.group]
  study <- colnames(TE)
  ##
  dat <- data.frame(TE.b, seTE.b, study)
  ##
  idx3 <- charmatch(tolower(rightcols), "pscore", nomatch = NA)
  sel3 <- !is.na(idx3) & idx3 == 1
  if (any(sel3)) {
    dat$Pscore <- meta:::format.NA(Pscore, digits = digits.Pscore)
    rightcols[sel3] <- "Pscore"
  }
  ##
  idx4 <- charmatch(tolower(leftcols), "pscore", nomatch = NA)
  sel4 <- !is.na(idx4) & idx4 == 1
  if (any(sel4)) {
    dat$Pscore <- meta:::format.NA(Pscore, digits = digits.Pscore)
    leftcols[sel4] <- "Pscore"
  }
  
  
  m1 <- metagen(TE.b, seTE.b, data = dat,
                sm = x$sm,
                studlab = colnames(TE), warn = FALSE)
  
  
  forest(m1,
         comb.fixed = FALSE, comb.random = FALSE,
         hetstat = FALSE,
         leftcols = leftcols,
         leftlabs = leftlabs,
         rightcols = rightcols,
         rightlabs = rightlabs,
         smlab = smlab,
         ...)
  
  invisible(NULL)
}
