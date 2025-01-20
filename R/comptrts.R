comptrts <- function(x,
                     upper = TRUE, reference.group = "", baseline.reference,
                     order = NULL, sep.trts = x$sep.trts, quote.trts = "") {
  
  
  seq.comps <- rownames(x$Cov.common)
  ##
  dat.trts <- matrix(unlist(compsplit(seq.comps, x$sep.trts)),
                     ncol = 2, byrow = TRUE)
  dat.trts <- as.data.frame(dat.trts)
  names(dat.trts) <- c("treat1", "treat2")
  ##
  if (!upper) {
    ##
    ## Comparison names are column:row (and must be switched)
    ##
    t1 <- dat.trts$treat1
    dat.trts$treat1 <- dat.trts$treat2
    dat.trts$treat2 <- t1
  }
  ##  
  if (is.null(order)) {
    ##
    ## Change treatment order if
    ## - reference group is specified, i.e., unequal to ""
    ## - reference group is first treatment
    ##   (argument 'baseline.reference' is TRUE)
    ## - reference group is second treatment
    ##   (argument 'baseline.reference' is FALSE)
    ##
    wo <- rep_len(FALSE, length(seq.comps))
    ##
    if (reference.group != "") {
      reference.group <- setref(reference.group, colnames(x$TE.common))
      ##
      if (baseline.reference)
        wo <- dat.trts$treat1 == reference.group
      else
        wo <- dat.trts$treat2 == reference.group
    }
    ##
    if (any(wo)) {
      t1.wo <- dat.trts$treat1[wo]
      dat.trts$treat1[wo] <- dat.trts$treat2[wo]
      dat.trts$treat2[wo] <- t1.wo
    }
  }
  else {
    treat1.pos <- as.numeric(factor(dat.trts$treat1, levels = order))
    treat2.pos <- as.numeric(factor(dat.trts$treat2, levels = order))
    ##
    wo <- treat1.pos > treat2.pos
    ##
    if (any(wo)) {
      ttreat1 <- dat.trts$treat1
      dat.trts$treat1[wo] <- dat.trts$treat2[wo]
      dat.trts$treat2[wo] <- ttreat1[wo]
      ##
      ttreat1.pos <- treat1.pos
      treat1.pos[wo] <- treat2.pos[wo]
      treat2.pos[wo] <- ttreat1.pos[wo]
    }
    ##
    o <- order(treat1.pos, treat2.pos)
    dat.trts <- dat.trts[o, ]
  }
  
  
  comparison <-
    as.character(
      interaction(
        paste0(quote.trts, dat.trts$treat1, quote.trts),
        paste0(quote.trts, dat.trts$treat2, quote.trts),
        sep = sep.trts))
  
  
  res <- cbind(comparison, dat.trts)
  res
}

