sortprop <- function(prop, dat.trts, sep.trts) {
  res <- rep_len(NA, length(prop))
  seq.comps <- names(prop)
  trts <-
    matrix(unlist(compsplit(seq.comps, sep.trts)),
           ncol = 2, byrow = TRUE)
  trts <- as.data.frame(trts)
  names(trts) <- c("treat1", "treat2")
  ##
  for (i in seq_len(nrow(dat.trts))) {
    sel.i <-
      (trts$treat1 == dat.trts$treat1[i] &
       trts$treat2 == dat.trts$treat2[i]) |
      (trts$treat1 == dat.trts$treat2[i] &
       trts$treat2 == dat.trts$treat1[i])
    ##
    res[i] <- prop[sel.i]
  }
  ##
  res
}
