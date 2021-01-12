designs <- function(treat1, treat2, studlab, sep.trts = ":") {
  
  
  o <- order(studlab, treat1, treat2)
  ##
  if (any(o != seq_along(studlab))) {
    treat1 <- treat1[o]
    treat2 <- treat2[o]
    studlab <- studlab[o]
  }
  
  
  studies <- unique(studlab)
  n.study <- length(unique(studies))
  ##
  designs <- data.frame(studlab = "", design = rep_len("", n.study))
  
  
  for (i in seq_len(n.study)) {
    designs$studlab[i] <- studies[i]
    designs$design[i] <-
      paste(sort(unique(c(treat1[studlab == studies[i]],
                          treat2[studlab == studies[i]]))),
            collapse = sep.trts)
  }
  
  
  dat <- data.frame(studlab, treat1, treat2)
  ##
  res <- merge(dat, designs, by = "studlab")
  ##
  res
}
