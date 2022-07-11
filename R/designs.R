designs <- function(treat1, treat2, studlab, sep.trts = ":") {
  
  
  id <- seq_along(studlab)
  o <- order(studlab, treat1, treat2)
  ##
  if (any(o != id)) {
    treat1 <- treat1[o]
    treat2 <- treat2[o]
    studlab <- studlab[o]
    id <- id[o]
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
  
  
  dat <- data.frame(studlab, treat1, treat2, o)
  ##
  res <- merge(dat, designs, by = "studlab")
  res <- res[order(res$o), ]
  res$treat1 <- NULL
  res$treat2 <- NULL
  res$o <- NULL
  ##
  res
}
