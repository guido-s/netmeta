treats <- function(x, nchar.trts = 8, row = TRUE) {
  
  meta:::chknumeric(nchar.trts, min = 1, single = TRUE)
  meta:::chklogical(row)
  
  if (is.matrix(x)) {
    if (row)
      trts <- rownames(x)
    else
      trts <- colnames(x)
  }
  else
    trts <- x
  ##
  ## Default: first 'nchar.trts' character of treatment names
  ##
  res <- substring(trts, 1, nchar.trts)
  ##
  ## Use abbreviated treatment names if necessary
  ##
  if (length(unique(res)) != length(unique(trts)))
    res <- as.character(abbreviate(trts, nchar.trts))
  ##
  res
}
