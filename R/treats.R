treats <- function(x, nchar, row = TRUE) {
  if (row)
    res <- as.character(abbreviate(rownames(x), nchar))
  else
    res <- as.character(abbreviate(colnames(x), nchar))
  ##
  res
}
