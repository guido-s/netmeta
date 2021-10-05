compsplit <- function(x, split) {
  
  if (split %in% .special.characters)
    split <- paste("\\", split, sep = "")

  res <- strsplit(x, split)

  if (is.list(res)) {
    withspace <- any(unlist(lapply(res, grepl, pattern = "^\\s+|\\s+$")))
    res <- lapply(res, gsub, pattern = "^\\s+|\\s+$", replacement = "")
  }
  else {
    withspace <- any(grepl("^\\s+|\\s+$", "", res))
    res <- gsub("^\\s+|\\s+$", "", res)
  }
  ##
  attr(res, "withspace") <- withspace
  
  res
}
