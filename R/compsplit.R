compsplit <- function(x, split) {
  
  if (split %in% .special.characters)
    split <- paste("\\", split, sep = "")

  res <- strsplit(x, split)
  
  if (is.list(res))
    res <- lapply(res, gsub, pattern = "^\\s+|\\s+$", replacement = "")
  else
    res <- gsub("^\\s+|\\s+$", "", res)
  
  res
}
