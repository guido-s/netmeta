compsplit <- function(x, split) {
  
  if (split %in% .special.characters)
    split <- paste("\\", split, sep = "")
  
  res <- strsplit(x, split)
  
  res
}
