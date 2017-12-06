compmatch <- function(x, split) {
  
  if (split %in% .special.characters)
    split <- paste("\\", split, sep = "")
  
  res <- any(grepl(split, x))
  
  res
}
