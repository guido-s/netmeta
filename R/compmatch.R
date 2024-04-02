compmatch <- function(x, split) {
  
  if (split %in% .special.characters)
    split <- paste0("\\", split)
  
  res <- any(grepl(split, x))
  
  res
}
