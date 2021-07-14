setref <- function(reference.group, levs, text = "reference.group") {
  if (length(reference.group) != 1)
    stop("Argument '", text,
         "' must be a numeric or a character string.",
         call. = FALSE)
  ##
  if (is.numeric(reference.group)) {
    if (is.na(reference.group))
      stop("Missing value not allowed in argument '", text, "'.",
           call. = FALSE)
    if (!(reference.group %in% (1:length(levs))))
      stop(paste("Argument '", text,
                 "' must be any of the integers from 1 to ",
                 length(levs), ".", sep = ""), call. = FALSE)
    res <- levs[reference.group]
  }
  else if (is.character(reference.group)) {
    if (length(unique(levs)) == length(unique(tolower(levs))))
      idx <- charmatch(tolower(reference.group), tolower(levs), nomatch = NA)
    else
      idx <- charmatch(reference.group, levs, nomatch=NA)
    if (anyNA(idx) || any(idx == 0))
      stop(paste("Argument '", text,
                 "' must be any of following values:\n  ",
                 paste(paste("'", levs, "'", sep = ""),
                       collapse = " - "), sep = ""), call. = FALSE)
    res <- levs[idx]
  }
  
  res
}
