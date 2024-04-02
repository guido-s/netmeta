setseq <- function(seq, levs, text, equal.length = TRUE) {
  name <- deparse(substitute(seq))
  if (missing(text))
      text <- paste0("Argument '", name, "'")
  ##
  if (length(levs) != length(seq) & equal.length)
    stop("Length of argument '", name,
         "' different from number of treatments.", call. = FALSE)
  ##
  if (length(unique(seq)) != length(seq))
    stop("Values for argument '", name,
         "' must all be disparate.", call. = FALSE)
  ##
  if (is.numeric(seq)) {
    if (anyNA(seq))
      stop("Missing values not allowed in argument '",
           name, "'.", call. = FALSE)
    if (any(!(seq %in% seq_len(length(levs)))))
      stop(paste0("Argument '", name,
                  "' must be a permutation of the integers from 1 to ",
                  length(levs), "."),
           call. = FALSE)
    res <- levs[seq]
  }
  else if (is.character(seq)) {
    if (length(unique(levs)) == length(unique(tolower(levs))))
      idx <- charmatch(tolower(seq), tolower(levs), nomatch = NA)
    else
      idx <- charmatch(seq, levs, nomatch = NA)
    ##
    if (equal.length && (anyNA(idx) || any(idx == 0)))
      stop(paste0(text,
                  " must be a permutation of the following values:\n  ",
                  paste(paste0("'", levs, "'"), collapse = " - ")),
           call. = FALSE)
    res <- levs[idx]
    if (!equal.length)
      res <- res[!is.na(res)]
  }
  else
    stop("Argument '", name, "' must be either a numeric or character vector.",
         call. = FALSE)
  
  res
}
