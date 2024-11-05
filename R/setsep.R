setsep <- function(x, sep, type = "treatment label",
                   argname = deparse(substitute(sep)),
                   missing = TRUE) {
  labels <- sort(unique(x))
  #
  if (compmatch(labels, sep)) {
    if (!missing)
      warning("Argument '", argname, "': ",
              "separator '", sep, "' used in at least one ", type, ". ",
              "Trying to use predefined separators: ",
              "':', '-', '_', '/', '+', '.', '|', '*'.",
              call. = FALSE)
    #
    if (!compmatch(labels, ":"))
      sep <- ":"
    else if (!compmatch(labels, "-"))
      sep <- "-"
    else if (!compmatch(labels, "_"))
      sep <- "_"
    else if (!compmatch(labels, "/"))
      sep <- "/"
    else if (!compmatch(labels, "+"))
      sep <- "+"
    else if (!compmatch(labels, "."))
      sep <- "."
    else if (!compmatch(labels, "|"))
      sep <- "|"
    else if (!compmatch(labels, "*"))
      sep <- "*"
    else
      stop("All predefined separators (':', '-', '_', '/', '+', ",
           "'.', '|', '*') are used in at least one ", type, ". ",
           "Please specify a different character that should be ",
           "used as separator (argument '", argname, "').",
           call. = FALSE)
  }
  #
  sep
}
