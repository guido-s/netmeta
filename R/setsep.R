setsep <- function(x, sep, type = "treatment",
                   argname = deparse(substitute(sep)),
                   missing.sep) {
  labels <- sort(unique(x))
  #
  if (compmatch(labels, sep)) {
    if (!missing.sep)
      warning("Separator '", sep,
              "' used in at least one ",
              type, " label. ",
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
      sep <- "-"
    else if (!compmatch(labels, "|"))
      sep <- "|"
    else if (!compmatch(labels, "*"))
      sep <- "*"
    else
      stop("All predefined separators (':', '-', '_', '/', '+', ",
           "'.', '|', '*') are used in at least one ",
           type, " label. ",
           "Please specify a different character that should be ",
           "used as separator (argument '", argname, "').",
           call. = FALSE)
  }
  #
  sep
}
