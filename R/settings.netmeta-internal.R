catarg <- function(x, newline = TRUE, end = "") {
  xname <- x
  x <- gsub(" ", "", x)
  ##
  if (newline)
    cat("- ")
  ##
  if (is.null(gs(x)))
    cat(paste0(xname, ' = NULL', end, '\n'))
  else if (is.character(gs(x)))
    cat(paste0(xname, ' = "', gs(x), '"', end, '\n'))
  else
    cat(paste0(xname, ' = ', gs(x), end, "\n"))
  ##
  invisible(NULL)
}


setcharacter <- function(argname, args, set = NULL, length = 1,
                         NULL.ok = FALSE, ignore.other = FALSE,
                         logical.ok = FALSE) {
  id <- argid(names(args), argname)
  ##
  if (!is.na(id)) {
    val <- args[[id]]
    ##
    if (NULL.ok & is.null(val)) {
      l <- list(val)
      names(l) <- argname
      settings.meta(l)
      return(invisible(NULL))
    }
    ##
    if (logical.ok & is.logical(val)) {
      l <- list(val)
      names(l) <- argname
      settings.meta(l)
      return(invisible(NULL))
    }
    ##
    if (!is.character(val) & ignore.other)
      return(invisible(id))
    ##
    if (!is.null(set))
      val <- setchar(val, set, name = argname)
    else
      chkchar(val, length = length, name = argname)
    ##
    l <- list(val)
    names(l) <- argname
    settings.meta(l)
  }
  ##
  invisible(id)
}


setcolor <- function(argname, args) {
  id <- argid(names(args), argname)
  ##
  if (!is.na(id)) {
    val <- args[[id]]
    chkcolor(val, name = argname)
    #
    l <- list(val)
    names(l) <- argname
    settings.meta(l)
  }
  ##
  invisible(id)
}


setlevel <- function(argname, args) {
  id <- argid(names(args), argname)
  ##
  if (!is.na(id)) {
    val <- args[[id]]
    chklevel(val, name = argname)
    #
    l <- list(val)
    names(l) <- argname
    settings.meta(l)
  }
  ##
  invisible(id)
}


setlogical <- function(argname, args, NULL.ok = FALSE,
                       ignore.other = FALSE) {
  id <- argid(names(args), argname)
  ##
  if (!is.na(id)) {
    val <- args[[id]]
    ##
    if (NULL.ok & is.null(val)) {
      l <- list(val)
      names(l) <- argname
      settings.meta(l)
      return(invisible(NULL))
    }
    ##
    if (!is.logical(val) & ignore.other)
      return(invisible(id))
    ##
    chklogical(val, name = argname)
    #
    l <- list(val)
    names(l) <- argname
    settings.meta(l)
  }
  ##
  invisible(id)
}


setnumeric <- function(argname, args, NULL.ok = FALSE) {
  id <- argid(names(args), argname)
  ##
  if (!is.na(id)) {
    val <- args[[id]]
    ##
    if (NULL.ok & is.null(val)) {
      l <- list(val)
      names(l) <- argname
      settings.meta(l)
      return(invisible(NULL))
    }
    ##
    chknumeric(val, min = 0, length = 1, name = argname)
    #
    l <- list(val)
    names(l) <- argname
    settings.meta(l)
  }
  ##
  invisible(id)
}
