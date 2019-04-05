##
##
## Definition of auxiliary functions for forest plots
##
##


setLab <- function(lab, col, match, value) {
  idx <- charmatch(tolower(col), tolower(match), nomatch = NA)
  sel <- !is.na(idx) & idx == 1
  if (any(sel) & length(lab) == length(col))
    lab[sel] <- value
  ##
  lab
}


setCol <- function(col, lab, match) {
  idx <- charmatch(tolower(col), tolower(match), nomatch = NA)
  sel <- !is.na(idx) & idx == 1
  if (any(sel))
    col[sel] <- match
  ##
  col
}


matchVar <- function(x, match) {
  idx <- charmatch(tolower(x), tolower(match), nomatch = NA)
  sel <- !is.na(idx) & idx == 1
  ##
  sel
}  
