lowertri <- function(x) {
  if (is.null(x))
    return(NULL)
  ##
  x[lower.tri(x)]
}
