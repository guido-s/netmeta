uppertri <- function(x) {
  if (is.null(x))
    return(NULL)
  ##
  x <- t(x)
  x[lower.tri(x)]
}
