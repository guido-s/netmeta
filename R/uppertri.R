uppertri <- function(x) {
  x <- t(x)
  x[lower.tri(x)]
}
