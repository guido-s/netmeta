.onAttach <-
function (libname, pkgname) 
{
  msg <- paste("Loading 'netmeta' package (version ",
               packageDescription("netmeta")$Version,
               ").",
               "\nType 'help(\"netmeta-package\")' for a brief overview.",
               sep = "")
  packageStartupMessage(msg)
}


.special.characters <- c("+", ".", "&", "$", "#", "|", "*", "^")


is.zero <- function(x, n = 10)
  abs(x) < n * .Machine$double.eps


invmat <- function(X) {
  n <- nrow(X)
  J <- matrix(1, nrow = n, ncol = n)
  ##
  res <- solve(X - J / n) + J / n
  ##
  res
}


calcV <- function(x, sm) {
  p2 <- (x$event2[1] + x$incr[1]) / (x$n2[1] + 2 * x$incr[1])
  n2 <- x$n2[1] + 2 * x$incr[1]
  ##
  if (sm == "OR")
    V <- matrix(1 / (x$event2[1] + x$incr) +
                1 / (x$n2[1] - x$event2[1] + x$incr[1]),
                nrow = nrow(x), ncol = nrow(x))
  else if (sm == "RR")
    V <- matrix((1 - p2) / (n2 * p2),
                nrow = nrow(x), ncol = nrow(x))
  else if (sm == "RD")
    V <- matrix(p2 * (1 - p2) / n2,
                nrow = nrow(x), ncol = nrow(x))
  else if (sm == "ASD")
    V <- matrix(0.25 * (1 / x$n2[1]),
                nrow = nrow(x), ncol = nrow(x))
  ##
  diag(V) <- x$seTE^2
  ##
  V
}
