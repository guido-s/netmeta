vikor_internal <- function(x, weights, v) {
  
  n.outcomes <- ncol(x)
  n.treatments <- nrow(x)
  #
  if (is.null(weights)) {
    # Assume equal outcome weights if argument 'weights' is NULL
    weights <- rep(1 / n.outcomes, n.outcomes)
  }
  else {
    if (length(weights) != n.outcomes)
      stop("Number of weights is different from the number of outcomes.",
           call. = FALSE)
    #
    if (!is_zero(sum(weights) - 1)) {
      weights.orig <- weights
      weights <- weights / sum(weights)
      #
      warning("Weights must sum up to 1. Accordingly, the given weights (",
              paste(weights.orig, collapse = ", "), ") ",
              "are normalized to sum to 1: (",
              paste(round(weights, 3), collapse = ", "), ").",
              call. = FALSE)
    }
  }
  #
  dist <- 1 - x
  #
  wnm <- t(t(dist) * weights)
  #
  Q <- R <- S <- vector("numeric", n.treatments)
  #
  R <- apply(wnm, 1, max)
  S <- apply(wnm, 1, sum, na.rm = TRUE)
  #
  min.R <- min(R, na.rm = TRUE)
  max.R <- max(R, na.rm = TRUE)
  #
  min.S <- min(S, na.rm = TRUE)
  max.S <- max(S, na.rm = TRUE)
  #
  Q <-
    v * (S - min.S) / (max.S - min.S) + (1 - v) * (R - min.R) / (max.R - min.R)
  #
  res <- data.frame(Q, S, R, row.names = row.names(x)) %>% arrange(Q)
  class(res) <- c("vikor", class(res))
  #
  attr(res, "weights") <- weights
  #
  res
}
