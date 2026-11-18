vikor_internal <- function(x, weights, v) {
  
  # Get rid of warning 'no visible binding for global variable'
  Q <- R <- S <- NULL
  
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
    chknumeric(weights, min = 0, zero = TRUE)
    #
    if (is_zero(sum(weights)))
      stop("Sum of weights must be larger than 0.",
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
  res <- data.frame(Q = NA, S = apply(wnm, 1, sum), R = apply(wnm, 1, max),
                    row.names = row.names(x))
  #
  dropped_treatments <- rownames(res %>% filter(if_any(c(S, R), is.na)))
  #
  res %<>% drop_na(S, R)
  #
  if (nrow(res) < 2)
    stop("VIKOR method requires at least two treatments with ",
         "complete rankings.",
         call. = FALSE)
  #
  min.R <- min(res$R, na.rm = TRUE)
  max.R <- max(res$R, na.rm = TRUE)
  #
  min.S <- min(res$S, na.rm = TRUE)
  max.S <- max(res$S, na.rm = TRUE)
  #
  if (is_zero(min.R - max.R))
    stop("VIKOR method not applicable as all values of ranking statistic R ",
         "are identical.",
         call. = FALSE)
  #
  if (is_zero(min.S - max.S))
    stop("VIKOR method not applicable as all values of ranking statistic S ",
         "are identical.",
         call. = FALSE)
  #
  res$Q <-
    v * (res$S - min.S) / (max.S - min.S) +
    (1 - v) * (res$R - min.R) / (max.R - min.R)
  #
  res %<>% arrange(Q)
  #
  class(res) <- c("vikor", class(res))
  #
  attr(res, "weights") <- weights
  attr(res, "dropped_treatments") <- dropped_treatments
  #
  res
}
