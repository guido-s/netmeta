netrank <- function(x, small.values = "good") {
  
  ## Check for netmeta object
  ##
  meta:::chkclass(x, "netmeta")
  small.values <- meta:::setchar(small.values, c("good", "bad"))
  
  
  TE.random <- x$TE.random
  pval.random <- x$pval.random
  ##
  if (!is.null(x$seq)) {
    TE.random <- TE.random[x$seq, x$seq]
    pval.random <- pval.random[x$seq, x$seq]
  }
  
  
  ## Calculate one-sided p-values
  ##
  w <- (1 + sign(TE.random)) / 2
  p <- pval.random
  ##
  if (small.values == "good")
    P <- w * p / 2       + (1 - w) * (1 - p / 2)
  else
    P <- w * (1 - p / 2) + (1 - w) * p / 2
  
  
  ## Row means provide P-scores
  ##
  Pscore <- rowMeans(P[, ], na.rm = TRUE)
  
  
  res <- list(Pscore = Pscore,
              Pmatrix = P,
              small.values = small.values,
              x = x,
              title = x$title,
              version = packageDescription("netmeta")$Version)

  class(res) <- "netrank"
  
  res
}
