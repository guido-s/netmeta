netrank <- function(x, small.values = "good") {
  
  ## Check for netmeta object
  ##
  meta:::chkclass(x, "netmeta")
  small.values <- meta:::setchar(small.values, c("good", "bad"))
  
  
  TE.fixed <- x$TE.fixed
  pval.fixed <- x$pval.fixed
  ##
  TE.random <- x$TE.random
  pval.random <- x$pval.random
  
  
  ## Calculate one-sided p-values
  ##
  w.fixed <- (1 + sign(TE.fixed)) / 2
  p.fixed <- pval.fixed
  ##
  if (small.values == "good")
    P.fixed <- w.fixed * p.fixed / 2 + (1 - w.fixed) * (1 - p.fixed / 2)
  else
    P.fixed <- w.fixed * (1 - p.fixed / 2) + (1 - w.fixed) * p.fixed / 2
  ##
  w.random <- (1 + sign(TE.random)) / 2
  p.random <- pval.random
  ##
  if (small.values == "good")
    P.random <- w.random * p.random / 2 + (1 - w.random) * (1 - p.random / 2)
  else
    P.random <- w.random * (1 - p.random / 2) + (1 - w.random) * p.random / 2
  
  
  ## Row means provide P-scores
  ##
  Pscore.fixed <- rowMeans(P.fixed, na.rm = TRUE)
  ##
  if (!all(is.na(TE.random)))
    Pscore.random <- rowMeans(P.random, na.rm = TRUE)
  else
    Pscore.random <- NA
  
  
  ##
  if (!x$comb.random) {
    Pscore <- Pscore.fixed
    Pmatrix <- P.fixed
  }
  else {
    Pscore <- Pscore.random
    Pmatrix <- P.random
  }
  
  
  res <- list(Pscore.fixed = Pscore.fixed,
              Pmatrix.fixed = P.fixed,
              Pscore.random = Pscore.random,
              Pmatrix.random = P.random,
              small.values = small.values,
              x = x,
              title = x$title,
              version = packageDescription("netmeta")$Version,
              Pscore = "'Pscore' replaced by 'Pscore.fixed' and 'Pscore.random'.",
              Pmatrix = "'Pmatrix' replaced by 'Pmatrix.fixed' and 'Pmatrix.random'.")

  class(res) <- "netrank"
  
  res
}
