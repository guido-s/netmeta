print.netleague <- function(x,
                            comb.fixed = x$comb.fixed,
                            comb.random = x$comb.random,
                            ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  meta:::chkclass(x, "netleague")
  ##
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  
  
  ##
  ##
  ## (2) Print league table for fixed effect model
  ##
  ##
  if (comb.fixed) {
    cat("League table (fixed effect model):\n")
    ##
    prmatrix(x$fixed, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(x$fixed)),
             collab = rep("", ncol(x$fixed)))
    if (comb.random)
      cat("\n")
  }
  
  
  ##
  ##
  ## (3) Print league table for random effects model
  ##
  ##
  if (comb.random) {
    cat("League table (random effects model):\n")
    ##
    prmatrix(x$rando, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(x$random)),
             collab = rep("", ncol(x$random)))
  }
  
  
  invisible(NULL)
}
