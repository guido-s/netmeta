netcomb <- function(x,
                    inactive = NULL,
                    sep.components = "+",
                    C.matrix, ...) {
  
  
  meta:::chkclass(x, "netmeta")
  ##
  meta:::chkchar(sep.components)
  
  
  n <- x$n # number of treatments / combinations
  m <- x$m # number of comparisons
  ##
  B.matrix <- x$B.matrix
  trts <- colnames(B.matrix) # list of treatments
  
  
  if (missing(C.matrix)) {
    ##
    ## Create C-matrix from netmeta object
    ##
    C.matrix <- as.matrix(createC(x, sep.components, inactive))
    C.matrix <- C.matrix[trts, , drop = FALSE]
  }
  else {
    ##
    ## Check if appropriate C.matrix is provided
    ##
    if (!(is.matrix(C.matrix) | is.data.frame(C.matrix)))
      stop("Argument 'C.matrix' must be a matrix or data frame.",
           call. = FALSE)
    ##
    if (nrow(C.matrix) != n)
      stop("Argument 'C.matrix' has wrong number of rows",
           " (must be equal to number of treatments).",
           call. = FALSE)
    ##
    ## Row names must be a permutation of treatments
    ##
    wrong.labels <- FALSE
    ##
    if (is.null(rownames(C.matrix)))
      wrong.labels <- TRUE
    else {
      if (length(unique(trts)) == length(unique(tolower(trts))))
        idx <- charmatch(tolower(rownames(C.matrix)),
                         tolower(trts), nomatch = NA)
      else
        idx <- charmatch(rownames(C.matrix), trts, nomatch = NA)
      ##
      if (any(is.na(idx)) || any(idx == 0))
        wrong.labels <- TRUE
    }
    ##
    if (wrong.labels)
      stop(paste("Row names of argument 'C.matrix' must be a ",
                 "permutation of treatment names:\n  ",
                 paste(paste("'", trts, "'", sep = ""),
                       collapse = " - "), sep = ""),
           call. = FALSE)
    ##
    C.matrix <- C.matrix[trts, , drop = FALSE]
    ##
    if (is.data.frame(C.matrix))
      C.matrix <- as.matrix(C.matrix)
  }
  ##
  c <- ncol(C.matrix) # number of components
  
  
  ##
  ## Design matrix based on treatment components
  ##
  X <- B.matrix %*% C.matrix
  colnames(X) <- colnames(C.matrix)
  rownames(X) <- x$studlab
  
  
  ##
  ## Fixed and random effects models
  ##
  df.Q.comp <- x$df.Q + x$n - c
  df.Q.diff <- x$n - c
  ##
  res.f <- nma.additive(x$TE, x$w.fixed, x$studlab,
                        x$treat1, x$treat2, x$level.comb,
                        X, C.matrix,
                        x$Q, df.Q.comp, df.Q.diff)
  ##
  res.r <- nma.additive(x$TE, x$w.random, x$studlab,
                        x$treat1, x$treat2, x$level.comb,
                        X, C.matrix,
                        x$Q, df.Q.comp, df.Q.diff)
  
  
  res <- list(k = x$k, n = n, m = m, c = c,
              ##
              comparisons.fixed = res.f$comparisons,
              components.fixed = res.f$components,
              combinations.fixed = res.f$combinations,
              ##
              comparisons.random = res.r$comparisons,
              components.random = res.r$components,
              combinations.random = res.r$combinations,
              ##
              sm = x$sm,
              level.comb = x$level.comb,
              comb.fixed = x$comb.fixed,
              comb.random = x$comb.random,
              ##
              Q = x$Q,
              df.Q = x$df.Q,
              pval.Q = x$pval.Q,
              ##
              Q.comp.fixed = res.f$Q.comp,
              Q.comp.random = res.r$Q.comp,
              df.Q.comp = df.Q.comp,
              pval.Q.comp.fixed = res.f$pval.Q.comp,
              pval.Q.comp.random = res.r$pval.Q.comp,
              ##
              Q.diff.fixed = res.f$Q.diff,
              Q.diff.random = res.r$Q.diff,
              df.Q.diff = df.Q.diff,
              pval.Q.diff.fixed = res.f$pval.Q.diff,
              pval.Q.diff.random = res.r$pval.Q.diff,
              ##
              C.matrix = C.matrix, B.matrix = B.matrix, X = X,
              ##
              backtransf = x$backtransf,
              nchar.trts = x$nchar.trts,
              ##
              title = x$title,
              ##
              x = x,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- "netcomb"
  ##
  res
}
