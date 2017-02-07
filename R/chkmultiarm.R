chkmultiarm <- function(treat1, treat2, TE, seTE, studlab) {
  tabnarms <- table(studlab)
  sel.multi <- tabnarms > 1
  ##
  if (any(sel.multi)) {
    studlab.multi <- names(tabnarms)[sel.multi]
    ##
    inconsistent.TE <- inconsistent.varTE <- rep_len(NA, sum(sel.multi))
    ##
    s.idx <- 0
    ##
    for (s in studlab.multi) {
      s.idx <- s.idx + 1
      sel <- studlab == s
      ##
      TE.s <- TE[sel]
      seTE.s <- seTE[sel]
      varTE.s <- seTE.s^2
      ##
      m <- length(TE.s)
      n <- (1 + sqrt(8 * m + 1)) / 2
      B <- matrix(0, nrow = m, ncol = n)
      ##
      k <- 0
      for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
          k <- k + 1
          B[k, i] <-  1
          B[k, j] <- -1
        }
      }
      ##
      ## Check treatment estimates
      ##
      TE.diff <- TE.s - B %*% as.vector(ginv(B) %*% TE.s)
      ##
      inconsistent.TE[s.idx] <- any(abs(TE.diff) > .Machine$double.eps^0.5)
      ##
      ## Check standard errors
      ##
      A <- abs(B)
      sigma <- as.vector(ginv(A) %*% varTE.s)
      ##
      varTE.diff <- varTE.s - A %*% sigma
      ##
      inconsistent.varTE[s.idx] <- any(abs(varTE.diff) > .Machine$double.eps^0.5)
    }
    ##
    if (sum(inconsistent.TE) > 0 | sum(inconsistent.varTE) > 0) {
      iTE <- sum(inconsistent.TE)
      ivarTE <- sum(inconsistent.varTE)
      ##
      if (iTE > 0)
        msgTE <- paste("  ",
                       if (iTE == 1) "Study " else "Studies ",
                       "with inconsistent treatment estimates: ",
                       paste(paste("'", studlab.multi[inconsistent.TE], "'", sep = ""),
                             collapse = ", "),
                       "\n",
                       sep = "")
      else
        msgTE <- ""
      ##
      if (ivarTE > 0)
        msgvarTE <- paste("  ",
                          if (ivarTE == 1) "Study " else "Studies ",
                          "with inconsistent variances: ",
                          paste(paste("'", studlab.multi[inconsistent.varTE], "'", sep = ""),
                                collapse = ", "),
                          "\n",
                          sep = "")
      else
        msgvarTE <- ""

      ##
      errmsg <- paste("Inconsistent ",
                      if (iTE > 0) "treatment effects ",
                      if (iTE > 0 & ivarTE > 0) "and ",
                      if (ivarTE > 0) "variances ",
                      "in multi-arm studies!\n",
                      msgTE, msgvarTE,
                      "  Please check original data used as input to netmeta().",
                      sep ="")
      ##
      stop(errmsg, call. = FALSE)
    }
  }
  
  
  invisible(NULL)
}
