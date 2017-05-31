chkmultiarm <- function(treat1, treat2, TE, seTE, studlab,
                        tol = .Machine$double.eps^0.5,
                        details = FALSE) {
  tabnarms <- table(studlab)
  sel.multi <- tabnarms > 1
  ##
  if (any(sel.multi)) {
    studlab.multi <- names(tabnarms)[sel.multi]
    ##
    dat.TE <- data.frame(studlab = "", treat1 = "", treat2 = "",
                         TE = NA, resid = NA,
                         stringsAsFactors = FALSE)
    ##
    dat.varTE <- data.frame(studlab = "", treat1 = "", treat2 = "",
                            varTE = NA, resid = NA,
                            stringsAsFactors = FALSE)
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
      studlab.s <- studlab[sel]
      treat1.s <- treat1[sel]
      treat2.s <- treat2[sel]
      ##
      m <- length(TE.s)
      n <- (1 + sqrt(8 * m + 1)) / 2
      ##
      ## Create full edge-vertex incidence matrix
      ##
      B <- createB(ncol = n)
      ##
      ## Check treatment estimates
      ##
      TE.diff <- TE.s - B %*% as.vector(ginv(B) %*% TE.s)
      ##
      inconsistent.TE[s.idx] <- any(abs(TE.diff) > tol)
      ##
      if (any(abs(TE.diff) > tol))
        dat.TE <- rbind(dat.TE,
                        data.frame(studlab = studlab.s,
                                   treat1 = treat1.s,
                                   treat2 = treat2.s,
                                   TE = round(TE.s, 8),
                                   resid = round(TE.diff, 8),
                                   stringsAsFactors = FALSE))
      ##
      ## Check standard errors
      ##
      A <- abs(B)
      sigma <- as.vector(ginv(A) %*% varTE.s)
      ##
      varTE.diff <- varTE.s - A %*% sigma
      ##
      inconsistent.varTE[s.idx] <- any(abs(varTE.diff) > tol)
      ##
      if (any(abs(varTE.diff) > tol))
        dat.varTE <- rbind(dat.varTE,
                           data.frame(studlab = studlab.s,
                                      treat1 = treat1.s,
                                      treat2 = treat2.s,
                                      varTE = round(varTE.s, 8),
                                      resid = round(varTE.diff, 8),
                                      stringsAsFactors = FALSE))
    }
    ##
    iTE <- sum(inconsistent.TE)
    ivarTE <- sum(inconsistent.varTE)
    ##
    studlab.inconsistent.TE <- character()
    studlab.inconsistent.varTE <- character()
    if (iTE > 0 )
      studlab.inconsistent.TE <- studlab.multi[inconsistent.TE]
    if (ivarTE > 0 )
      studlab.inconsistent.varTE <- studlab.multi[inconsistent.varTE]
    ##
    studlab.inconsistent <- unique(c(studlab.inconsistent.TE, studlab.inconsistent.varTE))
    
    
    ##
    ## Print information on deviations from consistency assumption in
    ## multi-arm studies
    ##
    if (details & (sum(inconsistent.TE) > 0 | sum(inconsistent.varTE) > 0)) {
      if (length(dat.TE$studlab) > 1) {
        dat.TE <- dat.TE[-1, ]
        cat("\nMulti-arm studies with inconsistent treatment effects:\n\n")
        prmatrix(dat.TE, quote = FALSE, right = TRUE,
                 rowlab = rep("", dim(dat.TE)[1]))
      }
      if (length(dat.varTE$studlab) > 1) {
        dat.varTE <- dat.varTE[-1, ]
        cat("\nMulti-arm studies with inconsistent variances:\n\n")
        prmatrix(dat.varTE, quote = FALSE, right = TRUE,
                 rowlab = rep("", dim(dat.TE)[1]))
      }
      ##
      cat("Legend:\n")
      cat(" resid - residual deviation (observed minus expected)\n")
      if (sum(inconsistent.TE) > 0)
        cat(" TE    - treatment estimate\n")
      if (sum(inconsistent.varTE) > 0)
        cat(" varTE - variance of treatment estimate\n")
      cat("\n")
    }
    
    
    ##
    ## Generate error message
    ##
    if (sum(inconsistent.TE) > 0 | sum(inconsistent.varTE) > 0) {
      ##
      if (iTE > 0)
        msgTE <- paste("  ",
                       if (iTE == 1) "- Study " else "- Studies ",
                       "with inconsistent treatment estimates: ",
                       paste(paste("'", studlab.inconsistent.TE, "'", sep = ""),
                             collapse = ", "),
                       "\n",
                       sep = "")
      else
        msgTE <- ""
      ##
      if (ivarTE > 0)
        msgvarTE <- paste("  ",
                          if (ivarTE == 1) "- Study " else "- Studies ",
                          "with inconsistent variances: ",
                          if (iTE > 0) "        ",
                          if ((iTE > 1 & ivarTE > 1) |
                              (iTE == 1 & ivarTE == 1)) "  ",
                          if (iTE > 1 & ivarTE == 1) "  ",
                          paste(paste("'", studlab.inconsistent.varTE, "'", sep = ""),
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
                      "in multi-arm ",
                      if (length(studlab.inconsistent) > 1) "studies" else "study",
                      "!\n",
                      msgTE,
                      msgvarTE,
                      "  - Please check original data used as input to netmeta().\n",
                      if (!details) "  - You may re-run netmeta() command with argument\n",
                      if (!details) "    details.tol.multiarm=TRUE to inspect deviations.\n",
                      "  - Argument tol.multiarm in netmeta() can be used to relax\n",
                      "    consistency assumption for multi-arm studies (if appropriate).",
                      sep ="")
      ##
      stop(errmsg, call. = FALSE)
    }
  }
  
  
  invisible(NULL)
}
