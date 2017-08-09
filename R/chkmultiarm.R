chkmultiarm <- function(treat1, treat2, TE, seTE, studlab,
                        tol = .Machine$double.eps^0.5,
                        details = FALSE) {
  tabnarms <- table(studlab)
  sel.multi <- tabnarms > 1
  ##
  if (any(sel.multi)) {
    ##
    msgdetails <- c("  - For more details, re-run netmeta() with argument details.chkmultiarm=TRUE.\n")
    ##
    studlab.multi <- names(tabnarms)[sel.multi]
    ##
    ## Check duplicate and incomplete comparisons
    ##
    dat.duplicate <- dat.incomplete <-
      data.frame(studlab = "", treat1 = "", treat2 = "",
                 stringsAsFactors = FALSE)
    ##
    incomplete <- rep_len(NA, sum(sel.multi))
    duplicate <- rep_len(NA, sum(sel.multi))
    ##
    s.idx <- 0
    ##
    for (s in studlab.multi) {
      s.idx <- s.idx + 1
      sel <- studlab == s
      ##
      TE.s <- TE[sel]
      studlab.s <- studlab[sel]
      treat1.s <- treat1[sel]
      treat2.s <- treat2[sel]
      treats.s <- unique(c(treat1.s, treat2.s))
      ##
      n <- (1 + sqrt(8 * length(TE.s) + 1)) / 2
      ##
      incomplete[s.idx] <- length(treats.s) != n
      duplicate[s.idx] <- any(table(interaction(treat1.s, treat2.s)) > 1)
      ##
      if (incomplete[s.idx])
        dat.incomplete <- rbind(dat.incomplete,
                               data.frame(studlab = studlab.s,
                                          treat1 = treat1.s,
                                          treat2 = treat2.s,
                                          stringsAsFactors = FALSE))
      ##
      if (duplicate[s.idx])
        dat.duplicate <- rbind(dat.duplicate,
                               data.frame(studlab = studlab.s,
                                          treat1 = treat1.s,
                                          treat2 = treat2.s,
                                          stringsAsFactors = FALSE))
    }
    ##
    if (details & any(incomplete)) {
      dat.incomplete <- dat.incomplete[-1, ]
      cat("\nMulti-arm studies with incomplete treatment comparisons:\n\n")
      prmatrix(dat.incomplete, quote = FALSE, right = TRUE,
               rowlab = rep("", dim(dat.incomplete)[1]))
      cat("\n")
    }
    ##
    if (details & any(duplicate)) {
      dat.duplicate <- dat.duplicate[-1, ]
      cat("\nMulti-arm studies with duplicate treatment comparisons:\n\n")
      prmatrix(dat.duplicate, quote = FALSE, right = TRUE,
               rowlab = rep("", dim(dat.duplicate)[1]))
      cat("\n")
    }
    ##
    if (any(incomplete)) {
      studlabs <- unique(dat.incomplete$studlab[-1])
      ##
      if (length(studlabs) == 1)
        errmsg.incomplete <-
          paste("  - Study '", studlabs,
                "' has an incomplete set of comparisons.\n",
                sep = "")
      else
        errmsg.incomplete <-
          paste("  - Studies with incomplete set of comparisons: ",
                paste(paste("'", studlabs, "'", sep = ""),
                      collapse = ", "),
                "\n", sep = "")
    }
    else
      errmsg.incomplete <- ""
    ##
    if (any(duplicate)) {
      studlabs <- unique(dat.duplicate$studlab[-1])
      ##
      if (length(studlabs) == 1)
        errmsg.duplicate <-
          paste("  - Duplicate comparison in study '", studlabs, "'.\n",
                sep = "")
      else
        errmsg.duplicate <-
          paste("  - Studies with duplicate comparisons: ",
                paste(paste("'", studlabs, "'", sep = ""),
                      collapse = ", "),
                "\n", sep = "")
    }
    else
      errmsg.duplicate <- ""
    ##
    if (any(incomplete) | any(duplicate))
      stop("Problem",
           if ((sum(incomplete) + sum(duplicate)) > 1) "s",
           " in multi-arm studies!\n",
           errmsg.incomplete, errmsg.duplicate,
           if (!details) msgdetails, call. = FALSE)
    ##
    ## Additional checks:
    ## - consistency of TE and varTE and
    ## - negative or zero treatment arm variance
    ##
    dat.TE <- data.frame(studlab = "", treat1 = "", treat2 = "",
                         TE = NA, resid = NA,
                         stringsAsFactors = FALSE)
    ##
    dat.varTE <- data.frame(studlab = "", treat1 = "", treat2 = "",
                            varTE = NA, resid = NA,
                            stringsAsFactors = FALSE)
    ##
    dat.zero <- data.frame(studlab = "", treat = "", var.treat = NA,
                           stringsAsFactors = FALSE)
    ##
    dat.negative <- data.frame(studlab = "", treat = "", var.treat = NA,
                               stringsAsFactors = FALSE)
    ##
    inconsistent.TE <- inconsistent.varTE <-
      zero.sigma2 <- negative.sigma2 <- rep_len(NA, sum(sel.multi))
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
      n <- (1 + sqrt(8 * length(TE.s) + 1)) / 2
      ##
      treats.s <- unique(c(treat1.s, treat2.s))
      studlab.s.arms <- rep_len(s, n)
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
      sigma2 <- as.vector(ginv(A) %*% varTE.s)
      ##
      varTE.diff <- varTE.s - A %*% sigma2
      ##
      inconsistent.varTE[s.idx] <- any(abs(varTE.diff) > tol)
      ##
      is.zero <- abs(sigma2) < .Machine$double.eps^0.75
      zero.sigma2[s.idx] <- any(is.zero)
      negative.sigma2[s.idx] <- any(sigma2[!is.zero] < 0)
      ##
      if (inconsistent.varTE[s.idx])
        dat.varTE <- rbind(dat.varTE,
                           data.frame(studlab = studlab.s,
                                      treat1 = treat1.s,
                                      treat2 = treat2.s,
                                      varTE = round(varTE.s, 8),
                                      resid = round(varTE.diff, 8),
                                      stringsAsFactors = FALSE))
      ##
      if (zero.sigma2[s.idx])
        dat.zero <- rbind(dat.zero,
                          data.frame(studlab = studlab.s.arms,
                                     treat = treats.s,
                                     var.treat = round(sigma2, 8),
                                     stringsAsFactors = FALSE))
      ##
      if (negative.sigma2[s.idx])
        dat.negative <- rbind(dat.negative,
                              data.frame(studlab = studlab.s.arms,
                                         treat = treats.s,
                                         var.treat = sigma2,
                                         stringsAsFactors = FALSE))
    }
    ##
    iTE <- sum(inconsistent.TE)
    ivarTE <- sum(inconsistent.varTE)
    inconsistent <- iTE > 0 | ivarTE > 0
    ##
    inegative.sigma2 <- sum(negative.sigma2)
    izero.sigma2 <- sum(zero.sigma2)
    zero <- izero.sigma2 > 0
    negative <- inegative.sigma2 > 0
    ##
    studlab.inconsistent.TE <- character()
    studlab.inconsistent.varTE <- character()
    studlab.zero.sigma2 <- character()
    studlab.negative.sigma2 <- character()
    if (iTE > 0)
      studlab.inconsistent.TE <- studlab.multi[inconsistent.TE]
    if (ivarTE > 0)
      studlab.inconsistent.varTE <- studlab.multi[inconsistent.varTE]
    if (zero)
      studlab.zero.sigma2 <- studlab.multi[zero.sigma2]
    if (negative)
      studlab.negative.sigma2 <- studlab.multi[negative.sigma2]
    ##
    studlab.inconsistent <- unique(c(studlab.inconsistent.TE, studlab.inconsistent.varTE,
                                     studlab.zero.sigma2, studlab.negative.sigma2))
    
    
    ##
    ## Print information on deviations from consistency assumption in
    ## multi-arm studies
    ##
    if (details & (inconsistent | zero | negative)) {
      if (length(dat.TE$studlab) > 1) {
        dat.TE <- dat.TE[-1, ]
        dat.TE$TE <- format(dat.TE$TE)
        dat.TE$resid <- format(dat.TE$resid)
        cat("\nMulti-arm studies with inconsistent treatment effects:\n\n")
        prmatrix(dat.TE, quote = FALSE, right = TRUE,
                 rowlab = rep("", dim(dat.TE)[1]))
        cat("\n")
      }
      if (length(dat.varTE$studlab) > 1 & ivarTE > 0) {
        dat.varTE <- dat.varTE[-1, ]
        dat.varTE$varTE <- format(dat.varTE$varTE)
        dat.varTE$resid <- format(dat.varTE$resid)
        cat("\nMulti-arm studies with inconsistent variances:\n\n")
        prmatrix(dat.varTE, quote = FALSE, right = TRUE,
                 rowlab = rep("", dim(dat.varTE)[1]))
        cat("\n")
      }
      ##
      ## Zero variance
      ##
      if (length(dat.zero$studlab) > 1 & zero) {
        dat.zero <- dat.zero[-1, ]
        cat("\nMulti-arm studies with zero treatment arm variance:\n\n")
        prmatrix(dat.zero, quote = FALSE, right = TRUE,
                 rowlab = rep("", dim(dat.zero)[1]))
        cat("\n")
      }
      ##
      ## Negative variance
      ##
      if (length(dat.negative$studlab) > 1 & negative) {
        dat.negative <- dat.negative[-1, ]
        cat("\nMulti-arm studies with negative treatment arm variance:\n\n")
        prmatrix(dat.negative, quote = FALSE, right = TRUE,
                 rowlab = rep("", dim(dat.negative)[1]))
        cat("\n")
      }
      ##
      cat("Legend:\n")
      if (inconsistent)
        cat(" resid - residual deviation (observed minus expected)\n")
      if (iTE > 0)
        cat(" TE    - treatment estimate\n")
      if (ivarTE > 0)
        cat(" varTE - variance of treatment estimate\n")
      if (negative | zero)
        cat(" var.treat - treatment arm variance\n")
      cat("\n")
    }
    
    
    ##
    ## Generate error message
    ##
    if (inconsistent | zero | negative) {
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
                          paste(paste("'", studlab.inconsistent.varTE, "'", sep = ""),
                                collapse = ", "),
                          "\n",
                          sep = "")
      else
        msgvarTE <- ""
      ##
      if (zero)
        msg0sigma2 <- paste("  ",
                           if (izero.sigma2 == 1) "- Study " else "- Studies ",
                           "with zero treatment arm variance: ",
                           paste(paste("'", studlab.zero.sigma2, "'", sep = ""),
                                 collapse = ", "),
                           "\n",
                           sep = "")
      else
        msg0sigma2 <- ""
      ##
      if (negative)
        msgsigma2 <- paste("  ",
                           if (inegative.sigma2 == 1) "- Study " else "- Studies ",
                           "with negative treatment arm variance: ",
                           paste(paste("'", studlab.negative.sigma2, "'", sep = ""),
                                 collapse = ", "),
                           "\n",
                           sep = "")
      else
        msgsigma2 <- ""
      ##
      errmsg <- paste("Problem",
                      if ((iTE + ivarTE + inegative.sigma2 + izero.sigma2) > 1) "s",
                      " in multi-arm studies!\n",
                      msgTE, msgvarTE, msg0sigma2, msgsigma2,
                      "  - Please check original data used as input to netmeta().\n",
                      if (!details) msgdetails,
                      if (inconsistent)
                        "  - Argument tol.multiarm in netmeta() can be used to relax consistency\n",
                      if (inconsistent)
                        "    assumption for multi-arm studies (if appropriate).",
                      sep ="")
      ##
      stop(errmsg, call. = FALSE)
    }
  }
  
  
  invisible(NULL)
}
