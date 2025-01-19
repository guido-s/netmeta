chkmultiarm <- function(TE, seTE, treat1, treat2, studlab, correlated,
                        tol.multiarm = 0.001,
                        tol.multiarm.se = NULL,
                        details = FALSE, debug = FALSE) {
  
  
  #
  # Ordering dataset (if necessary)
  #
  o <- order(studlab, treat1, treat2)
  #
  if (any(o != seq_along(studlab))) {
    TE <- TE[o]
    seTE <- seTE[o]
    treat1 <- treat1[o]
    treat2 <- treat2[o]
    studlab <- studlab[o]
    correlated <- correlated[o]
  }
  
  
  tabnarms <- table(studlab)
  sel.multi <- tabnarms > 1
  #
  if (any(sel.multi)) {
    #
    msgdetails <-
      paste0("  - For more details, re-run netmeta() with argument ",
             "details.chkmultiarm = TRUE.\n")
    #
    studlab.multi <- names(tabnarms)[sel.multi]
    #
    # Check duplicate and incomplete comparisons
    #
    dat.duplicate <- dat.incomplete <-
      data.frame(studlab = "", treat1 = "", treat2 = "",
                 stringsAsFactors = FALSE)
    #
    incomplete <- rep_len(NA, sum(sel.multi))
    duplicate <- rep_len(NA, sum(sel.multi))
    #
    s.idx <- 0
    
    #
    # (1) Check for multi-arm studies with incomplete or duplicate treatment
    #     comparisons
    #
    for (s in studlab.multi) {
      s.idx <- s.idx + 1
      sel <- studlab == s
      #
      TE.s <- TE[sel]
      studlab.s <- studlab[sel]
      treat1.s <- treat1[sel]
      treat2.s <- treat2[sel]
      treats.s <- unique(c(treat1.s, treat2.s))
      #
      n <- (1 + sqrt(8 * length(TE.s) + 1)) / 2
      #
      incomplete[s.idx] <- length(treats.s) != n
      duplicate[s.idx] <- any(table(interaction(treat1.s, treat2.s)) > 1)
      #
      if (incomplete[s.idx])
        dat.incomplete <- rbind(dat.incomplete,
                               data.frame(studlab = studlab.s,
                                          treat1 = treat1.s,
                                          treat2 = treat2.s,
                                          stringsAsFactors = FALSE))
      #
      if (duplicate[s.idx])
        dat.duplicate <- rbind(dat.duplicate,
                               data.frame(studlab = studlab.s,
                                          treat1 = treat1.s,
                                          treat2 = treat2.s,
                                          stringsAsFactors = FALSE))
    }
    #
    if (details & any(incomplete)) {
      dat.incomplete <- dat.incomplete[-1, ]
      cat("\nMulti-arm studies with incomplete treatment comparisons:\n\n")
      prmatrix(dat.incomplete, quote = FALSE, right = TRUE,
               rowlab = rep("", dim(dat.incomplete)[1]))
      cat("\n")
    }
    #
    if (details & any(duplicate)) {
      dat.duplicate <- dat.duplicate[-1, ]
      cat("\nMulti-arm studies with duplicate treatment comparisons:\n\n")
      prmatrix(dat.duplicate, quote = FALSE, right = TRUE,
               rowlab = rep("", dim(dat.duplicate)[1]))
      cat("\n")
    }
    #
    if (any(incomplete)) {
      studlabs <- unique(dat.incomplete$studlab[-1])
      #
      if (length(studlabs) == 1)
        errmsg.incomplete <-
          paste0("  - Study '", studlabs,
                "' has an incomplete set of comparisons.\n")
      else
        errmsg.incomplete <-
          paste0("  - Studies with incomplete set of comparisons: ",
                paste(paste0("'", studlabs, "'"), collapse = ", "),
                "\n")
    }
    else
      errmsg.incomplete <- ""
    #
    if (any(duplicate)) {
      studlabs <- unique(dat.duplicate$studlab[-1])
      #
      if (length(studlabs) == 1)
        errmsg.duplicate <-
          paste0("  - Duplicate comparison in study '", studlabs, "'.\n")
      else
        errmsg.duplicate <-
          paste0("  - Studies with duplicate comparisons: ",
                paste(paste0("'", studlabs, "'"), collapse = ", "),
                "\n")
    }
    else
      errmsg.duplicate <- ""
    #
    if (any(incomplete) | any(duplicate))
      stop("Problem",
           if ((sum(incomplete) + sum(duplicate)) > 1) "s",
           " in multi-arm studies!\n",
           errmsg.incomplete, errmsg.duplicate,
           if (!details) msgdetails, call. = FALSE)
    
    #
    # (2) Check for (i) consistency of TE and varTE or (2) negative or zero
    #     treatment arm variance (zero variance only results in a warning)
    #
    dat.TE <- data.frame(studlab = "", treat1 = "", treat2 = "",
                         TE = NA, resid = NA,
                         stringsAsFactors = FALSE)
    #
    dat.varTE <- data.frame(studlab = "", treat1 = "", treat2 = "",
                            varTE = NA, resid.var = NA,
                            seTE = NA, resid.se = NA,
                            stringsAsFactors = FALSE)
    #
    dat.negative <- data.frame(studlab = "", treat = "", var.treat = NA,
                               stringsAsFactors = FALSE)
    #
    dat.zero <- data.frame(studlab = "", treat = "", var.treat = NA,
                           stringsAsFactors = FALSE)
    #
    inconsistent.TE <- inconsistent.varTE <-
      zero.sigma2 <- negative.sigma2 <- rep_len(NA, sum(sel.multi))
    #
    s.idx <- 0
    #
    for (s in studlab.multi) {
      s.idx <- s.idx + 1
      sel <- studlab == s
      #
      TE.s <- TE[sel]
      seTE.s <- seTE[sel]
      varTE.s <- seTE.s^2
      studlab.s <- studlab[sel]
      treat1.s <- treat1[sel]
      treat2.s <- treat2[sel]
      #
      correlated.s <- unique(correlated[sel])
      #
      n <- (1 + sqrt(8 * length(TE.s) + 1)) / 2
      #
      treats.s <- unique(c(treat1.s, treat2.s))
      studlab.s.arms <- rep_len(s, n)
      #
      # Create full edge-vertex incidence matrix
      #
      B <- createB(ncol = n)
      #
      # Check treatment estimates
      #
      TE.diff <- TE.s - B %*% as.vector(ginv(B) %*% TE.s)
      if (debug) {
        cat("*** TE.diff = TE - B %*% as.vector(ginv(B) %*% TE) ***\n")
        print(data.frame(TE = TE.s,
                         TE.calc =  B %*% as.vector(ginv(B) %*% TE.s),
                         TE.diff))
      }
      #
      inconsistent.TE[s.idx] <- any(abs(TE.diff) > tol.multiarm)
      #
      if (any(abs(TE.diff) > tol.multiarm))
        dat.TE <- rbind(dat.TE,
                        data.frame(studlab = studlab.s,
                                   treat1 = treat1.s,
                                   treat2 = treat2.s,
                                   TE = round(TE.s, 8),
                                   resid = round(TE.diff, 8),
                                   stringsAsFactors = FALSE))
      #
      # Check standard errors
      #
      if (correlated.s) {
        inconsistent.varTE[s.idx] <- 0
        negative.sigma2[s.idx] <- 0
        zero.sigma2[s.idx] <- 0
      }
      else {
        A <- abs(B)
        #
        sigma2 <- as.vector(ginv(A) %*% varTE.s)
        #
        varTE.diff <- varTE.s - A %*% sigma2
        #
        if (debug) {
          cat("*** varTE.diff = varTE - A %*% sigma2 ***\n")
          cat("*** with sigma2 = as.vector(ginv(A) %*% varTE) ***\n")
          print(data.frame(varTE = varTE.s,
                           varTE.calc = A %*% sigma2,
                           varTE.diff))
          print(data.frame(sigma2))
        }
        #
        if (!is.null(tol.multiarm.se))
          inconsistent.varTE[s.idx] <- any(abs(varTE.diff) > tol.multiarm.se^2)
        else
          inconsistent.varTE <- rep_len(FALSE, length(inconsistent.varTE))
        #
        is.negative <- sigma2 < 0
        negative.sigma2[s.idx] <- any(is.negative)
        zero.sigma2[s.idx] <- any(is_zero(sigma2[!is.negative]))
        #
        if (inconsistent.varTE[s.idx])
          dat.varTE <- rbind(dat.varTE,
                             data.frame(studlab = studlab.s,
                                        treat1 = treat1.s,
                                        treat2 = treat2.s,
                                        varTE = round(varTE.s, 8),
                                        resid.var = round(varTE.diff, 8),
                                        seTE = round(sqrt(varTE.s), 8),
                                        resid.se = sign(varTE.diff) *
                                          round(sqrt(abs(varTE.diff)), 8),
                                        stringsAsFactors = FALSE))
        #
        if (negative.sigma2[s.idx])
          dat.negative <- rbind(dat.negative,
                                data.frame(studlab = studlab.s.arms,
                                           treat = treats.s,
                                           var.treat = sigma2,
                                           stringsAsFactors = FALSE))
        #
        if (zero.sigma2[s.idx])
          dat.zero <- rbind(dat.zero,
                            data.frame(studlab = studlab.s.arms,
                                       treat = treats.s,
                                       var.treat = round(sigma2, 8),
                                       stringsAsFactors = FALSE))
      }
    }
    #
    iTE <- sum(inconsistent.TE)
    ivarTE <- sum(inconsistent.varTE)
    inconsistent <- iTE > 0 | ivarTE > 0
    #
    inegative.sigma2 <- sum(negative.sigma2)
    izero.sigma2 <- sum(zero.sigma2)
    zero <- izero.sigma2 > 0
    negative <- inegative.sigma2 > 0
    #
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
    #
    studlab.inconsistent <-
      unique(c(studlab.inconsistent.TE, studlab.inconsistent.varTE,
               studlab.zero.sigma2, studlab.negative.sigma2))
    
    #
    # Print information on deviations from consistency assumption in
    # multi-arm studies
    #
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
        dat.varTE$resid.var <- format(dat.varTE$resid.var)
        dat.varTE$seTE <- format(dat.varTE$seTE)
        dat.varTE$resid.se <- format(dat.varTE$resid.se)
        cat("\nMulti-arm studies with inconsistent",
            "variances / standard errors:\n\n")
        prmatrix(dat.varTE, quote = FALSE, right = TRUE,
                 rowlab = rep("", dim(dat.varTE)[1]))
        cat("\n")
      }
      #
      # Negative variance
      #
      if (length(dat.negative$studlab) > 1 & negative) {
        dat.negative <- dat.negative[-1, ]
        cat("\nMulti-arm studies with negative treatment arm variance:\n\n")
        prmatrix(dat.negative, quote = FALSE, right = TRUE,
                 rowlab = rep("", dim(dat.negative)[1]))
        cat("\n")
      }
      #
      # Zero variance
      #
      if (length(dat.zero$studlab) > 1 & zero) {
        dat.zero <- dat.zero[-1, ]
        cat("\nMulti-arm studies with zero treatment arm variance:\n\n")
        prmatrix(dat.zero, quote = FALSE, right = TRUE,
                 rowlab = rep("", dim(dat.zero)[1]))
        cat("\n")
        #
        warning(
          paste0("Note, a zero treatment arm variance has been calculated ",
                 "for the following multi-arm ",
                 if (izero.sigma2 == 1) "study" else "studies",
                 ": ",
                 paste(paste0("'", studlab.zero.sigma2, "'"), collapse = ", ")),
          call. = FALSE)
      }
      #
      cat("Legend:\n")
      if (inconsistent)
        cat(" resid - residual deviation (observed minus expected)\n")
      if (iTE > 0)
        cat(" TE", if (ivarTE > 0) "   ", " - treatment estimate\n", sep = "")
      if (ivarTE > 0) {
        cat(" varTE - variance of treatment estimate\n")
        cat(" seTE  - standard error of treatment estimate\n")
      }
      if (negative | zero)
        cat(" var.treat - treatment arm variance\n")
      cat("\n")
    }
    
    
    #
    # Generate error message
    #
    if (inconsistent | negative) {
      #
      if (iTE > 0)
        msgTE <-
          paste0("  ",
                 if (iTE == 1) "- Study " else "- Studies ",
                 "with inconsistent treatment estimates: ",
                 paste(paste0("'", studlab.inconsistent.TE, "'"),
                       collapse = ", "),
                 "\n")
      else
        msgTE <- ""
      #
      if (ivarTE > 0)
        msgvarTE <-
          paste0("  ",
                 if (ivarTE == 1) "- Study " else "- Studies ",
                 "with inconsistent variances: ",
                 paste(paste0("'", studlab.inconsistent.varTE, "'"),
                       collapse = ", "),
                 "\n")
      else
        msgvarTE <- ""
      #
      if (negative)
        msgsigma2 <-
          paste0("  ",
                 if (inegative.sigma2 == 1) "- Study " else "- Studies ",
                 "with negative treatment arm variance: ",
                 paste(paste0("'", studlab.negative.sigma2, "'"),
                       collapse = ", "),
                 "\n",
                 "    Potential solutions:\n",
                 "    1. Use argument 'func.inverse' to specify a different ",
                 "function for matrix inversion;\n",
                 "    2. Use argument 'correlated' to identify studies with ",
                 "correlated treatment arms, e.g., due to body-split design;\n",
                 "    3. Fix data errors.\n")
      else
        msgsigma2 <- ""
      #
      errmsg <-
        paste0("Problem",
               if ((iTE + ivarTE + inegative.sigma2 + izero.sigma2) > 1) "s",
               " in multi-arm studies!\n",
               msgTE, msgvarTE, msgsigma2,
               "  - Please check original data used as input to netmeta().\n",
               if (!details) msgdetails,
               if (inconsistent)
                 paste0("  - Argument",
                        if (iTE & ivarTE) "s",
                        if (iTE) " 'tol.multiarm'",
                        if (iTE & ivarTE) " and",
                        if (ivarTE) " 'tol.multiarm.se'",
                        " in netmeta() can be used to",
                        if (iTE & ivarTE) "\n   ",
                        " relax consistency",
                        if (!(iTE & ivarTE)) "\n   ",
                        " assumption for multi-arm studies (if appropriate)."))
      #
      stop(errmsg, call. = FALSE)
    }
  }
  
  invisible(NULL)
}
