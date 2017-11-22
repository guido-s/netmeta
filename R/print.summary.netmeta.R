print.summary.netmeta <- function(x,
                                  comb.fixed = x$comb.fixed,
                                  comb.random = x$comb.random,
                                  prediction = x$prediction,
                                  reference.group = x$reference.group,
                                  baseline.reference = x$baseline.reference,
                                  all.treatments = x$all.treatments,
                                  backtransf = x$backtransf,
                                  nchar.trts = x$nchar.trts,
                                  header = TRUE,
                                  digits = max(3, .Options$digits - 3),
                                  ...) {
  
  
  meta:::chkclass(x, "summary.netmeta")
  
  
  if (is.null(x$df.Q))
    oldversion <- TRUE
  else
    oldversion <- FALSE
  ##
  if (is.null(x$predict$lower))
    prediction <- FALSE
  ##
  if (is.null(x$backtransf))
    backtransf <- TRUE
  ##
  if (is.null(x$nchar.trts))
    nchar.trts <- 666
  
  
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  meta:::chklogical(prediction)
  meta:::chklogical(baseline.reference)
  ##
  meta:::chklogical(backtransf)
  meta:::chknumeric(nchar.trts, min = 1, single = TRUE)
  
  
  ##
  ## Additional arguments
  ##
  fun <- "print.summary.netmeta"
  addargs <- names(list(...))
  ##
  meta:::warnarg("logscale", addargs, fun, otherarg = "backtransf")
  
  
  k <- x$k
  m <- x$m
  n <- x$n
  sm <- x$sm
  
  sm.lab <- sm
  ##
  if (!backtransf & meta:::is.relative.effect(sm))
    sm.lab <- paste("log", sm, sep = "")

  ci.lab <- paste(round(100 * x$fixed$level, 1), "%-CI", sep = "")

  
  TE.fixed    <- x$fixed$TE
  seTE.fixed  <- x$fixed$seTE
  lowTE.fixed <- x$fixed$lower
  uppTE.fixed <- x$fixed$upper
  ##
  TE.random    <- x$random$TE
  seTE.random  <- x$random$seTE
  lowTE.random <- x$random$lower
  uppTE.random <- x$random$upper
  ##
  lowTE.predict <- x$predict$lower
  uppTE.predict <- x$predict$upper
  ##
  if (!is.null(x$seq)) {
    TE.fixed <- TE.fixed[x$seq, x$seq]
    seTE.fixed <- seTE.fixed[x$seq, x$seq]
    lowTE.fixed <- lowTE.fixed[x$seq, x$seq]
    uppTE.fixed <- uppTE.fixed[x$seq, x$seq]
    ##
    TE.random <- TE.random[x$seq, x$seq]
    seTE.random <- seTE.random[x$seq, x$seq]
    lowTE.random <- lowTE.random[x$seq, x$seq]
    uppTE.random <- uppTE.random[x$seq, x$seq]
    ##
    lowTE.predict <- lowTE.predict[x$seq, x$seq]
    uppTE.predict <- uppTE.predict[x$seq, x$seq]
  }
  
  
  noeffect <- 0
  ##
  if (backtransf & meta:::is.relative.effect(sm)) {
    noeffect <- 1
    ##
    TE.fixed    <- exp(TE.fixed)
    lowTE.fixed <- exp(lowTE.fixed)
    uppTE.fixed <- exp(uppTE.fixed)
    ##
    TE.random    <- exp(TE.random)
    lowTE.random <- exp(lowTE.random)
    uppTE.random <- exp(uppTE.random)
    ##
    if (prediction) {
      lowTE.predict <- exp(lowTE.predict)
      uppTE.predict <- exp(uppTE.predict)
    }
  }
  ##
  TE.fixed    <- round(TE.fixed, digits)
  lowTE.fixed <- round(lowTE.fixed, digits)
  uppTE.fixed <- round(uppTE.fixed, digits)
  pTE.fixed   <- x$fixed$p
  zTE.fixed   <- round(x$fixed$z, digits)
  ##
  TE.random    <- round(TE.random, digits)
  lowTE.random <- round(lowTE.random, digits)
  uppTE.random <- round(uppTE.random, digits)
  pTE.random   <- x$random$p
  zTE.random   <- round(x$random$z, digits)
  ##
  if (prediction) {
    lowTE.predict <- round(lowTE.predict, digits)
    uppTE.predict <- round(uppTE.predict, digits)
  }
  ##
  I2 <- x$I2
  
  
  if (header)
    matitle(x)
  
  
  if (reference.group != "" & missing(all.treatments))
    all.treatments <- FALSE
  ##
  if (reference.group != "")
    reference.group <- setref(reference.group, rownames(TE.fixed))
  
  
  if (comb.fixed | comb.random) {
    cat(paste("Number of studies: k = ", k, "\n", sep = ""))
    cat(paste("Number of treatments: n = ", n, "\n", sep = ""))
    cat(paste("Number of pairwise comparisons: m = ", m, "\n", sep = ""))
    if (!oldversion)
      cat(paste("Number of designs: d = ", x$d, "\n", sep = ""))
    
    
    if (reference.group != "")
      if (baseline.reference)
        comptext <- paste("comparison: ",
                          if (x$n == 2)
                            paste("'",
                                  treats(rownames(TE.fixed),
                                         nchar.trts)[rownames(TE.fixed)
                                                     != reference.group],
                                  "'", sep = "")
                          else
                            "other treatments",
                          " vs '", reference.group, "'", sep = "")
      else
        comptext <- paste("comparison: '",
                          reference.group, "' vs ",
                          if (x$n == 2)
                            paste("'",
                                  treats(rownames(TE.fixed),
                                         nchar.trts)[rownames(TE.fixed)
                                                     != reference.group],
                                  "'", sep = "")
                          else
                            "other treatments", sep = "")
    
    
    if (comb.fixed) {
      if (all.treatments | reference.group != "")
        cat("\nFixed effect model\n")
      if (all.treatments) {
        cat("\nTreatment estimate (sm = '", sm.lab, "'):\n", sep = "")
        ##
        TEf <- meta:::format.NA(TE.fixed, digits = digits)
        rownames(TEf) <- treats(TEf, nchar.trts)
        colnames(TEf) <- treats(TEf, nchar.trts, FALSE)
        ##
        if (all(diag(TE.fixed) == noeffect))
          diag(TEf) <- "."
        ##
        prmatrix(TEf, quote = FALSE, right = TRUE)
        ##
        cat("\nLower ", 100 * x$fixed$level, "%-confidence limit:\n", sep = "")
        ##
        lowTEf <- meta:::format.NA(lowTE.fixed, digits = digits)
        rownames(lowTEf) <- treats(lowTEf, nchar.trts)
        colnames(lowTEf) <- treats(lowTEf, nchar.trts, FALSE)
        ##
        if (all(diag(lowTE.fixed) == noeffect))
          diag(lowTEf) <- "."
        ##
        prmatrix(lowTEf, quote = FALSE, right = TRUE)
        ##
        cat("\nUpper ", 100 * x$fixed$level, "%-confidence limit:\n", sep = "")
        ##
        uppTEf <- meta:::format.NA(uppTE.fixed, digits = digits)
        rownames(uppTEf) <- treats(uppTEf, nchar.trts)
        colnames(uppTEf) <- treats(uppTEf, nchar.trts, FALSE)
        ##
        if (all(diag(uppTE.fixed) == noeffect))
          diag(uppTEf) <- "."
        ##
        prmatrix(uppTEf, quote = FALSE, right = TRUE)
        ##
        ## Print prediction intervals
        ##
        if (!comb.random & prediction & x$df.Q >= 2) {
          cat("\nPrediction intervals\n")
          ##
          cat("\nLower ", 100 * x$predict$level, "%-prediction limit:\n", sep = "")
          ##
          lowTEp <- meta:::format.NA(lowTE.predict, digits = digits)
          rownames(lowTEp) <- treats(lowTEp, nchar.trts)
          colnames(lowTEp) <- treats(lowTEp, nchar.trts, FALSE)
          ##
          if (all(diag(lowTE.predict) == noeffect))
            diag(lowTEp) <- "."
          ##
          prmatrix(lowTEp, quote = FALSE, right = TRUE)
          ##
          cat("\nUpper ", 100 * x$predict$level, "%-prediction limit:\n", sep = "")
          ##
          uppTEp <- meta:::format.NA(uppTE.predict, digits = digits)
          rownames(uppTEp) <- treats(uppTEp, nchar.trts)
          colnames(uppTEp) <- treats(uppTEp, nchar.trts, FALSE)
          ##
          if (all(diag(uppTE.predict) == noeffect))
            diag(uppTEp) <- "."
          ##
          prmatrix(uppTEp, quote = FALSE, right = TRUE)
        }
      }
      if (reference.group != "") {
        if (all(colnames(TE.fixed) != reference.group))
          stop(paste("Argument 'reference.group' must match any of the following values: ",
                     paste(paste("'", colnames(TE.fixed), "'", sep = ""),
                           collapse = " - "), sep = ""))
        ##
        if (baseline.reference) {
          TE.fixed.b <- TE.fixed[, colnames(TE.fixed) == reference.group]
          lowTE.fixed.b <- lowTE.fixed[, colnames(lowTE.fixed) == reference.group]
          uppTE.fixed.b <- uppTE.fixed[, colnames(uppTE.fixed) == reference.group]
        }
        else {
          TE.fixed.b <- TE.fixed[rownames(TE.fixed) == reference.group, ]
          lowTE.fixed.b <- lowTE.fixed[rownames(lowTE.fixed) == reference.group, ]
          uppTE.fixed.b <- uppTE.fixed[rownames(uppTE.fixed) == reference.group, ]
        }
        ##
        ## Add prediction interval (or not)
        ##
        if (!comb.random & prediction & x$df.Q >= 2) {
          if (baseline.reference) {
            lowTE.predict.b <- lowTE.predict[, colnames(lowTE.predict) == reference.group]
            uppTE.predict.b <- uppTE.predict[, colnames(uppTE.predict) == reference.group]
          }
          else {
            lowTE.predict.b <- lowTE.predict[rownames(lowTE.predict) == reference.group, ]
            uppTE.predict.b <- uppTE.predict[rownames(uppTE.predict) == reference.group, ]
          }
          ##
          pi.lab <- paste(round(100 * x$predict$level, 1), "%-PI", sep = "")
          ##
          res <- cbind(format.TE(TE.fixed.b, na = TRUE),
                       p.ci(meta:::format.NA(lowTE.fixed.b, digits = digits),
                            meta:::format.NA(uppTE.fixed.b, digits = digits)),
                       p.ci(meta:::format.NA(lowTE.predict.b, digits = digits),
                            meta:::format.NA(uppTE.predict.b, digits = digits)))
          dimnames(res) <- list(colnames(TE.fixed), c(sm.lab, ci.lab, pi.lab))
        }
        else {
          res <- cbind(format.TE(TE.fixed.b, na = TRUE),
                       p.ci(meta:::format.NA(lowTE.fixed.b, digits = digits),
                            meta:::format.NA(uppTE.fixed.b, digits = digits)))
          dimnames(res) <- list(colnames(TE.fixed), c(sm.lab, ci.lab))
        }
        ##
        if (TE.fixed.b[rownames(res) == reference.group] == noeffect)
          res[rownames(res) == reference.group, ] <- "."
        
        cat("\nTreatment estimate (sm = '", sm.lab,
            "', ", comptext, "):\n", sep = "")
        
        prmatrix(res, quote = FALSE, right = TRUE)
      }
    }
    
    
    if (comb.random) {
      if (all.treatments | reference.group != "")
        cat("\nRandom effects model\n")
      if (all.treatments) {
        cat("\nTreatment estimate (sm = '", sm.lab, "'):\n", sep = "")
        ##
        TEr <- meta:::format.NA(TE.random, digits = digits)
        rownames(TEr) <- treats(TEr, nchar.trts)
        colnames(TEr) <- treats(TEr, nchar.trts, FALSE)
        ##
        if (all(diag(TE.random) == noeffect))
          diag(TEr) <- "."
        ##
        prmatrix(TEr, quote = FALSE, right = TRUE)
        ##
        cat("\nLower ", 100 * x$random$level, "%-confidence limit:\n", sep = "")
        ##
        lowTEr <- meta:::format.NA(lowTE.random, digits = digits)
        rownames(lowTEr) <- treats(lowTEr, nchar.trts)
        colnames(lowTEr) <- treats(lowTEr, nchar.trts, FALSE)
        ##
        if (all(diag(lowTE.random) == noeffect))
          diag(lowTEr) <- "."
        ##
        prmatrix(lowTEr, quote = FALSE, right = TRUE)
        ##
        cat("\nUpper ", 100 * x$random$level, "%-confidence limit:\n", sep = "")
        ##
        uppTEr <- meta:::format.NA(uppTE.random, digits = digits)
        rownames(uppTEr) <- treats(uppTEr, nchar.trts)
        colnames(uppTEr) <- treats(uppTEr, nchar.trts, FALSE)
        ##
        if (all(diag(uppTE.random) == noeffect))
          diag(uppTEr) <- "."
        ##
        prmatrix(uppTEr, quote = FALSE, right = TRUE)
        ##
        ## Print prediction intervals
        ##
        if (prediction & x$df.Q >= 2) {
          cat("\nPrediction intervals\n")
          ##
          cat("\nLower ", 100 * x$predict$level, "%-prediction limit:\n", sep = "")
          ##
          lowTEp <- meta:::format.NA(lowTE.predict, digits = digits)
          rownames(lowTEp) <- treats(lowTEp, nchar.trts)
          colnames(lowTEp) <- treats(lowTEp, nchar.trts, FALSE)
          ##
          if (all(diag(lowTE.predict) == noeffect))
            diag(lowTEp) <- "."
          ##
          prmatrix(lowTEp, quote = FALSE, right = TRUE)
          ##
          cat("\nUpper ", 100 * x$predict$level, "%-prediction limit:\n", sep = "")
          ##
          uppTEp <- meta:::format.NA(uppTE.predict, digits = digits)
          rownames(uppTEp) <- treats(uppTEp, nchar.trts)
          colnames(uppTEp) <- treats(uppTEp, nchar.trts, FALSE)
          ##
          if (all(diag(uppTE.predict) == noeffect))
            diag(uppTEp) <- "."
          ##
          prmatrix(uppTEp, quote = FALSE, right = TRUE)
        }
      }
      if (reference.group != "") {
        if (all(colnames(TE.random) != reference.group))
          stop(paste("Argument 'reference.group' must match any of the following values: ",
                     paste(paste("'", colnames(TE.random), "'", sep = ""),
                           collapse = " - "), sep = ""))
        ##
        if (baseline.reference) {
          TE.random.b <- TE.random[, colnames(TE.random) == reference.group]
          lowTE.random.b <- lowTE.random[, colnames(lowTE.random) == reference.group]
          uppTE.random.b <- uppTE.random[, colnames(uppTE.random) == reference.group]
        }
        else {
          TE.random.b <- TE.random[colnames(TE.random) == reference.group]
          lowTE.random.b <- lowTE.random[rownames(lowTE.random) == reference.group, ]
          uppTE.random.b <- uppTE.random[rownames(uppTE.random) == reference.group, ]
        }
        ##
        ## Add prediction interval (or not)
        ##
        if (prediction & x$df.Q >= 2) {
          if (baseline.reference) {
            lowTE.predict.b <- lowTE.predict[, colnames(lowTE.predict) == reference.group]
            uppTE.predict.b <- uppTE.predict[, colnames(uppTE.predict) == reference.group]
          }
          else {
            lowTE.predict.b <- lowTE.predict[rownames(lowTE.predict) == reference.group, ]
            uppTE.predict.b <- uppTE.predict[rownames(uppTE.predict) == reference.group, ]
          }
          ##
          pi.lab <- paste(round(100 * x$predict$level, 1), "%-PI", sep = "")
          ##
          res <- cbind(format.TE(TE.random.b, na = TRUE),
                       p.ci(meta:::format.NA(lowTE.random.b, digits = digits),
                            meta:::format.NA(uppTE.random.b, digits = digits)),
                       p.ci(meta:::format.NA(lowTE.predict.b, digits = digits),
                            meta:::format.NA(uppTE.predict.b, digits = digits)))
          dimnames(res) <- list(colnames(TE.fixed), c(sm.lab, ci.lab, pi.lab))
        }
        else {
          res <- cbind(format.TE(TE.random.b, na = TRUE),
                       p.ci(meta:::format.NA(lowTE.random.b, digits = digits),
                            meta:::format.NA(uppTE.random.b, digits = digits)))
          dimnames(res) <- list(colnames(TE.fixed), c(sm.lab, ci.lab))
        }
        ##
        if (TE.random.b[rownames(res) == reference.group] == noeffect)
          res[rownames(res) == reference.group, ] <- "."
        
        cat("\nTreatment estimate (sm = '", sm.lab,
            "', ", comptext, "):\n", sep = "")
        
        prmatrix(res, quote = FALSE, right = TRUE)
      }
    }
    ##
    zlab <- "z"
    
    
    if (!is.null(x$tau.preset))
      tau <- x$tau.preset
    else
      tau <- x$tau
    ##
    if (!is.na(tau))
      cat(paste("\nQuantifying heterogeneity / inconsistency:\n",
                meta:::format.p(tau^2,
                                lab = TRUE, labval = "tau^2",
                                digits = 4,
                                lab.NA = "NA"),
                paste("; I^2 = ", round(I2, 1), "%",
                      "",
                      ##ifelse(FALSE,
                      ##       p.ci(paste(round(100 * lowI2, 1), "%", sep = ""),
                      ##            paste(round(100 * uppI2, 1), "%", sep = "")),
                      ##       ""),
                      sep = ""),
                "\n", sep = ""))
    
    
    if (m > 1) {
      
      if (oldversion) {
        Qdata <- cbind(round(x$Q, 2), x$df,
                       ifelse(x$df == 0, "--",
                              meta:::format.p(x$pval.Q)))
        
        dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
        ##
        cat("\nTest of heterogeneity / inconsistency:\n")
        prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
      }
      else {
        if (x$d == 1 | is.na(x$Q.heterogeneity) | is.na(x$Q.inconsistency)) {
          Qdata <- cbind(round(x$Q, 2), x$df.Q,
                         ifelse(x$df.Q == 0, "--",
                                meta:::format.p(x$pval.Q)))
          
          dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
          ##
          cat("\nTest of heterogeneity / inconsistency:\n")
          prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
        }
        else {
          Qs <- c(x$Q, x$Q.heterogeneity, x$Q.inconsistency)
          df.Qs <- c(x$df.Q, x$df.Q.heterogeneity, x$df.Q.inconsistency)
          pval.Qs <- c(x$pval.Q, x$pval.Q.heterogeneity, x$pval.Q.inconsistency)
          pval.Qs <- ifelse(df.Qs == 0, "--",
                            meta:::format.p(pval.Qs))
          cat("\nTests of heterogeneity (within designs) and inconsistency (between designs):\n")
          Qdata <- data.frame(Q = round(Qs, 2),
                              df = df.Qs,
                              pval = pval.Qs)
          names(Qdata) <- c("Q", "d.f.", "p-value")
          rownames(Qdata) <- c("Total",
                               "Within designs",
                               "Between designs")
          print(Qdata)
        }
      }
    }
    
    
    if (!is.null(x$tau.preset)) {
      cat("\nDetails:")
      ##
      tau2 <- x$tau.preset^2
      tau2 <- meta:::format.p(tau2, lab = TRUE, labval = "tau^2",
                              digits = 4,
                              lab.NA = "NA")
      ##
      cat(paste("\n- Preset between-study variance: ",
                tau2, "\n", sep = ""))
    }
  }
  
  
  if (any(rownames(TE.fixed) != treats(TE.fixed, nchar.trts))) {
    abbr <- unique(treats(TE.fixed, nchar.trts))
    full <- unique(rownames(TE.fixed))
    ##
    tmat <- data.frame(abbr, full)
    names(tmat) <- c("Abbreviation", "Treatment name")
    tmat <- tmat[order(tmat$Abbreviation), ]
    ##
    cat("\nLegend:\n")
    prmatrix(tmat, quote = FALSE, right = TRUE,
             rowlab = rep("", length(abbr))) 
  }
  
  
  invisible(NULL)
}
