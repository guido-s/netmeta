print.summary.netmeta <- function(x,
                                  digits=max(3, .Options$digits - 3),
                                  logscale=FALSE,
                                  all.treatments=x$all.treatments,
                                  reference.group=x$reference.group,
                                  comb.fixed=x$comb.fixed,
                                  comb.random=x$comb.random,
                                  header=TRUE,
                                  ...){
  
  
  if (!inherits(x, "summary.netmeta"))
    stop("Argument 'x' must be an object of class \"summary.netmeta\"")
  
  
  rmSpace <- function(x, end=FALSE, pat=" "){
    ##
    if ( !end ){
      while (any(substring(x, 1, 1) == pat, na.rm=TRUE)){
        sel <- substring(x, 1, 1) == pat
        x[sel] <- substring(x[sel], 2)
      }
    }
    else{
      last <- nchar(x)
      ##
      while ( any(substring(x, last, last) == pat, na.rm=TRUE) ){
        sel <- substring(x, last, last) == pat
          x[sel] <- substring(x[sel], 1, last[sel]-1)
        last <- nchar(x)
      }
    }
    ##
    x
  }

  p.ci <- function(lower, upper){
    lower <- rmSpace(lower)
    upper <- rmSpace(upper)
    ##
    ifelse (lower=="NA" & upper=="NA",
            "",
            paste(" [", format(lower, justify="right"),
                  "; ", format(upper, justify="right"), "]", sep=""))
  }

  format.TE <- function(TE, na=FALSE){
    TE <- rmSpace(TE)
    if (na) res <- format(TE)
    else res <- ifelse(is.na(TE), "", format(TE))
    res
  }
  
  m <- x$m
  k <- x$k
  sm <- x$sm
  
  sm.lab <- sm
  ##
  if (logscale & (sm == "RR" | sm == "OR" | sm == "HR"))
    sm.lab <- paste("log", sm, sep="")

  ci.lab <- paste(round(100*x$level.comb, 1), "%-CI", sep="")

  
  TE.fixed    <- x$fixed$TE
  seTE.fixed  <- x$fixed$seTE
  lowTE.fixed <- x$fixed$lower
  uppTE.fixed <- x$fixed$upper
  ##
  TE.random    <- x$random$TE
  seTE.random  <- x$random$seTE
  lowTE.random <- x$random$lower
  uppTE.random <- x$random$upper
  
  
  if (!logscale & (sm == "RR" | sm == "OR" | sm == "HR")){
    TE.fixed    <- exp(TE.fixed)
    lowTE.fixed <- exp(lowTE.fixed)
    uppTE.fixed <- exp(uppTE.fixed)
    ##
    TE.random <- exp(TE.random)
    lowTE.random <- exp(lowTE.random)
    uppTE.random <- exp(uppTE.random)
  }
  
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
  I2 <- x$I2
  
  
  if (header)
    matitle(x)
  
  if (reference.group!="" & missing(all.treatments))
    all.treatments <- FALSE
  
  if (comb.fixed|comb.random){
    cat(paste("Number of studies: k=", k, "\n", sep=""))
    cat(paste("Number of pairwise comparisons: m=", m, "\n", sep=""))
    
    if (comb.fixed){
      if (all.treatments | reference.group!="")
        cat("\nFixed effect model\n")
      if (all.treatments){
        cat("\nTreatment estimate (sm='", sm.lab, "'):\n", sep="")
        print(TE.fixed)
        cat("\nLower ", 100*x$fixed$level, "%-confidence limit:\n", sep="")
        print(lowTE.fixed)
        cat("\nUpper ", 100*x$fixed$level, "%-confidence limit:\n", sep="")
        print(uppTE.fixed)
      }
      if (reference.group!=""){
        if (all(colnames(TE.fixed)!=reference.group))
          stop(paste("Argument 'reference.group' must match any of the following values: ",
                     paste(paste("'", colnames(TE.fixed), "'", sep=""),
                           collapse=" - "), sep=""))
        ##
        TE.fixed.b <- TE.fixed[,colnames(TE.fixed)==reference.group]
        lowTE.fixed.b <- lowTE.fixed[,colnames(lowTE.fixed)==reference.group]
        uppTE.fixed.b <- uppTE.fixed[,colnames(uppTE.fixed)==reference.group]
        ##
        res <- cbind(format.TE(TE.fixed.b, na=TRUE),
                     p.ci(format(lowTE.fixed.b), format(uppTE.fixed.b)))
        dimnames(res) <-
          list(colnames(TE.fixed), c(sm.lab, ci.lab))
        
        cat("\nTreatment estimate (sm='", sm.lab,
            "', reference.group='", reference.group, "'):\n", sep="")
        
        prmatrix(res, quote=FALSE, right=TRUE)
      }
    }
    
    if (comb.random){
      cat("\nRandom effects model\n\nTreatment estimate (sm='", sm.lab, "'):\n", sep="")
      print(TE.random)
      cat("\nLower ", 100*x$random$level, "%-confidence limit:\n", sep="")
      print(lowTE.random)
      cat("\nUpper ", 100*x$random$level, "%-confidence limit:\n", sep="")
      print(uppTE.random)
    }
    zlab <- "z"
    
    
    if (!is.na(x$tau))
      cat(paste("\nQuantifying heterogeneity/inconsistency:\n",
                if (x$tau^2 < 0.0001)
                "tau^2 < 0.0001"
                else
                paste("tau^2 = ",
                      format(round(x$tau^2, 4), 4, nsmall=4, scientific=FALSE), sep="")
                ,
                paste("; I^2 = ", round(I2, 1), "%",
                      "",
                      ##ifelse(FALSE,
                      ##       p.ci(paste(round(100*lowI2, 1), "%", sep=""),
                      ##            paste(round(100*uppI2, 1), "%", sep="")),
                      ##       ""),
                      sep=""),
                "\n", sep=""))
    

    
    if (m > 1){
      
      Qdata <- cbind(round(x$Q, 2), x$df,
                     meta:::format.p(1-pchisq(x$Q, df=x$df)))
      
      dimnames(Qdata) <- list("", c("Q", "d.f.", "p.value"))
      ##
      cat("\nTest of heterogeneity/inconsistency:\n")
      prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
      ##
      ##
    }
  }
  
  invisible(NULL)
}
