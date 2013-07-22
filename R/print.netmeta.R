print.netmeta <- function(x,
                          sortvar,
                          level=x$level, level.comb=x$level.comb,
                          comb.fixed=x$comb.fixed, comb.random=x$comb.random,
                          reference.group=x$reference.group,
                          all.treatments=x$all.treatments,
                          details=TRUE, ma=TRUE, logscale=FALSE,
                          digits=max(4, .Options$digits - 3),
                          ...
                          ){
  
  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")
  
  
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
  
  
  k.all <- length(x$TE)
  ##
  if (missing(sortvar)) sortvar <- 1:k.all
  ##
  if (length(sortvar) != k.all)
    stop("'x' and 'sortvar' have different length")
  ##
  ci.lab <- paste(round(100*level, 1), "%-CI", sep="")

  sm <- x$sm
  
  sm.lab <- sm
  ##
  if (logscale & (sm == "RR" | sm == "OR" | sm == "HR"))
    sm.lab <- paste("log", sm, sep="")

  
  matitle(x)
  
  
  if (details){

    cat(paste("Original data",
              ifelse(any(x$narms>2),
                     " (with adjusted standard errors for multi-arm studies)",
                     ""),
              ":\n\n", sep=""))
    
    res <- data.frame(treat1=x$treat1,
                      treat2=x$treat2,
                      TE=format.TE(round(x$TE, digits), na=TRUE),
                      seTE=format(round(x$seTE, digits)),
                      studlab=x$studlab)
    ##
    if (any(x$narms>2)){
      tdata1 <- data.frame(studlab=as.character(x$studies),
                           narms=x$narms)
      res$OrDeR <- 1:dim(res)[[1]]
      res <- merge(res, tdata1,
                   by="studlab", all.x=TRUE, all.y=FALSE,
                   sort=FALSE)
      res <- res[order(res$OrDeR),]
      res$multiarm <- ifelse(res$narms>2, "*", "")
      res$OrDeR <- NULL
      res$studlab <- NULL
      res <- as.matrix(res)
    }
    ##
    dimnames(res)[[1]] <- x$studlab
    
    prmatrix(res[order(sortvar),], quote=FALSE, right=TRUE)
    cat("\n")
    
    cat("Number of treatment arms (by study):\n")
    prmatrix(data.frame(narms=x$narms, row.names=x$studies),
             quote=FALSE, right=TRUE)
    cat("\n")
  }
  
  
  tsum <- summary(x, level=level, level.comb=level.comb, warn=FALSE)
  ##
  TE.f    <- tsum$comparison.nma.fixed$TE
  lowTE.f <- tsum$comparison.nma.fixed$lower
  uppTE.f <- tsum$comparison.nma.fixed$upper
  ##
  if (!logscale & (sm == "RR" | sm == "OR" | sm == "HR")){
    TE.f    <- exp(TE.f)
    lowTE.f <- exp(lowTE.f)
    uppTE.f <- exp(uppTE.f)
  }
  ##
  TE.f <- round(TE.f, digits)
  lowTE.f <- round(lowTE.f, digits)
  uppTE.f <- round(uppTE.f, digits)
  ##
  ##
  TE.r    <- tsum$comparison.nma.random$TE
  lowTE.r <- tsum$comparison.nma.random$lower
  uppTE.r <- tsum$comparison.nma.random$upper
  ##
  if (!logscale & (sm == "RR" | sm == "OR" | sm == "HR")){
    TE.r    <- exp(TE.r)
    lowTE.r <- exp(lowTE.r)
    uppTE.r <- exp(uppTE.r)
  }
  ##
  TE.r <- round(TE.r, digits)
  lowTE.r <- round(lowTE.r, digits)
  uppTE.r <- round(uppTE.r, digits)
  
  if (comb.fixed)
    if (sum(x$w.fixed)>0)
      w.fixed.p <- 100*round(x$w.fixed/sum(x$w.fixed, na.rm=TRUE), 4)
    else w.fixed.p <- x$w.fixed
  
  if (comb.random)
    if (sum(x$w.random)>0)
      w.random.p <- 100*round(x$w.random/sum(x$w.random, na.rm=TRUE), 4)
    else w.random.p <- x$w.random
    

  res.f <- cbind(x$treat1, x$treat2,
                 format.TE(TE.f, na=TRUE),
                 p.ci(format(lowTE.f), format(uppTE.f)),
                 if (comb.fixed) format(w.fixed.p),
                 if (comb.fixed) format(round(x$Q.fixed, 2)),
                 if (comb.fixed) format(round(x$leverage.fixed, 2)))
  dimnames(res.f) <-
    list(x$studlab, c("treat1", "treat2",
                      sm.lab, ci.lab,
                      if (comb.fixed) "%W",
                      if (comb.fixed) "Q",
                      if (comb.fixed) "leverage"))
  
  res.r <- cbind(x$treat1, x$treat2,
                 format.TE(TE.r, na=TRUE),
                 p.ci(format(lowTE.r), format(uppTE.r)),
                 if (comb.random) format(w.random.p))
  dimnames(res.r) <-
    list(x$studlab, c("treat1", "treat2",
                      sm.lab, ci.lab,
                      if (comb.random) "%W"))
  
  
  if (comb.fixed){
    cat("Data utilised in network meta-analysis (fixed effect model):\n\n")
    
    prmatrix(res.f[order(sortvar),], quote=FALSE, right=TRUE)
    
    cat("\n")
  }
  
  if (comb.random){
    cat("Data utilised in network meta-analysis (random effects model):\n\n")
    
    prmatrix(res.r[order(sortvar),], quote=FALSE, right=TRUE)
    
    cat("\n")
  }
  
  if (reference.group!="" & missing(all.treatments))
    all.treatments <- FALSE
  
  if (ma)
    print(tsum, digits=digits,
          comb.fixed=comb.fixed, comb.random=comb.random,
          logscale=logscale,
          all.treatments=all.treatments,
          reference.group=reference.group,
          header=FALSE)
  
  invisible(NULL)
}
