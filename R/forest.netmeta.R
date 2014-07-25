forest.netmeta <- function(x,
                           pooled=ifelse(x$comb.random, "random", "fixed"),
                           reference.group=x$reference.group,
                           leftcols="studlab",
                           leftlabs="Treatment",
                           smlab=NULL,
                           sortvar=x$seq,
                           ...){

  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")
  
  ipool <- charmatch(tolower(pooled), c("fixed", "random"), nomatch = NA)
  ##
  if (is.na(ipool)) 
        stop("Argument 'pooled' should be \"fixed\" or \"random\"")
  ##
  pooled <- c("fixed", "random")[ipool]
  
  if (pooled=="fixed"){
    TE   <- x$TE.fixed
    seTE <- x$seTE.fixed
    if (is.null(smlab))
      smlab <- "Fixed Effect Model"
  }
  ##
  if (pooled=="random"){
    TE   <- x$TE.random
    seTE <- x$seTE.random
    if (is.null(smlab))
      smlab <- "Random Effects Model"
  }

  if (!is.null(sortvar)){
    if (is.character(sortvar)){
      tlevs <- rownames(TE)
      ##
      if (length(tlevs)!=length(sortvar))
        stop("Length of argument 'sortvar' different from number of treatments.")
      else{
        if (length(unique(sortvar)) != length(sortvar))
          stop("Values for argument 'sortvar' must all be disparate.")
        if (any(!(sortvar %in% tlevs)))
          stop(paste("Argument 'sortvar' must be a permutation of the following values:\n  ",
                     paste(paste("'", tlevs, "'", sep=""),
                           collapse=" - "), sep=""))
      }
      TE <- TE[sortvar, sortvar]
      seTE <- seTE[sortvar, sortvar]
    }
    else{
      o <- order(sortvar)
      TE <- TE[o,o]
      seTE <- seTE[o,o]
    }
  }
  
  
  if (all(colnames(TE)!=reference.group))
    stop(paste("Argument 'reference.group' must match any of the following values: ",
               paste(paste("'", colnames(TE), "'", sep=""),
                     collapse=" - "), sep=""))
  ##
  TE.b <- TE[,colnames(TE)==reference.group]
  seTE.b <- seTE[,colnames(seTE)==reference.group]

  m1 <- metagen(TE.b, seTE.b, sm=x$sm,
                studlab=colnames(TE), warn=FALSE)
  
  forest(m1,
         comb.fixed=FALSE, comb.random=FALSE,
         hetstat=FALSE,
         leftcols=leftcols,
         leftlabs=leftlabs,
         smlab=smlab,
         ...)
  
  invisible(NULL)
}
