forest.netmeta <- function(x,
                           reference.group=x$reference.group,
                           pooled="fixed",
                           leftcols="studlab",
                           leftlabs="Treatment",
                           ...){

  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")

  if (pooled=="fixed"){
    TE    <- x$TE.fixed
    seTE  <- x$seTE.fixed
  }
  
  if (all(colnames(TE)!=reference.group))
    stop(paste("Argument 'reference.group' must match any of the following values: ",
               paste(paste("'", colnames(TE), "'", sep=""),
                     collapse=" - "), sep=""))
  ##
  TE.b <- TE[,colnames(TE)==reference.group]
  seTE.b <- seTE[,colnames(seTE)==reference.group]

  m1 <- metagen(TE.b, seTE.b, sm=x$sm,
                studlab=colnames(TE))

  meta:::forest(m1,
                comb.fixed=FALSE, comb.random=FALSE,
                hetstat=FALSE,
                leftcols=leftcols,
                leftlabs=leftlabs,
                ...)

  invisible(NULL)
}
