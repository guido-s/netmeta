netgraph <- function(x, seq=NULL,
                     labels=NULL, plastic=TRUE, thickness=TRUE,
                     col="slateblue", lwd=1,
                     highlight=NULL,
                     col.highlight="red2", lwd.highlight=2,
                     ...){
  
  
  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")
  
  
  nmak <- nma.krahn(x)
  ##
  comparisons <- nmak$direct$comparison
  ##
  X <- nmak$X.full[comparisons,]
  S <- nmak$direct$seTE
  
  
  if (!is.null(seq) & length(seq) != length(nmak$trts))
    stop("Length of argument 'seq' different from number of treatments.")
  if (is.numeric(seq) && any(is.na(seq)))
    stop("Missing values not allowed in argument 'seq'.")
  if (length(unique(seq)) != length(seq))
    stop("Values for argument 'seq' must all be disparate.")
  if (is.numeric(seq) && any(!(seq %in% (1:length(nmak$trts)))))
    stop(paste("Argument 'seq' must be a permutation of the integers from 1 to ",
               length(nmak$trts), ".", sep=""))
  if (is.character(seq) && any(!(seq %in% nmak$trts)))
    stop(paste("Argument 'seq' must be a permutation of the following values:\n  ",
               paste(paste("'", nmak$trts, "'", sep=""),
                     collapse=" - "), sep=""))
  
  
  if (is.null(labels))
    labels <- nmak$trts
  ##
  g <- ncol(X)+1
  if (is.null(seq))
    seq1 <- 1:g
  ##
  if (is.numeric(seq)){
    seq1 <- seq
    labels <- labels[seq]
  }
  if (is.character(seq)){
    seq1 <- match(seq, nmak$trts)
    labels <- labels[seq1]
  }
  ##
  X1 <- cbind(-rowSums(X), X)
  X1 <- X1[, seq1]
  
  
  xc <- cos((1:g)*2*pi/g)
  yc <- sin((1:g)*2*pi/g)
  
  
  plot(1.4*c(-1,1), 1.4*c(-1,1),
       asp=1, type="n", axes=FALSE, xlab="", ylab="")
  ##
  lines(xc, yc, type="p")
  
  
  if (plastic){
    lwds <- (min(S)/S)^lwd*20
    for (ii in (n <- 30):1){
      for (i in 1:nrow(X)){
        sel <- X1[i,]!=0
        lines(xc[sel], yc[sel],
              lwd=lwds[i]*sin(pi*ii/2/n),
              col=paste("gray",round(100*(1-ii/n)),sep=""))
      }
    }
  }
  ##
  if (plastic==FALSE & thickness){
    lwds <- (min(S)/S)^lwd*20
    for (ii in (n <- 30):1){
      for (i in 1:nrow(X)){
        sel <- X1[i,]!=0
        lines(xc[sel], yc[sel],
              lwd=lwds[i]*sin(pi*ii/2/n),
              col=col)
      }
    }
  }
  ##
  if (plastic==FALSE & thickness==FALSE){
    for (i in 1:nrow(X)){
      sel <- X1[i,]!=0
      lines(xc[sel], yc[sel],
            lwd=lwd*10, col=col)
    }
  }
  
  
  if (length(highlight)>0){
    for (i in 1:length(highlight)){
      if (highlight[i] %in% comparisons){
        lines(xc[X1[highlight[i],]!=0],
              yc[X1[highlight[i],]!=0],
              col=col.highlight, lwd=lwd.highlight)
      }
      else{
        h <- unlist(lapply(strsplit(highlight[i], ":"),
                           function(x) paste(x[2], x[1],sep=":")
                           )
                    )
        lines(xc[X1[h,]!=0], yc[X1[h,]!=0],
              col=col.highlight, lwd=lwd.highlight)
      }
    }
  }
  
  
  text(xc*1.25, yc*1.25, labels, cex=1.3)
  
  
  invisible(NULL)
}
