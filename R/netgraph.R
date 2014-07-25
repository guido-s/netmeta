netgraph <- function(x, seq=x$seq,
                     labels=NULL, plastic, thickness,
                     col="slateblue", lwd=5,
                     highlight=NULL,
                     col.highlight="red2", lwd.highlight=lwd,
                     highlight.split=":",
                     ##
                     start.layout="circle",
                     iterate,
                     tol=0.0001, maxit=500, allfigures=FALSE,
                     eig1=2, eig2=3,
                     A.matrix=x$A.matrix,
                     N.matrix=sign(A.matrix),
                     ##
                     scale=1.10, offset=0.0175,
                     ##
                     cex.labels=1,
                     ##
                     lwd.min=lwd/2,
                     lwd.max=lwd*4,
                     ##
                     comparisons=NULL,
                     cex=1,
                     ##
                     multiarm=any(x$narms>2),
                     col.multiarm=NULL,
                     multiarm.colfunc=colorspace::sequential_hcl,
                     ##
                     points=FALSE,
                     col.point="red", cex.point=1, pch.point=20,
                     ##
                     arrows=FALSE,
                     
                     ##
                     xpos=NULL, ypos=NULL,
                     ...){
  
  
  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")
  
  
  isl <- charmatch(tolower(start.layout),
                   c("eigen", "prcomp", "circle", "random"), nomatch = NA)
  ##
  if(is.na(isl))
    stop("Argument 'start.layout' should be \"eigen\", \"prcomp\", \"circle\", or \"random\"")
  ##
  start.layout <- c("eigen", "prcomp", "circle", "random")[isl]
  
  
  if (missing(iterate))
    iterate <- ifelse(start.layout=="circle", FALSE, TRUE)
  
  
  if (missing(plastic))
    if (start.layout=="circle" & iterate==FALSE)
      plastic <- TRUE
    else
      plastic <- FALSE
  
  
  if (missing(thickness)){
    if (start.layout=="circle" & iterate==FALSE & plastic==TRUE){
      thick <- "se.fixed"
      thickness <- "se.fixed"
    }
    else{
      thick <- "equal"
      thickness <- "equal"
    }
  }
  else{
    if (!is.matrix(thickness)){
      if (length(thickness)==1 & is.character(thickness)){
        ithick <- charmatch(tolower(thickness),
                            c("equal", "number.of.studies", "se.fixed", "se.random", "w.fixed", "w.random"), nomatch = NA)
        ##
        if(is.na(ithick))
          stop("Argument 'thickness' should be \"equal\", \"number.of.studies\", \"se.fixed\", \"se.random\", \"w.fixed\", or \"w.random\"")
        ##
        thick <- c("equal", "number.of.studies", "se.fixed", "se.random", "w.fixed", "w.random")[ithick]
      }
      ##
      else if (length(thickness)==1 & is.logical(thickness)){
        if (thickness)
          thick <- "w.fixed"
        else
          thick <- "equal"
      }
    }
    else{
      if ((dim(thickness)[1] != dim(A.matrix)[1]) |
          (dim(thickness)[2] != dim(A.matrix)[2]))
        stop("Dimension of argument 'A.matrix' and 'thickness' are different.")
      else{
        W.matrix <- thickness
        thick <- "matrix"
      }
    }
  }
  
  
  treatments <- dimnames(x$A.matrix)[[1]]
  ##
  if (is.null(seq) | !(start.layout=="circle" & iterate==FALSE)){
    seq1 <- 1:length(treatments)
    if (!missing(seq) & (is.null(xpos) & is.null(ypos)))
      warning("Argument 'seq' only considered if start.layout=\"circle\" and iterate=FALSE.")
  }
  else{
    if (!is.null(seq) & length(seq) != length(treatments))
      stop("Length of argument 'seq' different from number of treatments.")
    if (is.numeric(seq) && any(is.na(seq)))
      stop("Missing values not allowed in argument 'seq'.")
    if (length(unique(seq)) != length(seq))
      stop("Values for argument 'seq' must all be disparate.")
    if (is.numeric(seq) && any(!(seq %in% (1:length(treatments)))))
      stop(paste("Argument 'seq' must be a permutation of the integers from 1 to ",
                 length(treatments), ".", sep=""))
    if (is.character(seq) && any(!(seq %in% treatments)))
      stop(paste("Argument 'seq' must be a permutation of the following values:\n  ",
                 paste(paste("'", treatments, "'", sep=""),
                       collapse=" - "), sep=""))
    ##
    else if (is.numeric(seq))
      seq1 <- seq
    else if (is.character(seq))
      seq1 <- match(seq, treatments)
  }
  ##
  A.matrix <- A.matrix[seq1, seq1]
  N.matrix <- N.matrix[seq1, seq1]
  
  
  A.sign <- sign(A.matrix)
  
  
  if (is.null(xpos) & is.null(ypos)){
    stressdata <- stress(x,
                         ##
                         start.layout=start.layout,
                         iterate=iterate,
                         eig1=eig1, eig2=eig2,
                         tol=tol, maxit=maxit,
                         ##
                         seq=seq,
                         ##
                         A.matrix=A.matrix,
                         N.matrix=N.matrix,
                         ##
                         thickness=thickness,
                         plastic=plastic,
                         ##
                         scale=scale, offset=offset,
                         ##
                         labels=labels,
                         cex.labels=cex.labels,
                         ##
                         lwd=lwd,
                         lwd.min=lwd.min,
                         lwd.max=lwd.max,
                         ##
                         comparisons=comparisons,
                         col=col,
                         cex=cex,
                         ##
                         col.multiarm=col.multiarm,
                         multiarm.colfunc=multiarm.colfunc,
                         ##
                         highlight=highlight,
                         col.highlight=col.highlight,
                         lwd.highlight=lwd.highlight,
                         highlight.split=highlight.split,
                         ##
                         points=points,
                         col.point=col.point, cex.point=cex.point,
                         pch.point=pch.point,
                         ##
                         arrows=arrows,
                         ##
                         allfigures=allfigures,
                         ...)
    ##
    xpos <- stressdata$x
    ypos <- stressdata$y
    ##
    if (is.null(labels))
      labels <- stressdata$labels
  }
  else
    if (is.null(labels))
      labels <- treatments
    else
      labels <- labels[seq1]
  
  
  if (allfigures)
    return(invisible(NULL))
  
  
  n <- dim(A.matrix)[1]
  d <- scale*max(abs(c(min(c(xpos, ypos), na.rm=TRUE),
                       max(c(xpos, ypos), na.rm=TRUE))))
  
  
  ## Generate dataset for plotting
  ##
  pd <- data.frame(xpos, ypos, labels)
  pd$adj1 <- NA
  pd$adj2 <- NA
  ##
  pd$adj1[pd$xpos>0] <- 0
  pd$adj1[pd$xpos<0] <- 1
  ##
  pd$adj2[pd$ypos>0] <- 0
  pd$adj2[pd$ypos<0] <- 1
  ##
  offset <- offset*2*d
  ##
  pd$xpos.labels <- pd$xpos - offset + 2*(pd$adj1==0)*offset
  pd$ypos.labels <- pd$ypos - offset + 2*(pd$adj2==0)*offset
  
  
  ## Plot graph
  ##
  oldpar <- par(xpd=TRUE, pty="s")
  on.exit(par(oldpar))
  ##
  range <- c(-d, d)
  ##
  plot(xpos, ypos,
       xlim=range, ylim=range,
       type="n", bty="n", xlab="", ylab="", axes=FALSE,
       ...)
  
  
  ## Add coloured regions for multi-arm studies
  ##
  if (multiarm){
    td1 <- data.frame(studies=x$studies, narms=x$narms)
    td1 <- td1[rev(order(td1$narms)),]
    td1 <- td1[td1$narms>2,]
    multiarm.studies <- td1$studies
    ##
    n.multi <- length(multiarm.studies)
    ##
    if (is.null(col.multiarm)){
      ## Check for R package colorspace & use various gray values if
      ## not installed
      if (missing(multiarm.colfunc) &
          !any(as.data.frame(installed.packages())$Package=="colorspace")){
        col.multiarm <- paste("gray", round(seq.int(95, 30, length.out=n.multi)), sep="")
      }
      ##
      else{
        mcname <- deparse(substitute(multiarm.colfunc))
        ##
        csfun <- function(fcall, fname){
          is.cs <- length(grep(fname, fcall))>0
          if (is.cs)
            meta:::is.installed.package("colorspace")
          is.cs
        }
        ##
        if (csfun(mcname, "rainbow_hcl"))
          col.multiarm <- colorspace::rainbow_hcl(n.multi, start=240, end=60, alpha=0.5)
        else if (csfun(mcname, "sequential_hcl"))
          col.multiarm <- colorspace::sequential_hcl(n.multi, alpha=0.5)
        else if (csfun(mcname, "diverge_hcl"))
          col.multiarm <- colorspace::diverge_hcl(n.multi, alpha=0.5)
        else if (csfun(mcname, "heat_hcl"))
          col.multiarm <- colorspace::heat_hcl(n.multi, alpha=0.5)
        else if (csfun(mcname, "terrain_hcl"))
          col.multiarm <- colorspace::terrain_hcl(n.multi, alpha=0.5)
        else if (csfun(mcname, "diverge_hsv"))
          col.multiarm <- colorspace::diverge_hsv(n.multi, alpha=0.5)
        else if (csfun(mcname, "choose_palette")){
          fcolm <- colorspace::choose_palette(n=n.multi)
          col.multiarm <- fcolm(n=n.multi)
        }
        else 
          col.multiarm <- sapply(n.multi, multiarm.colfunc)
        ##
        if (csfun(mcname, "sequential_hcl")|
            csfun(mcname, "diverge_hcl")|
            csfun(mcname, "heat_hcl"))
          col.multiarm <- rev(col.multiarm)
      }
    }
    else
      if (length(col.multiarm) > 1 & length(col.multiarm) != n.multi)
        stop("Length of argument 'col.multiarm' must be equal to one or the number of multi-arm studies: ", n.multi)
    ##
    if (n.multi>0){
      multiarm.labels <- vector("list", n.multi)
      if (length(col.multiarm)==1)
        col.multiarm <- rep(col.multiarm, n.multi)
      for (i in 1:n.multi){
        treat1 <- x$treat1[x$studlab %in% multiarm.studies[i]]
        treat2 <- x$treat2[x$studlab %in% multiarm.studies[i]]
        multiarm.labels[[i]] <- sort(unique(c(treat2, treat1)))
        ##
        pdm <- pd[pd$labels %in% multiarm.labels[[i]],]
        ##
        ## Clockwise ordering of polygon coordinates
        ##
        polysort <- function(x, y) {
          xnorm <- (x-mean(x))/sd(x) # Normalise coordinate x
          ynorm <- (y-mean(y))/sd(y) # Normalise coordinate y
          r <- sqrt(xnorm^2 + ynorm^2) # Calculate polar coordinates
          cosphi <- xnorm/r
          sinphi <- ynorm/r
          s <- as.numeric(sinphi>0) # Define angles to lie in [0, 2*pi]
          phi <- acos(cosphi)
          alpha <- s*phi + (1-s)*(2*pi-phi)
          ##
          res <- order(alpha)
          res
        }
        ##
        pdm <- pdm[polysort(pdm$xpos, pdm$ypos),]
        ##
        polygon(pdm$xpos, pdm$ypos,
                col=col.multiarm[i], border=NA)
      }
    }
  }
  
  
  if (thick=="number.of.studies"){
    W.matrix <- lwd.max*A.matrix/max(A.matrix)
    W.matrix[W.matrix<lwd.min&W.matrix!=0] <- lwd.min
    ##
  }
  else if (thick=="equal"){
    W.matrix <- lwd*A.sign
  }
  else if (thick=="se.fixed"){
    IV.matrix <- x$seTE.direct.fixed
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max*min(IV.matrix, na.rm=TRUE)/IV.matrix
    W.matrix[W.matrix<lwd.min&W.matrix!=0] <- lwd.min
  }
  else if (thick=="se.random"){
    IV.matrix <- x$seTE.direct.random
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max*min(IV.matrix, na.rm=TRUE)/IV.matrix
    W.matrix[W.matrix<lwd.min&W.matrix!=0] <- lwd.min
  }
  else if (thick=="w.fixed"){
    IV.matrix <- 1/x$seTE.direct.fixed^2
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max*IV.matrix/max(IV.matrix, na.rm=TRUE)
    W.matrix[W.matrix<lwd.min&W.matrix!=0] <- lwd.min
  }
  else if (thick=="w.random"){
    IV.matrix <- 1/x$seTE.direct.random^2
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max*IV.matrix/max(IV.matrix, na.rm=TRUE)
    W.matrix[W.matrix<lwd.min&W.matrix!=0] <- lwd.min
  }
  else if (thick=="matrix"){
    W.matrix[is.infinite(W.matrix)] <- NA
    if (min(W.matrix[W.matrix!=0], na.rm=TRUE)==max(W.matrix[W.matrix!=0], na.rm=TRUE))
      W.matrix <- lwd*W.matrix
    else
      W.matrix <- lwd.max*W.matrix/max(W.matrix, na.rm=TRUE)
    ##
    W.matrix[W.matrix<lwd.min&W.matrix!=0] <- lwd.min
  }
  
  
  if (is.null(comparisons) || !is.matrix(comparisons))
    comparisons <- matrix(NA, nrow=n, ncol=n)
  
  
  ## Draw lines
  ##
  if (plastic){
    n.plastic <- 30
    lwd.multiply <- rep(NA, n.plastic)
    cols <- rep("", n.plastic)
    j <- 0
    for (i in n.plastic:1){
      j <- j+1
      lwd.multiply[j] <- sin(pi*i/2/n.plastic)
      cols[j] <- paste("gray", round(100*(1-i/n.plastic)), sep="")
    }
  }
  else{
    lwd.multiply <- 1
    cols <- col
  }
  ##
  for (n.plines in 1:length(lwd.multiply)){
    for (i in 1:(n-1)){
      for (j in (i+1):n){
        if (A.sign[i,j]>0)
          lines(c(xpos[i], xpos[j]), c(ypos[i], ypos[j]),
                lwd=W.matrix[i,j]*lwd.multiply[n.plines],
                col=cols[n.plines])
      }
    }
  }
  
  
  ## Draw arrows / print treatment labels
  ##
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if (A.sign[i,j]>0){
        ##
        if (arrows)
          arrows(xpos[i], ypos[i], xpos[j], ypos[j],
                 lwd=2, col=col,
                 code=2, angle=20, length=0.4)
        ##
        text((xpos[i]+xpos[j])/2, (ypos[i]+ypos[j])/2,
             labels=comparisons[i,j],
             cex=cex)
      }
    }
  }
  
  
  ## Add highlighted comparisons
  ##
  if (!is.null(highlight)){
    for (high in highlight){
      highs <- unlist(strsplit(high, split=highlight.split))
      if (length(highs)!=2)
        stop("Wrong format for argument 'highlight' (see helpfile of plotgraph command).")
      ##
      if (sum(pd$labels %in% highs)!=2)
        stop(paste("Argument 'highlight' must contain two of the following values (separated by \":\"):\n  ",
                 paste(paste("'", pd$labels, "'", sep=""),
                       collapse=" - "), sep=""))
      ##
      pdh <- pd[pd$labels %in% highs,]
      ##
      lines(pdh$xpos, pdh$ypos,
            lwd=W.matrix[labels==highs[1], labels==highs[2]],
            col=col.highlight)
    }
  }
  
  
  ## Add points for labels
  ##
  if (points)
    points(xpos, ypos,
           pch=pch.point, cex=cex.point, col=col.point)               
  

  ## Print treatment labels
  ##
  if (!is.null(labels))
    for (i in 1:n)
      text(pd$xpos.labels[i], pd$ypos.labels[i],
           labels=pd$labels[i],
           cex=cex,
           adj=c(pd$adj1[i], pd$adj2[i]))
  
  
  invisible(NULL)
}
