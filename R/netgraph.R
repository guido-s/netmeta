netgraph <- function(x, seq = x$seq,
                     labels = x$trts,
                     cex = 1, adj = NULL,
                     offset = if (!is.null(adj) && all(unique(adj) == 0.5)) 0 else 0.0175,
                     scale = 1.10,
                     col = "slateblue", plastic, thickness,
                     lwd = 5, lwd.min = lwd / 2.5, lwd.max = lwd * 4,
                     dim = "2d",
                     ##
                     highlight = NULL, col.highlight = "red2",
                     lwd.highlight = lwd,
                     ##
                     multiarm = any(x$narms > 2),
                     col.multiarm = NULL,
                     alpha.transparency = 0.5,
                     ##
                     points = FALSE, col.points = "red",
                     cex.points = 1, pch.points = 20,
                     ##
                     number.of.studies = FALSE,
                     cex.number.of.studies = cex,
                     col.number.of.studies = "white",
                     bg.number.of.studies = "black",
                     pos.number.of.studies = 0.5,
                     ##
                     start.layout = ifelse(dim == "2d", "circle", "eigen"),
                     eig1 = 2, eig2 = 3, eig3 = 4,
                     iterate,
                     tol = 0.0001, maxit = 500, allfigures = FALSE,
                     A.matrix = x$A.matrix,
                     N.matrix = sign(A.matrix),
                     D.matrix = netdistance(N.matrix),
                     ##
                     xpos = NULL, ypos = NULL, zpos = NULL,
                     ...) {
  
  
  meta:::chkclass(x, "netmeta")
  ##
  n.edges <- sum(x$A.matrix[upper.tri(x$A.matrix)] > 0)
  n.trts <- length(x$trts)
  
  
  dim <- meta:::setchar(dim, c("2d", "3d"))
  is_2d <- dim == "2d"
  is_3d <- !is_2d
  ##
  if (is_3d & !meta:::is.installed.package("rgl", stop = FALSE)) {
    warning(paste("2-D plot generated as package 'rgl' is missing.",
                  "\n  ",
                  "Please install library 'rgl' in order to produce 3-D plots",
                  "\n  ",
                  "(R command: 'install.packages(\"rgl\")').",
                  if (length(grep("darwin", R.Version()$os)) == 1)
                    "\n  Note, macOS users have to install XQuartz, see https://www.xquartz.org/.",
                  sep = ""))
    dim <- "2d"
    is_2d <- TRUE
    is_3d <- FALSE
  }
  ##
  start.layout <- meta:::setchar(start.layout, c("eigen", "prcomp", "circle", "random"))
  ##
  if (!missing(seq) & is.null(seq))
    stop("Argument 'seq' must be not NULL.")
  ##
  if (!missing(labels) & is.null(labels))
    stop("Argument 'labels' must be not NULL.")
  ##
  n.col <- length(col)
  if (n.col == 1)
    col <- rep(col, n.edges)
  else if (n.col != n.edges)
    stop("Length of argument 'col' (",
         n.col, ") is different from the number of direct pairwise comparisons (",
         n.edges, ")")
  ##
  n.pos <- length(pos.number.of.studies)
  if (n.pos == 1)
    pos.number.of.studies <- rep(pos.number.of.studies, n.edges)
  else if (n.pos != n.edges)
    stop("Length of argument 'pos.number.of.studies' (",
         n.pos, ") is different from the number of direct pairwise comparisons (",
         n.edges, ")")
  ##
  if (length(cex.points) == 1)
    cex.points <- rep(cex.points, n.trts)
  else if (length(cex.points) != n.trts)
    stop("Length of argument 'cex.points' must be equal to the number of treatments.")
  ##
  if (length(col.points) == 1)
    col.points <- rep(col.points, n.trts)
  else if (length(col.points) != n.trts)
    stop("Length of argument 'col.points' must be equal to the number of treatments.")
  ##
  if (length(pch.points) == 1)
    pch.points <- rep(pch.points, n.trts)
  else if (length(pch.points) != n.trts)
    stop("Length of argument 'pch.points' must be equal to number of treatments.")
  
  
  if (missing(iterate))
    iterate <- ifelse(start.layout == "circle", FALSE, TRUE)
  
  
  addargs <- names(list(...))
  if ("highlight.split" %in% addargs)
    warning("Argument 'highlight.split' has been removed from R function netgraph.\n",
            "  This argument has been replaced by argument 'sep.trts' in R function netmeta.")
  ##
  highlight.split <- x$sep.trts
  ##
  if (is.null(highlight.split))
    highlight.split <- ":"
  
  
  if (missing(plastic))
    if (start.layout == "circle" & iterate == FALSE & is_2d)
      plastic <- TRUE
    else
      plastic <- FALSE
  
  
  if (missing(thickness)) {
    if (start.layout == "circle" & iterate == FALSE & plastic == TRUE) {
      thick <- "se.fixed"
      thickness <- "se.fixed"
    }
    else {
      thick <- "equal"
      thickness <- "equal"
    }
  }
  else {
    if (!is.matrix(thickness)) {
      if (length(thickness) == 1 & is.character(thickness))
        thick <- meta:::setchar(thickness,
                                c("equal", "number.of.studies",
                                  "se.fixed", "se.random", "w.fixed", "w.random"))
      ##
      else if (length(thickness) == 1 & is.logical(thickness)) {
        if (thickness)
          thick <- "se.fixed"
        else
          thick <- "equal"
      }
    }
    else {
      if ((dim(thickness)[1] != dim(A.matrix)[1]) |
          (dim(thickness)[2] != dim(A.matrix)[2]))
        stop("Dimension of argument 'A.matrix' and 'thickness' are different.")
      if (is.null(dimnames(thickness)))
        stop("Matrix 'thickness' must have row and column names identical to argument 'A.matrix'.")
      else {
        if (any(rownames(thickness) != rownames(A.matrix)))
          stop("Row names of matrix 'thickness' must be identical to argument 'A.matrix'.")
        if (any(colnames(thickness) != colnames(A.matrix)))
          stop("Column names of matrix 'thickness' must be identical to argument 'A.matrix'.")
      }
      ##
      W.matrix <- thickness
      thick <- "matrix"
    }
  }
  
  
  if (allfigures & is_3d) {
    warning("Argument 'allfigures' set to FALSE for 3-D network plot.")
    allfigures <- FALSE
  }
  
  
  if (is.null(seq) | !(start.layout == "circle" & iterate == FALSE)) {
    seq1 <- 1:length(labels)
    if (!missing(seq) & !is.null(seq) & (is.null(xpos) & is.null(ypos)))
      warning("Argument 'seq' only considered if start.layout=\"circle\" and iterate=FALSE.")
  }
  else {
    rn <- rownames(x$TE.fixed)
    seq1 <- charmatch(setseq(seq, rn), rn)
  }
  ##
  col.matrix <- matrix("", nrow = n.trts, ncol = n.trts)
  dimnames(col.matrix) <- dimnames(A.matrix)
  col.matrix <- t(col.matrix)
  col.matrix[lower.tri(col.matrix) & t(A.matrix) > 0] <- col
  tcm <- col.matrix
  col.matrix <- t(col.matrix)
  col.matrix[lower.tri(col.matrix)] <- tcm[lower.tri(tcm)]
  ##
  pos.matrix <- matrix(NA, nrow = n.trts, ncol = n.trts)
  dimnames(pos.matrix) <- dimnames(A.matrix)
  pos.matrix <- t(pos.matrix)
  pos.matrix[lower.tri(pos.matrix) & t(A.matrix) > 0] <- pos.number.of.studies
  tam <- pos.matrix
  pos.matrix <- t(pos.matrix)
  pos.matrix[lower.tri(pos.matrix)] <- tam[lower.tri(tam)]
  ##
  A.matrix <- A.matrix[seq1, seq1]
  N.matrix <- N.matrix[seq1, seq1]
  D.matrix <- D.matrix[seq1, seq1]
  ##
  col.matrix <- col.matrix[seq1, seq1]
  pos.matrix <- pos.matrix[seq1, seq1]
  ##
  if (thick == "matrix")
    W.matrix <- W.matrix[seq1, seq1]
  ##
  labels <- labels[seq1]
  
  
  A.sign <- sign(A.matrix)
  
  
  if ((is_2d & (is.null(xpos) & is.null(ypos))) |
      (is_3d & (is.null(xpos) & is.null(ypos) & is.null(zpos)))) {
    stressdata <- stress(x,
                         A.matrix = A.matrix,
                         N.matrix = N.matrix,
                         D.matrix = D.matrix,
                         ##
                         dim = dim,
                         start.layout = start.layout,
                         iterate = iterate,
                         eig1 = eig1, eig2 = eig2, eig3 = eig3,
                         tol = tol,
                         maxit = maxit,
                         ##
                         allfigures = allfigures,
                         ##
                         seq = seq,
                         ##
                         labels = labels,
                         cex = cex,
                         col = col,
                         adj = adj,
                         offset = offset,
                         scale = scale,
                         ##
                         plastic = plastic,
                         thickness = thickness,
                         lwd = lwd,
                         lwd.min = lwd.min,
                         lwd.max = lwd.max,
                         ##
                         highlight = highlight,
                         col.highlight = col.highlight,
                         lwd.highlight = lwd.highlight,
                         ## multiarm
                         col.multiarm = col.multiarm,
                         alpha.transparency = alpha.transparency,
                         ##
                         points = points, col.points = col.points,
                         cex.points = cex.points, pch.points = pch.points,
                         ##
                         number.of.studies = number.of.studies,
                         cex.number.of.studies = cex.number.of.studies,
                         col.number.of.studies = col.number.of.studies,
                         bg.number.of.studies = bg.number.of.studies,
                         pos.number.of.studies = pos.number.of.studies,
                         ##
                         ...)
    ##
    xpos <- stressdata$x
    ypos <- stressdata$y
    if (is_3d)
      zpos <- stressdata$z
  }
  
  
  if (allfigures)
    return(invisible(NULL))
  
  
  n <- dim(A.matrix)[1]
  d <- scale * max(abs(c(min(c(xpos, ypos), na.rm = TRUE),
                         max(c(xpos, ypos), na.rm = TRUE))))
  

  ##
  ##
  ## Generate datasets for plotting
  ##
  ##
  ##
  ## Dataset for nodes
  ##
  dat.nodes <- data.frame(labels, seq,
                          xpos, ypos, zpos = NA,
                          xpos.labels = NA, ypos.labels = NA,
                          cex = cex.points,
                          col = col.points,
                          pch = pch.points,
                          stringsAsFactors = FALSE)
  if (is_2d)
    dat.nodes$zpos <- NULL
  else {
    dat.nodes$zpos <- zpos
    dat.nodes$zpos.labels <- NA
  }
  ##
  if (is.null(adj)) {
    dat.nodes$adj.x <- NA
    dat.nodes$adj.y <- NA
    if (!is_2d)
      dat.nodes$adj.z <- NA
    ##
    dat.nodes$adj.x[dat.nodes$xpos >= 0] <- 0
    dat.nodes$adj.x[dat.nodes$xpos <  0] <- 1
    ##
    dat.nodes$adj.y[dat.nodes$ypos >  0] <- 0
    dat.nodes$adj.y[dat.nodes$ypos <= 0] <- 1
    ##
    if (!is_2d) {
      dat.nodes$adj.z[dat.nodes$zpos >  0] <- 0
      dat.nodes$adj.z[dat.nodes$zpos <= 0] <- 1
    }
  }
  else {
    dat.nodes$adj.x <- NA
    dat.nodes$adj.y <- NA
    ##
    if (length(adj) == 1) {
      dat.nodes$adj.x <- adj
      dat.nodes$adj.y <- adj
      if (!is_2d)
        dat.nodes$adj.z <- adj
    }
    else if (length(adj) == 2) {
      dat.nodes$adj.x <- adj[1]
      dat.nodes$adj.y <- adj[2]
      if (!is_2d)
        dat.nodes$adj.z <- 0.5
    }
    else if (length(adj) == 3 & !is_2d) {
      dat.nodes$adj.x <- adj[1]
      dat.nodes$adj.y <- adj[2]
      dat.nodes$adj.z <- adj[3]
    }
    else if (is.vector(adj)) {
        if (length(adj) != length(labels))
          stop("Length of vector 'adj' must be equal to number of treatments.")
        ##
        names(adj) <- x$trts
        dat.nodes$adj.x <- adj[seq1]
        dat.nodes$adj.y <- adj[seq1]
        ##
        if (!is_2d)
          dat.nodes$adj.z <- adj[seq1]
    }
    else if (is.matrix(adj)) {
      if (nrow(adj) != length(labels))
        stop("Number of rows of matrix 'adj' must be equal to number of treatments.")
      rownames(adj) <- x$trts
      dat.nodes$adj.x <- adj[seq1, 1]
      dat.nodes$adj.y <- adj[seq1, 2]
      ##
      if (!is_2d & ncol(adj) >= 3)
        dat.nodes$adj.z <- adj[seq1, 3]
    }
  }
  ##
  if (is_2d) {
    offset <- offset * 2 * d
    ##
    if (length(offset) == 1) {
      offset.x <- offset
      offset.y <- offset
    }
    else if (length(offset) == 2) {
      offset.x <- offset[1]
      offset.y <- offset[2]
    }
    else if (is.vector(offset)) {
      if (length(offset) != length(labels))
        stop("Length of vector 'offset' must be equal to number of treatments.")
      ##
      rownames(offset) <- x$trts
      offset.x <- offset[seq1]
      offset.y <- offset[seq1]
    }
    else if (is.matrix(adj)) {
      if (nrow(offset) != length(labels))
        stop("Number of rows of matrix 'offset' must be equal to number of treatments.")
      ##
      rownames(offset) <- x$trts
      offset.x <- offset[seq1, 1]
      offset.y <- offset[seq1, 2]
    }
    ##
    dat.nodes$xpos.labels <- dat.nodes$xpos - offset.x +
      2 * (dat.nodes$adj.x == 0) * offset.x
    dat.nodes$ypos.labels <- dat.nodes$ypos - offset.y +
      2 * (dat.nodes$adj.y == 0) * offset.y
  }
  else {
    dat.nodes$xpos.labels <- dat.nodes$xpos
    dat.nodes$ypos.labels <- dat.nodes$ypos
    dat.nodes$zpos.labels <- dat.nodes$zpos
  }
  ##
  ## Dataset for edges
  ##
  dat.edges <- data.frame(treat1 = rep("", n.edges),
                          treat2 = "",
                          n.stud = NA,
                          xpos = NA, ypos = NA,
                          adj = NA, pos.number.of.studies,
                          col = "",
                          stringsAsFactors = FALSE)
  ##
  comp.i <- 1
  ##
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (A.sign[i, j] > 0) {
        ##
        dat.edges$treat1[comp.i] <- rownames(A.matrix)[i]
        dat.edges$treat2[comp.i] <- colnames(A.matrix)[j]
        dat.edges$n.stud[comp.i] <- A.matrix[i, j]
        dat.edges$adj[comp.i] <- lambda <- pos.matrix[i, j]
        dat.edges$xpos[comp.i] <- lambda * xpos[i] + (1 - lambda) * xpos[j]
        dat.edges$ypos[comp.i] <- lambda * ypos[i] + (1 - lambda) * ypos[j]
        dat.edges$col[comp.i] <- col.matrix[i, j]
        ##
        comp.i <- comp.i + 1
      }
    }
  }
  
  
  ##
  ## Define coloured regions for multi-arm studies
  ##
  if (multiarm) {
    td1 <- data.frame(studies = x$studies, narms = x$narms)
    td1 <- td1[rev(order(td1$narms)), ]
    td1 <- td1[td1$narms > 2, ]
    multiarm.studies <- td1$studies
    ##
    n.multi <- length(multiarm.studies)
    ##
    missing.col.multiarm <- missing(col.multiarm)
    ##
    if (missing.col.multiarm | is.null(col.multiarm)) {
      ## Check for R package colorspace & use various gray values if
      ## not installed packages
      if (!any(as.data.frame(installed.packages())$Package == "colorspace"))
        col.polygon <- grDevices::rainbow(n.multi, alpha = alpha.transparency)
      else
        col.polygon <- colorspace::sequential_hcl(n.multi, alpha = alpha.transparency)
    }
    else {
      ##
      if (is.function(col.multiarm)) {
        mcname <- deparse(substitute(col.multiarm))
        ##
        csfun <- function(fcall, fname) {
          is.cs <- length(grep(fname, fcall)) > 0
          if (is.cs)
            meta:::is.installed.package("colorspace")
          is.cs
        }
        ##
        if (csfun(mcname, "rainbow_hcl"))
          col.polygon <- colorspace::rainbow_hcl(n.multi, start = 240, end = 60, alpha = alpha.transparency)
        else if (csfun(mcname, "sequential_hcl"))
          col.polygon <- colorspace::sequential_hcl(n.multi, alpha = alpha.transparency)
        else if (csfun(mcname, "diverge_hcl"))
          col.polygon <- colorspace::diverge_hcl(n.multi, alpha = alpha.transparency)
        else if (csfun(mcname, "heat_hcl"))
          col.polygon <- colorspace::heat_hcl(n.multi, alpha = alpha.transparency)
        else if (csfun(mcname, "terrain_hcl"))
          col.polygon <- colorspace::terrain_hcl(n.multi, alpha = alpha.transparency)
        else if (csfun(mcname, "diverge_hsv"))
          col.polygon <- colorspace::diverge_hsv(n.multi, alpha = alpha.transparency)
        else if (csfun(mcname, "choose_palette")) {
          fcolm <- colorspace::choose_palette(n = n.multi)
          col.polygon <- fcolm(n = n.multi)
        }
        else
          col.polygon <- sapply(n.multi, col.multiarm, alpha = alpha.transparency)
        ##
        if (csfun(mcname, "sequential_hcl") |
            csfun(mcname, "diverge_hcl") |
            csfun(mcname, "heat_hcl"))
          col.polygon <- rev(col.polygon)
      }
    }
    ##
    if (!missing.col.multiarm & is.character(col.multiarm)) {
      if (length(col.multiarm) > 1 & length(col.multiarm) != n.multi)
        stop("Length of argument 'col.multiarm' must be equal to one or the number of multi-arm studies: ", n.multi)
      col.polygon <- col.multiarm
    }
  }
  ##
  ## Define line width
  ##
  if (thick == "number.of.studies") {
    W.matrix <- lwd.max * A.matrix / max(A.matrix)
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "equal") {
    W.matrix <- lwd * A.sign
  }
  else if (thick == "se.fixed") {
    IV.matrix <- x$seTE.direct.fixed[seq1, seq1]
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max * min(IV.matrix, na.rm = TRUE) / IV.matrix
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "se.random") {
    IV.matrix <- x$seTE.direct.random[seq1, seq1]
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max * min(IV.matrix, na.rm = TRUE) / IV.matrix
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "w.fixed") {
    IV.matrix <- 1 / x$seTE.direct.fixed[seq1, seq1]^2
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max * IV.matrix / max(IV.matrix, na.rm = TRUE)
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "w.random") {
    IV.matrix <- 1 / x$seTE.direct.random[seq1, seq1]^2
    IV.matrix[is.infinite(IV.matrix)] <- NA
    W.matrix <- lwd.max * IV.matrix / max(IV.matrix, na.rm = TRUE)
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  else if (thick == "matrix") {
    W.matrix[is.infinite(W.matrix)] <- NA
    if (min(W.matrix[W.matrix != 0], na.rm = TRUE) == max(W.matrix[W.matrix != 0], na.rm = TRUE))
      W.matrix <- lwd * W.matrix
    else
      W.matrix <- lwd.max * W.matrix / max(W.matrix, na.rm = TRUE)
    W.matrix[W.matrix < lwd.min & W.matrix != 0] <- lwd.min
  }
  
  
  
  
  
  ##
  ##
  ## Plot graph
  ##
  ##
  range <- c(-d, d)
  ##
  if (is_2d) {
    oldpar <- par(xpd = TRUE, pty = "s")
    on.exit(par(oldpar))
    ##
    plot(xpos, ypos,
         xlim = range, ylim = range,
         type = "n", axes = FALSE, bty = "n",
         xlab = "", ylab = "",
         ...)
    ##
    ## Add coloured regions for multi-arm studies
    ##
    if (multiarm) {
      ##
      if (n.multi > 0) {
        multiarm.labels <- vector("list", n.multi)
        if (length(col.polygon) == 1)
          col.polygon <- rep(col.polygon, n.multi)
        for (i in 1:n.multi) {
          treat1 <- x$treat1[x$studlab %in% multiarm.studies[i]]
          treat2 <- x$treat2[x$studlab %in% multiarm.studies[i]]
          multiarm.labels[[i]] <- sort(unique(c(treat2, treat1)))
          ##
          dat.multi <- dat.nodes[dat.nodes$labels %in% multiarm.labels[[i]], ]
          if (nrow(dat.multi) == 0)
            dat.multi <- dat.nodes[dat.nodes$labels %in% multiarm.labels[[i]], ]
          ##
          ## Clockwise ordering of polygon coordinates
          ##
          polysort <- function(x, y) {
            xnorm <- (x - mean(x)) / sd(x) # Normalise coordinate x
            ynorm <- (y - mean(y)) / sd(y) # Normalise coordinate y
            r <- sqrt(xnorm^2 + ynorm^2)   # Calculate polar coordinates
            cosphi <- xnorm / r
            sinphi <- ynorm / r
            s <- as.numeric(sinphi > 0) # Define angles to lie in [0, 2 * pi]
            phi <- acos(cosphi)
            alpha <- s * phi + (1 - s) * (2 * pi - phi)
            ##
            res <- order(alpha)
            res
          }
          ##
          dat.multi <- dat.multi[polysort(dat.multi$xpos, dat.multi$ypos), ]
          ##
          polygon(dat.multi$xpos, dat.multi$ypos,
                  col = col.polygon[i], border = NA)
        }
      }
    }
    ##
    ## Draw lines
    ##
    if (plastic) {
      n.plastic <- 30
      lwd.multiply <- rep(NA, n.plastic)
      cols <- cols.highlight <- rep("", n.plastic)
      j <- 0
      for (i in n.plastic:1) {
        j <- j + 1
        lwd.multiply[j] <- sin(pi * i / 2 / n.plastic)
        cols[j] <- paste("gray", round(100 * (1 - i / n.plastic)), sep = "")
        cols.highlight[j] <- paste("gray", round(100 * (1 - i / n.plastic)), sep = "")
      }
      if (substring(col.highlight, nchar(col.highlight)) %in% 1:4)
        col.highlight <- substring(col.highlight, 1, nchar(col.highlight) - 1)
      cols.highlight[1:12] <- rep(paste(col.highlight, 4:1, sep = ""), rep(3, 4))
      cols.highlight[13:15] <- rep(col.highlight, 3)
    }
    else {
      lwd.multiply <- 1
      cols <- col
      cols.highlight <- col.highlight
    }
    ##
    comp.i <- 1
    ##
    for (n.plines in 1:length(lwd.multiply)) {
      for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
          if (A.sign[i, j] > 0) {
            lines(c(xpos[i], xpos[j]), c(ypos[i], ypos[j]),
                  lwd = W.matrix[i, j] * lwd.multiply[n.plines],
                  col = col.matrix[i, j])
            ##
            comp.i <- comp.i + 1
          }
        }
      }
    }
    ##
    ## Add highlighted comparisons
    ##
    if (!is.null(highlight)) {
      for (high in highlight) {
        highs <- unlist(compsplit(high, split = highlight.split))
        if (length(highs) != 2)
          stop("Wrong format for argument 'highlight' (see helpfile of plotgraph command).")
        ##
        if (sum(dat.nodes$labels %in% highs) != 2)
          stop(paste("Argument 'highlight' must contain two of the following values (separated by \":\"):\n  ",
                     paste(paste("'", dat.nodes$labels, "'", sep = ""),
                           collapse = " - "), sep = ""))
        ##
        dat.high <- dat.nodes[dat.nodes$labels %in% highs, ]
        ##
        if (is_2d)
          for (n.plines in 1:length(lwd.multiply))
            lines(dat.high$xpos, dat.high$ypos,
                  lwd = W.matrix[labels == highs[1], labels == highs[2]] * lwd.multiply[n.plines],
                  col = cols.highlight[n.plines])
      }
    }
    ##
    ## Add points for labels
    ##
    if (points)
      points(xpos, ypos,
             pch = pch.points, cex = cex.points, col = col.points)
    ##
    ## Print treatment labels
    ##
    if (!is.null(labels))
      for (i in 1:n)
        text(dat.nodes$xpos.labels[i], dat.nodes$ypos.labels[i],
             labels = dat.nodes$labels[i],
             cex = cex,
             adj = c(dat.nodes$adj.x[i], dat.nodes$adj.y[i]))
    ##
    ## Print number of treatments
    ##
    if (number.of.studies) {
      comp.i <- 1
      ##
      for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
          if (A.sign[i, j] > 0) {
            ##
            shadowtext(dat.edges$xpos[comp.i],
                       dat.edges$ypos[comp.i],
                       labels = dat.edges$n.stud[comp.i],
                       cex = cex.number.of.studies,
                       col = col.number.of.studies,
                       bg = bg.number.of.studies)
            ##
            comp.i <- comp.i + 1
          }
        }
      }
    }
  }
  else {
    rgl::plot3d(xpos, ypos, zpos,
                size = 10, col = col.points, cex = cex.points,
                axes = FALSE, box = FALSE,
                xlab = "", ylab = "", zlab = "")
    ##
    ## Add points for labels
    ##
    if (points)
      rgl::points3d(xpos, ypos, zpos,
                    pch = pch.points, cex = cex.points, col = col.points)
    ##
    ## Print treatment labels
    ##
    if (!is.null(labels))
      for (i in 1:n)
        rgl::text3d(dat.nodes$xpos.labels[i], dat.nodes$ypos.labels[i],
                    dat.nodes$zpos.labels[i],
                    texts = dat.nodes$labels[i],
                    cex = cex,
                    adj = c(dat.nodes$adj.x[i], dat.nodes$adj.y[i]))
    ##
    ## Add highlighted comparisons
    ##
    if (!is.null(highlight)) {
      for (high in highlight) {
        highs <- unlist(compsplit(high, split = highlight.split))
        if (length(highs) != 2)
          stop("Wrong format for argument 'highlight' (see helpfile of plotgraph command).")
        ##
        if (sum(dat.nodes$labels %in% highs) != 2)
          stop(paste("Argument 'highlight' must contain two of the following values (separated by \":\"):\n  ",
                     paste(paste("'", dat.nodes$labels, "'", sep = ""),
                           collapse = " - "), sep = ""))
        ##
        dat.high <- dat.nodes[dat.nodes$labels %in% highs, ]
        ##
        rgl::lines3d(dat.high$xpos * (1 + 1e-4), dat.high$ypos * (1 + 1e-4),
                     dat.high$zpos * (1 + 1e-4),
                     lwd = W.matrix[labels == highs[1], labels == highs[2]],
                     col = col.highlight)
      }
    }
    ##
    ## Add coloured regions for multi-arm studies
    ##
    if (multiarm) {
      ##
      morethan3 <- FALSE
      ##
      if (n.multi > 0) {
        multiarm.labels <- vector("list", n.multi)
        if (length(col.polygon) == 1)
          col.polygon <- rep(col.polygon, n.multi)
        for (i in 1:n.multi) {
          treat1 <- x$treat1[x$studlab %in% multiarm.studies[i]]
          treat2 <- x$treat2[x$studlab %in% multiarm.studies[i]]
          multiarm.labels[[i]] <- sort(unique(c(treat2, treat1)))
          ##
          dat.multi <- dat.nodes[dat.nodes$labels %in% multiarm.labels[[i]], ]
          if (nrow(dat.multi) == 0)
            dat.multi <- dat.nodes[dat.nodes$labels %in% multiarm.labels[[i]], ]
          if (nrow(dat.multi) == 3)
            rgl::triangles3d(dat.multi$xpos, dat.multi$ypos, dat.multi$zpos,
                             col = col.polygon[i])
          else
            morethan3 <- TRUE
        }
      }
      if (morethan3)
        warning("Multi-arm studies with more than three treatments not shown in 3-D plot.")
    }
    ##
    ## Draw lines
    ##
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (A.sign[i, j] > 0) {
          rgl::lines3d(c(xpos[i], xpos[j]), c(ypos[i], ypos[j]), c(zpos[i], zpos[j]),
                       lwd = W.matrix[i, j],
                       col = col)
        }
      }
    }
  }


  is.zero <- function(x) abs(x) < .Machine$double.eps^0.75
  ##
  dat.nodes$xpos[is.zero(dat.nodes$xpos)] <- 0
  dat.nodes$ypos[is.zero(dat.nodes$ypos)] <- 0
  ##
  if (!is_2d) {
    dat.nodes$zpos[is.zero(dat.nodes$zpos)] <- 0
    ##
    dat.nodes$xpos.labels <- NULL
    dat.nodes$ypos.labels <- NULL
  }
  else {
    dat.nodes$xpos.labels[is.zero(dat.nodes$xpos.labels)] <- 0
    dat.nodes$ypos.labels[is.zero(dat.nodes$ypos.labels)] <- 0
  }
  ##
  dat.nodes$zpos.labels <- NULL
  
  
  dat.edges$xpos[is.zero(dat.edges$xpos)] <- 0
  dat.edges$ypos[is.zero(dat.edges$ypos)] <- 0
  
  
  invisible(list(nodes = dat.nodes, edges = dat.edges))
}
