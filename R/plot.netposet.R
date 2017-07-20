plot.netposet <- function(x,
                          plottype = "scatter",
                          pooled = ifelse(x$comb.random, "random", "fixed"),
                          dim = "2d",
                          sel.x = 1, sel.y = 2, sel.z = 3,
                          cex.text = 1, col.text = "black",
                          adj.x = 0, adj.y = 1,
                          offset.x = 0.005, offset.y = -0.005,
                          col.lines = "black", lty.lines = 1, lwd.lines = 1,
                          arrows = FALSE,
                          length = 0.05,
                          pch = NULL, cex.points = 1, col.points = "black",
                          grid = TRUE,
                          col.grid = "gray", lty.grid = 2, lwd.grid = 1,
                          ...) {
  
  meta:::chkclass(x, "netposet")
  ##
  pooled <- meta:::setchar(pooled, c("fixed", "random"))
  
  
  if (pooled == "fixed") {
    p.matrix <- x$P.fixed
    M0 <- x$M0.fixed
  }
  else {
    p.matrix <- x$P.random
    M0 <- x$M0.random
  }
  ##
  outcomes <- colnames(p.matrix)
  treatments <- rownames(p.matrix)
  ##
  n.outcomes   <- length(outcomes)
  n.treatments <- length(treatments)
  
  
  dim <- meta:::setchar(dim, c("2d", "3d"))
  is_2d <- dim == "2d"
  is_3d <- !is_2d
  ##
  plottype <- meta:::setchar(plottype, c("scatter", "biplot"))
  is_biplot <- plottype == "biplot"
  ##
  if (is_biplot) {
    outcomes <- paste("Principal Component", 1:3)
    sel.x <- 2
    sel.y <- 1
    ##
    if (n.outcomes > 2 & is_3d) {
      sel.x <- 2
      sel.y <- 3
      sel.z <- 1
    }
    p.matrix <- prcomp(p.matrix, scale = TRUE)$x
  }
  else {
    sel.x <- as.numeric(meta:::setchar(sel.x, seq_len(n.outcomes)))
    sel.y <- as.numeric(meta:::setchar(sel.y, seq_len(n.outcomes)))
    ##
    if (n.outcomes > 2)
      sel.z <- as.numeric(meta:::setchar(sel.z, seq_len(n.outcomes)))
  }
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
  use_pch <- !is.null(pch)
  ##
  meta:::chknumeric(cex.text, min = 0, zero = 0)
  meta:::chknumeric(adj.x)
  meta:::chknumeric(adj.y)
  meta:::chknumeric(offset.x)
  meta:::chknumeric(offset.y)
  meta:::chknumeric(lty.lines, min = 0, zero = 0, single = TRUE)
  meta:::chknumeric(lwd.lines, min = 0, zero = 0, single = TRUE)
  meta:::chklogical(arrows)
  meta:::chknumeric(length, min = 0, zero = 0, single = TRUE)
  if (use_pch)
    meta:::chknumeric(cex.points, min = 0, zero = 0)
  meta:::chklogical(grid)
  meta:::chknumeric(lty.grid, min = 0, zero = 0, single = TRUE)
  meta:::chknumeric(lwd.grid, min = 0, zero = 0, single = TRUE)
  
  
  if (!(length(cex.text) %in% c(1, n.treatments)))
    stop("Argument 'cex.text' must be a single numeric or vector of length ",
         n.treatments, ".")
  ##
  if (!(length(col.text) %in% c(1, n.treatments)))
    stop("Argument 'col.text' must be a character string or vector of length ",
         n.treatments, ".")
  ##
  if (!(length(adj.x) %in% c(1, n.treatments)))
    stop("Argument 'adj.x' must be a single numeric or vector of length ",
         n.treatments, ".")
  ##
  if (!(length(adj.y) %in% c(1, n.treatments)))
    stop("Argument 'adj.y' must be a single numeric or vector of length ",
         n.treatments, ".")
  ##
  if (!(length(offset.x) %in% c(1, n.treatments)))
    stop("Argument 'offset.x' must be a single numeric or vector of length ",
         n.treatments, ".")
  ##
  if (!(length(offset.y) %in% c(1, n.treatments)))
    stop("Argument 'offset.y' must be a single numeric or vector of length ",
         n.treatments, ".")
  ##
  if (use_pch) {
    if (!(length(pch) %in% c(1, n.treatments)))
      stop("Argument 'pch' must be a single numeric or vector of length ",
           n.treatments, ".")
    ##
    if (!(length(cex.points) %in% c(1, n.treatments)))
      stop("Argument 'cex.points' must be a single numeric or vector of length ",
           n.treatments, ".")
    ##
    if (!(length(col.points) %in% c(1, n.treatments)))
      stop("Argument 'col.points' must be a character string or vector of length ",
           n.treatments, ".")
  }
  
  
  if (length(cex.text) == 1)
    cex.text <- rep(cex.text, n.treatments)
  ##
  if (length(col.text) == 1)
    col.text <- rep(col.text, n.treatments)
  ##
  if (length(adj.x) == 1)
    adj.x <- rep(adj.x, n.treatments)
  ##
  if (length(adj.y) == 1)
    adj.y <- rep(adj.y, n.treatments)
  ##
  if (length(offset.x) == 1)
    offset.x <- rep(offset.x, n.treatments)
  ##
  if (length(offset.y) == 1)
    offset.y <- rep(offset.y, n.treatments)
  ##
  if (use_pch) {
    if (length(pch) == 1)
      pch <- rep(pch, n.treatments)
    ##
    if (length(cex.points) == 1)
      cex.points <- rep(cex.points, n.treatments)
    ##
    if (length(col.points) == 1)
      col.points <- rep(col.points, n.treatments)
  }
  
  
  seq.treats <- seq_len(n.treatments)
  ##
  if (is_2d) {
    ##
    ## 2-D plot
    ##
    xvals <- p.matrix[, sel.x]
    yvals <- p.matrix[, sel.y]
    ##
    if (is_biplot) {
      plot(xvals, yvals,
           type = "n",
           xlab = outcomes[sel.x], ylab = outcomes[sel.y],
           ...)
    }
    else {
      ##
      plot(xvals, yvals,
           type = "n",
           xlab = outcomes[sel.x], ylab = outcomes[sel.y],
           xlim = c(0, 1), ylim = c(0, 1),
           ...)
    }
    
    if (grid & !is_biplot) {
      for (i in seq.treats) {
        lines(x = c(xvals[i], xvals[i]),
              y = c(0, yvals[i]),
              col = col.grid, lty = lty.grid, lwd = lwd.grid)
        ##
        lines(x = c(0, xvals[i]),
              y = c(yvals[i], yvals[i]),
              col = col.grid, lty = lty.grid, lwd = lwd.grid)
      }
    }
    ##
    if (use_pch)
      for (i in seq.treats)
        points(xvals[i], yvals[i],
               pch = pch[i], col = col.points[i], cex = cex.points[i])
    ##
    for (i in seq.treats) {
      for (j in seq.treats) {
        if (M0[i, j] == 1) {
          if (arrows)
            arrows(xvals[i], yvals[i],
                   xvals[j], yvals[j],
                   col = col.lines, lty = lty.lines, lwd = lwd.lines,
                   length = length)
          else
            lines(c(xvals[i], xvals[j]),
                  c(yvals[i], yvals[j]),
                   col = col.lines, lty = lty.lines, lwd = lwd.lines)
        }
      }
    }
    ##
    for (i in seq.treats)
      text(xvals[i] + offset.x[i],
           yvals[i] + offset.y[i],
           treatments[i],
           adj = c(adj.x[i], adj.y[i]),
           col = col.text[i], cex = cex.text[i])
  }
  else {
    ##
    ## 3-D plot
    ##
    xvals <- p.matrix[, sel.x]
    yvals <- p.matrix[, sel.y]
    zvals <- p.matrix[, sel.z]
    ##
    rgl::plot3d(xvals, yvals, zvals,
                xlab = outcomes[sel.x],
                ylab = outcomes[sel.y],
                zlab = outcomes[sel.z])
    ##
    if (use_pch)
      for (i in seq.treats)
        rgl::points3d(xvals[i], yvals[i], zvals[i],
                      pch = pch[i], col = col.points[i], cex = cex.points[i])
    ##
    for (i in seq.treats)
      for (j in seq.treats)
        if (M0[i, j] == 1)
          rgl::lines3d(x = c(xvals[i], xvals[j]),
                       y = c(yvals[i], yvals[j]),
                       z = c(zvals[i], zvals[j]),
                       lwd = 2)
    ##
    for (i in seq.treats)
      rgl::text3d(xvals[i], yvals[i], zvals[i],
                  treatments[i],
                  col = col.text[i], cex = cex.text[i])
  }
  
  invisible(NULL)
}
