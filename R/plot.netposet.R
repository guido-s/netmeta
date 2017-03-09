plot.netposet <- function(x,
                          pooled = ifelse(x$comb.random, "random", "fixed"),
                          sel.x = 1, sel.y = 2, sel.z = 3,
                          dim = "2d",
                          cex = 1, col = "black",
                          adj.x = 0, adj.y = 1,
                          offset.x = 0.005, offset.y = -0.005,
                          arrows = TRUE,
                          col.lines = "black", lty.lines = 1, lwd.lines = 1,
                          length = 0.05,
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
  
  
  sel.x <- as.numeric(meta:::setchar(sel.x, seq_len(n.outcomes)))
  sel.y <- as.numeric(meta:::setchar(sel.y, seq_len(n.outcomes)))
  if (n.outcomes > 2)
    sel.z <- as.numeric(meta:::setchar(sel.z, seq_len(n.outcomes)))
  ##  
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
  meta:::chknumeric(cex, min = 0, zero = 0)
  meta:::chknumeric(adj.x)
  meta:::chknumeric(adj.y)
  meta:::chknumeric(offset.x)
  meta:::chknumeric(offset.y)
  meta:::chklogical(arrows)
  meta:::chknumeric(lty.lines, min = 0, zero = 0, single = TRUE)
  meta:::chknumeric(lwd.lines, min = 0, zero = 0, single = TRUE)
  meta:::chknumeric(length, min = 0, zero = 0, single = TRUE)
  meta:::chklogical(grid)
  meta:::chknumeric(lty.grid, min = 0, zero = 0, single = TRUE)
  meta:::chknumeric(lwd.grid, min = 0, zero = 0, single = TRUE)
  
  
  if (!(length(cex) %in% c(1, n.treatments)))
    stop("Argument 'cex' must be a single numeric or vector of length ",
         n.treatments, ".")
  ##
  if (!(length(col) %in% c(1, n.treatments)))
    stop("Argument 'col' must be a character string or vector of length ",
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
  
  
  if (length(cex) == 1)
    cex <- rep(cex, n.treatments)
  ##
  if (length(col) == 1)
    col <- rep(col, n.treatments)
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
  
  
  seq.treats <- seq_len(n.treatments)
  ##
  if (is_2d) {
    plot(p.matrix[, sel.x], p.matrix[, sel.y],
         type = "n",
         xlab = outcomes[sel.x], ylab = outcomes[sel.y],
         xlim = c(0, 1), ylim = c(0, 1),
         ...)
    
    if (grid) {
      for (i in seq.treats) {
        lines(x = c(p.matrix[i, sel.x], p.matrix[i, sel.x]),
              y = c(0, p.matrix[i,sel.y]),
              col = col.grid, lty = lty.grid, lwd = lwd.grid)
        ##
        lines(x = c(0, p.matrix[i, sel.x]),
              y = c(p.matrix[i, sel.y], p.matrix[i, sel.y]),
              col = col.grid, lty = lty.grid, lwd = lwd.grid)
      }
    }
    ##
    for (i in seq.treats) {
      for (j in seq.treats) {
        if (M0[i, j] == 1) {
          if (arrows)
            arrows(p.matrix[i, sel.x], p.matrix[i, sel.y],
                   p.matrix[j, sel.x], p.matrix[j, sel.y],
                   col = col.lines, lty = lty.lines, lwd = lwd.lines,
                   length = length)
          else
            lines(c(p.matrix[i, sel.x], p.matrix[j, sel.x]),
                  c(p.matrix[i, sel.y], p.matrix[j, sel.y]),
                   col = col.lines, lty = lty.lines, lwd = lwd.lines)
        }
      }
    }
    ##
    for (i in seq.treats)
      text(p.matrix[i, sel.x] + offset.x[i],
           p.matrix[i, sel.y] + offset.y[i],
           rownames(p.matrix)[i],
           adj = c(adj.x[i], adj.y[i]),
           col = col[i], cex = cex[i])
  }
  else {
    rgl::plot3d(p.matrix[, sel.x], p.matrix[, sel.y], p.matrix[, sel.z], 
                xlab = outcomes[sel.x],
                ylab = outcomes[sel.y],
                zlab = outcomes[sel.z])
    ##
    for (i in seq.treats)
      for (j in seq.treats)
        if (M0[i, j] == 1)
          rgl::lines3d(x = c(p.matrix[i, sel.x], p.matrix[j, sel.x]),
                       y = c(p.matrix[i, sel.y], p.matrix[j, sel.y]),
                       z = c(p.matrix[i, sel.z], p.matrix[j, sel.z]),
                       lwd = 2)
  }
  
  invisible(NULL)
}
