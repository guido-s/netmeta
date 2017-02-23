plot.netposet <- function(x,
                          sel.x = 1, sel.y = 2, sel.z = 3,
                          dim = "2d",
                          col = "darkblue", col.lines = 8,
                          cex = 1, pch = 3, lty = 1, lwd = 2,
                          tcex = 0.8,
                          ...) {
  
  meta:::chkclass(x, "netposet")
  
  
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
  
  
  p.matrix <- x$pscore.matrix
  
  if (is_2d) {
    plot(p.matrix[, sel.x],
         p.matrix[, sel.y],
         pch = pch, cex = cex,
         xlab = x$outcomes[sel.x],
         ylab = x$outcomes[sel.y],
         xlim = c(0, 1), ylim = c(0, 1))
    
    for (i in seq_along(x$treatments)) {
      lines(x = c(p.matrix[i, sel.x], p.matrix[i, sel.x]),
            y = c(0, p.matrix[i,sel.y]),
            col = col.lines)
      ##
      lines(x = c(0, p.matrix[i, sel.x]),
            y = c(p.matrix[i, sel.y], p.matrix[i, sel.y]),
            col = col.lines)
      ##
      for (j in seq_along(x$treatments)) {
        if (x$M0[i, j] == 1) {
          lines(x = c(p.matrix[i, sel.x], p.matrix[j, sel.x]),
                y = c(p.matrix[i, sel.y], p.matrix[j, sel.y]),
                lwd = lwd)
        }
      }
    }
    ##
    text(p.matrix[, sel.x], p.matrix[, sel.y],
         rownames(p.matrix), col = col, cex = tcex)
  }
  else {
    rgl::plot3d(p.matrix[, sel.x], p.matrix[, sel.y], p.matrix[, sel.z], 
                xlab = x$outcomes[sel.x],
                ylab = x$outcomes[sel.y],
                zlab = x$outcomes[sel.z])
    ##
    for (i in seq_along(x$treatments))
      for (j in seq_along(x$treatments))
        if (x$M0[i, j] == 1)
          rgl::lines3d(x = c(p.matrix[i, sel.x], p.matrix[j, sel.x]),
                       y = c(p.matrix[i, sel.y], p.matrix[j, sel.y]),
                       z = c(p.matrix[i, sel.z], p.matrix[j, sel.z]),
                       lwd = 2)
  }
  
  invisible(NULL)
}
