va.image <- function(x,
                     design, Hpn, h1,
                     ##
                     size = abs(Hpn), labels = NULL,
                     cex = 1, col, bg = "grey", box.bg = TRUE,
                     min.size = 0.10,
                     xlab = NULL, ylab = NULL, cex.axis = 1.1) {
  
  
  design.comb <- rep(NA, length(as.character(design$comparison)))
  ##
  for (i in 1:(length(as.character(design$comparison)))) {
    if ((((design$narms)[h1$order])[i]) > 2)
      design.comb[i] <- paste(as.character(design$comparison)[h1$order][i],
                              as.character(design$design)[h1$order][i],
                              sep = "_")
    else
      design.comb[i] <-as.character(design$comparison)[h1$order][i]
  }
  ##
  sc <- max(strwidth(design.comb, "inches")) / par("cin")[2]
  oldpar <- par(mar = c(1, 1.5 + sc, 1.5 + sc,
                        1.5 + strwidth("-4", "inches") / par("cin")[2]))
  ##
  plot(0, type = "n", xlim = c(0.036, 0.963), ylim = c(0.036, 0.963),
       bty = "n", xlab = "", ylab = "", axes = FALSE)
  ##
  if (is.null(xlab) && !is.null(rownames(x))) xlab <- rownames(x)
  if (is.null(ylab) && !is.null(colnames(x))) ylab <- colnames(x)
  ##
  rect(0, 0, 1, 1, col = bg, border = NA)
  ##
  runit <- 1 / ncol(x)
  ##
  center.x <- seq(from = runit / 2, to = 1 - runit / 2, len = ncol(x))
  center.y <- seq(from = runit / 2, to = 1 - runit / 2, len = ncol(x))
  ##
  val.grid <- seq(from = min(x), to = max(x), len = length(col) + 1)
  val.grid <- val.grid[-length(val.grid)]
  if (min(x) == max(x)) val.grid <- val.grid[1]
  ##
  col[unlist(lapply(x[, 1], function(arg) sum(arg < val.grid)))]
  ##
  min.sv <- min(size)
  max.sv <- max(size)

  fg.col <- bg.col <- apply(col2rgb(col) / 255, 2,
                            function(arg) rgb(arg[1], arg[2], arg[3]))
  if (box.bg)
    fg.col <- apply(col2rgb(col) / 255 * 0.8, 2,
                    function(arg) rgb(arg[1], arg[2], arg[3]))
  
  
  for (i in 1:ncol(x)) {
    actual.size <-
      (sqrt((size[, i] - min.sv) / (max.sv - min.sv)) + min.size) /
      (1 + 4 * min.size) * runit / 2
    ##
    rect(center.x - runit / 2, center.y[i] - runit / 2,
         center.x + runit / 2, center.y[i] + runit / 2,
         col = bg.col[unlist(lapply(x[, i],
                                    function(arg) sum(arg >= val.grid)))],
         border = NA)
    ##
    if (box.bg)
      rect(center.x - actual.size, center.y[i] - actual.size,
           center.x + actual.size, center.y[i] + actual.size,
           col = fg.col[unlist(lapply(x[, i],
                                      function(arg) sum(arg >= val.grid)))],
           border = NA)
    if (!is.null(labels))
      text(center.x, center.y[i], labels[, i], cex = cex)
  }
  
  
  if (!is.null(ylab))
    axis(side = 2, at = center.y, labels = ylab, tick = FALSE,
         line = NA, las = 2, cex.axis = cex.axis)
  ##
  if (!is.null(xlab))
    axis(side = 3, at = center.x, labels = xlab, tick = FALSE,
         line = NA, las = 2, cex.axis = cex.axis)

  axis(2, at = center.x, lty = 0, las = 2,
       labels = design.comb[length(design.comb):1], cex.axis = 1)
  axis(3, at = center.y, lty = 0, las = 3,
       labels = design.comb, cex.axis = 1)
  
  oldpar
}


legend.col <- function(col, lev, oldpar) {
  n <- length(col)
  bx <- par()$usr
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  xx <- rep(box.cx, each = 2)
  oldpar <- c(par(xpd = TRUE), oldpar)
  ##
  for (i in 1:n) {
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * i),
            box.cy[1] + (box.sy * i),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
  }
  ##
  par(new = TRUE)
  ##
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "", xaxt = "n", xlab = "",
       frame.plot = FALSE, yaxs = "i")
  ##
  axis(side = 4, tick = FALSE, line = .25, las = 1)
  
  invisible(NULL)
}
