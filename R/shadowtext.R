shadowtext <- function(x, y = NULL, labels,
                       col = 'white', bg = 'black',
                       theta = seq(pi / 4, 2 * pi, length.out = 8),
                       r = 0.1, ...) {

  ##
  ## R function from TeachingDemos package by Greg Snow
  ##
  
  xy <- xy.coords(x, y)
  xo <- r * strwidth('A')
  yo <- r * strheight('A')
  
  for (i in theta)
    text(xy$x + cos(i) * xo, xy$y + sin(i) * yo, labels, col = bg, ...)
  
  text(xy$x, xy$y, labels, col = col, ...)
}
