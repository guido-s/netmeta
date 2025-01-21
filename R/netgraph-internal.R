multicols <- function(studies, narms, missing.col.multiarm,
                      col.multiarm, alpha.transparency) {
  ##
  ## Define coloured regions for multi-arm studies
  ##
  dat <- data.frame(studies = studies, narms = narms,
                    idx = seq_along(studies))
  dat <- dat[rev(order(dat$narms)), ]
  dat <- dat[dat$narms > 2, ]
  ##
  o <- order(dat$idx)
  multiarm.studies <- dat$studies
  ##
  n.multi <- length(multiarm.studies)
  ##
  missing.col.multiarm <- missing(col.multiarm)
  ##
  o.cols <- TRUE
  ##
  if (missing.col.multiarm | is.null(col.multiarm))
    cols <- sequential_hcl(n.multi, alpha = alpha.transparency)
  else {
    ##
    if (is.function(col.multiarm)) {
      mcname <- deparse(substitute(col.multiarm))
      ##
      colfun <- function(fcall, fname)
        length(grep(fname, fcall)) > 0
      ##
      if (colfun(mcname, "rainbow_hcl"))
        cols <-
        rainbow_hcl(n.multi, start = 240, end = 60, alpha = alpha.transparency)
      else if (colfun(mcname, "sequential_hcl"))
        cols <- sequential_hcl(n.multi, alpha = alpha.transparency)
      else if (colfun(mcname, "diverge_hcl"))
        cols <- diverge_hcl(n.multi, alpha = alpha.transparency)
      else if (colfun(mcname, "heat_hcl"))
        cols <- heat_hcl(n.multi, alpha = alpha.transparency)
      else if (colfun(mcname, "terrain_hcl"))
        cols <- terrain_hcl(n.multi, alpha = alpha.transparency)
      else if (colfun(mcname, "diverge_hsv"))
        cols <- diverge_hsv(n.multi, alpha = alpha.transparency)
      else if (colfun(mcname, "choose_palette")) {
        fcolm <- choose_palette(n = n.multi)
        cols <- fcolm(n = n.multi)
      }
      else
        cols <- sapply(n.multi, col.multiarm, alpha = alpha.transparency)
      ##
      if (colfun(mcname, "sequential_hcl") | colfun(mcname, "diverge_hcl") |
          colfun(mcname, "heat_hcl"))
        cols <- rev(cols)
    }
  }
  ##
  if (!missing.col.multiarm & is.character(col.multiarm)) {
    if (length(col.multiarm) > 1 & length(col.multiarm) != n.multi)
      stop("Length of argument 'col.multiarm' must be equal to one or ",
           "the number of multi-arm studies: ", n.multi)
    cols <- col.multiarm
    o.cols <- FALSE
  }
  
  
  res <- data.frame(studlab = multiarm.studies[o])
  ##
  if (o.cols)
    res$col <- cols[o]
  else
    res$col <- cols
  ##
  res
}
