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
  if (missing.col.multiarm | is.null(col.multiarm)) {
    ## Check for R package colorspace & use various gray values if
    ## not installed packages
    if (!any(as.data.frame(installed.packages())$Package == "colorspace"))
      cols <- grDevices::rainbow(n.multi, alpha = alpha.transparency)
    else
      cols <- colorspace::sequential_hcl(n.multi, alpha = alpha.transparency)
  }
  else {
    ##
    if (is.function(col.multiarm)) {
      mcname <- deparse(substitute(col.multiarm))
      ##
      csfun <- function(fcall, fname) {
        is.cs <- length(grep(fname, fcall)) > 0
        if (is.cs)
          is.installed.package("colorspace")
        is.cs
      }
      ##
      if (csfun(mcname, "rainbow_hcl"))
        cols <- colorspace::rainbow_hcl(n.multi,
                                        start = 240, end = 60,
                                        alpha = alpha.transparency)
      else if (csfun(mcname, "sequential_hcl"))
        cols <- colorspace::sequential_hcl(n.multi, alpha = alpha.transparency)
      else if (csfun(mcname, "diverge_hcl"))
        cols <- colorspace::diverge_hcl(n.multi, alpha = alpha.transparency)
      else if (csfun(mcname, "heat_hcl"))
        cols <- colorspace::heat_hcl(n.multi, alpha = alpha.transparency)
      else if (csfun(mcname, "terrain_hcl"))
        cols <- colorspace::terrain_hcl(n.multi, alpha = alpha.transparency)
      else if (csfun(mcname, "diverge_hsv"))
        cols <- colorspace::diverge_hsv(n.multi, alpha = alpha.transparency)
      else if (csfun(mcname, "choose_palette")) {
        fcolm <- colorspace::choose_palette(n = n.multi)
        cols <- fcolm(n = n.multi)
      }
      else
        cols <- sapply(n.multi, col.multiarm, alpha = alpha.transparency)
      ##
      if (csfun(mcname, "sequential_hcl") |
          csfun(mcname, "diverge_hcl") |
          csfun(mcname, "heat_hcl"))
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
