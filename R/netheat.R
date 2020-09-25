#' Net heat plot
#' 
#' @description
#' This function creates a net heat plot, a graphical tool for
#' locating inconsistency in network meta-analyses.
#' 
#' @param x An object of class \code{netmeta}.
#' @param random A logical indicating whether the net heat plot should
#'   be based on a random effects model.
#' @param tau.preset An optional value for the square-root of the
#'   between-study variance \eqn{\tau^2} for a random effects model on
#'   which the net heat plot will be based.
#' @param showall A logical indicating whether results should be shown
#'   for all designs or only a sensible subset, see Details.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param \dots Additional arguments.
#' 
#' @details
#' The net heat plot is a matrix visualization proposed by Krahn et
#' al. (2013) that highlights hot spots of inconsistency between
#' specific direct evidence in the whole network and renders
#' transparent possible drivers.
#' 
#' In this plot, the area of a gray square displays the contribution
#' of the direct estimate of one design in the column to a network
#' estimate in a row. In combination, the colors show the detailed
#' change in inconsistency when relaxing the assumption of consistency
#' for the effects of single designs. The colors on the diagonal
#' represent the inconsistency contribution of the corresponding
#' design. The colors on the off-diagonal are associated with the
#' change in inconsistency between direct and indirect evidence in a
#' network estimate in the row after relaxing the consistency
#' assumption for the effect of one design in the column. Cool colors
#' indicate an increase and warm colors a decrease: the stronger the
#' intensity of the color, the greater the difference between the
#' inconsistency before and after the detachment. So, a blue colored
#' element indicates that the evidence of the design in the column
#' supports the evidence in the row. A clustering procedure is applied
#' to the heat matrix in order to find warm colored hot spots of
#' inconsistency. In the case that the colors of a column
#' corresponding to design \eqn{d} are identical to the colors on the
#' diagonal, the detaching of the effect of design \eqn{d} dissolves
#' the total inconsistency in the network.
#' 
#' The pairwise contrasts corresponding to designs of three- or
#' multi-arm studies are marked by '_' following the treatments of the
#' design.
#' 
#' By default (\code{showall = FALSE}), designs where only one
#' treatment is involved in other designs of the network or where the
#' removal of corresponding studies would lead to a splitting of the
#' network do not contribute to the inconsistency assessment and are
#' not incorporated into the net heat plot.
#' 
#' In the case of \code{random = TRUE}, the net heat plot is based on
#' a random effects model generalised for multivariate meta-analysis
#' in which the between-study variance \eqn{\tau^2} is estimated by the
#' method of moments (see Jackson et al., 2012) and embedded in a full
#' design-by-treatment interaction model (see Higgins et al., 2012).
#'
#' @author Ulrike Krahn \email{ulrike.krahn@@bayer.com}
#'
#' @seealso \link{netmeta}
#' 
#' @references
#' Krahn U, Binder H, König J (2013):
#' A graphical tool for locating inconsistency in network meta-analyses.
#' \emph{BMC Medical Research Methodology},
#' \bold{13}, 35
#' 
#' Jackson D, White IR and Riley RD (2012):
#' Quantifying the impact of between-study heterogeneity in
#' multivariate meta-analyses.
#' \emph{Statistics in Medicine},
#' \bold{31}, 3805--20
#' 
#' Higgins JPT, Jackson D, Barrett JK, Lu G, Ades AE, White IR (2012):
#' Consistency and inconsistency in network meta-analysis: concepts
#' and models for multi-arm studies.
#' \emph{Research Synthesis Methods},
#' \bold{3}, 98--110
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Senn2013)
#' 
#' # Generation of an object of class 'netmeta' with reference
#' # treatment 'plac', i.e. placebo
#' #
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'                 data = Senn2013, sm = "MD", reference = "plac")
#'         
#' # Generate a net heat plot based on a fixed effects model
#' #
#' netheat(net1) 
#' 
#' \dontrun{
#' # Generate a net heat plot based on a random effects model
#' #
#' netheat(net1, random = TRUE)
#' }
#'  
#' @export netheat


netheat <- function(x, random = FALSE, tau.preset = NULL,
                    showall = FALSE,
                    nchar.trts = x$nchar.trts,
                    ...) {
  
  
  meta:::chkclass(x, "netmeta")
  ##
  meta:::chklogical(random)
  meta:::chklogical(showall)
  meta:::chknumeric(nchar.trts, min = 1, length = 1)
  
  
  if (is.null(x$nchar.trts))
    nchar.trts <- 666
  ##
  trts <- colnames(x$A.matrix)
  trts.abbr <- treats(trts, nchar.trts)
  ##
  x$treat1 <- as.character(factor(x$treat1,
                                  levels = trts,
                                  labels = trts.abbr))
  ##
  x$treat2 <- as.character(factor(x$treat2,
                                  levels = trts,
                                  labels = trts.abbr))
  ##
  colnames(x$A.matrix) <- trts.abbr
  ##
  if (x$reference.group != "")
    x$reference.group <- trts.abbr[trts == x$reference.group]
  
  
  if (random == FALSE & length(tau.preset) == 0) {
    nmak <- nma.krahn(x)
    if (is.null(nmak)) {
      warning("Only a single design in network meta-analysis.",
              call. = FALSE)
      return(invisible(NULL))
    }
    decomp <- decomp.design(x)
    residuals <- decomp$residuals.inc.detach
    Q.inc.design <- decomp$Q.inc.design
  }
  ##
  if (length(tau.preset) == 1) {
    nmak <- nma.krahn(x, tau.preset = tau.preset)
    if (is.null(nmak)) {
      warning("Only a single design in network meta-analysis.",
              call. = FALSE)
      return(invisible(NULL))
    }
    decomp <- decomp.design(x, tau.preset = tau.preset)
    residuals <- decomp$residuals.inc.detach.random.preset
    Q.inc.design <- decomp$Q.inc.design.random.preset
  }
  ##
  if (random == TRUE & length(tau.preset) == 0) {
    tau.within <- tau.within(x)
    nmak <- nma.krahn(x, tau.preset = tau.within)
    if (is.null(nmak)) {
      warning("Only a single design in network meta-analysis.",
              call. = FALSE)
      return(invisible(NULL))
    }
    decomp <- decomp.design(x, tau.preset = tau.within)
    residuals <- decomp$residuals.inc.detach.random.preset
    Q.inc.design <- decomp$Q.inc.design.random.preset
  }
  
  
  if (nmak$d <= 2) {
    warning("Net heat plot not available due to small number of designs: ",
            nmak$d,
            call. = FALSE)
    return(invisible(NULL))
  }
  
  
  if (!showall)
    drop.designs <-
      names(Q.inc.design)[abs(Q.inc.design) <= .Machine$double.eps^0.5]
  
  
  H <- nmak$H
  V <- nmak$V
  ##
  if (!showall)
    design <- nmak$design[!(nmak$design$design %in% drop.designs), ]
  else
    design <- nmak$design
  
  
  df.Q.between.designs <- decomp$Q.decomp["Between designs", "df"]
  ##
  if (df.Q.between.designs == 0) {
    warning("Net heat plot not available because the network ",
            "does not contain any loop",
            if (any(x$narms > 2)) " of independent studies." else ".",
            call. = FALSE)
    return(invisible(NULL))
  }
  else if (df.Q.between.designs == 1) {
    warning("Net heat plot not available because detaching single designs ",
            "leads to a network without loops.",
            call. = FALSE)
    return(invisible(NULL))
  }
  else if (df.Q.between.designs > 1 & !any(!is.na(residuals))) {
    warning("Net heat plot not available because detaching ",
            "single designs leads to a insufficiently connected network.",
            call. = FALSE)
    return(invisible(NULL))
  }
  ##
  Q.inc.design.typ <- apply(residuals, 2, function(x) t(x) %*% solve(V) * x)
  inc <- matrix(Q.inc.design,
                nrow = nrow(Q.inc.design.typ),
                ncol = ncol(Q.inc.design.typ))
  diff <- inc - Q.inc.design.typ
  colnames(diff) <- colnames(Q.inc.design.typ)
  rownames(diff) <- rownames(residuals)
  ##
  if (!showall)
    diff <- diff[!(rownames(diff) %in% drop.designs),
                 !(colnames(diff) %in% drop.designs), drop = FALSE]
  ##
  if (all(is.na(diff))) {
    warning("Net heat plot not available as ",
            "no between-design heterogeneity exists.",
            call. = FALSE)
    return(invisible(NULL))
  }
  ##
  if (length(diff) == 1) {
    warning("Net heat plot not available due to small number of ",
            "informative designs.",
            call. = FALSE)
    return(invisible(NULL))
  }
  
  
  t1 <- -diff
  wi <- which(apply(t1, 2, function(x) sum(is.na(x)) == nrow(t1)))
  if (length(wi) > 0) {
    t1 <- t1[-wi, -wi, drop = FALSE]
  }
  dmat <- t1
  d1 <- dist(dmat, method = "manhattan")
  d1 <- d1 + dist(t(dmat), method = "manhattan")
  h1 <- hclust(d1)
  t1 <- t1[h1$order, h1$order]
  
  
  tn <- matrix(NA, nrow = nrow(t1), ncol = nrow(t1))
  
  
  for (i in 1:(nrow(t1)))
    tn[i, ] <- t1[nrow(t1) - (i - 1), ]
  
  
  tn[tn < -8] <- -8
  tn[tn >  8] <-  8
  
  
  ncolo <- 40
  sat <- min(max(abs(tn)) / 8, 1)
  nblue <- ifelse((max(max(tn), 0) - min(min(tn), 0)) != 0,
                  round(max(max(tn), 0) /
                        (max(max(tn), 0) - min(min(tn), 0)) * ncolo),
                  0)
  bf <- min(max(max(tn), 0) / abs(min(min(tn), 0)), 1) * sat
  clev <- seq(from = 1, to = 1 - bf, len = nblue)
  heat.vec <- heat.colors(round((ncolo - nblue - 1) / max(sat, 0.00001)))
  mycol <- c(heat.vec[(length(heat.vec) -
                       (ncolo - nblue - 1) + 1):length(heat.vec)],
             rgb(1, 1, 1), rgb(clev, clev, 1))
  
  
  if (length(wi) > 0)
    Hp <- H[(as.character(design$comparison)[-wi]), -wi]
  else
    Hp <- H[as.character(design$comparison), ]
  ##
  if (!showall)
    Hp <- Hp[!(rownames(Hp) %in% drop.designs),
             !(colnames(Hp) %in% drop.designs), drop = FALSE]
  ##
  Hp <- Hp[h1$order, h1$order]
  
  
  Hpn <- matrix(NA, nrow = nrow(Hp), ncol = ncol(Hp))
  ##
  for (i in 1:(nrow(Hp)))
    Hpn[i, ] <- Hp[ncol(Hp) - (i - 1), ]
  ##
  Hpn <- t(Hpn)
  
  
  va.image <- function(x, size = abs(Hpn), labels = NULL,
                       cex = 1, col = mycol, bg = "grey", box.bg = TRUE,
                       min.size = 0.10,
                       xlab = NULL, ylab = NULL, cex.axis = 1.1, ...) {
    
    
    if (length(wi) > 0)
      design2 <- design[-wi, ]
    else
      design2 <- design
    ##
    design.comb <- rep(NA, length(as.character(design2$comparison)))
    ##
    for (i in 1:(length(as.character(design2$comparison)))) {
      if ((((design2$narms)[h1$order])[i]) > 2)
        design.comb[i] <- paste(as.character(design2$comparison)[h1$order][i],
                                as.character(design2$design)[h1$order][i],
                                sep = "_")
      else
        design.comb[i] <-as.character(design2$comparison)[h1$order][i]
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
  
  
  oldpar <- va.image(x = t(tn) + max(abs(tn)))
  
  
  on.exit(par(oldpar))
  
  
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
  
  
  if (min(round(t1, 10)) != max(round(t1, 10)))
    legend.col(rev(mycol),
               seq(from = max(-max(t1), -8), to = min(-min(t1), 8),
                   len = length(mycol)),
               oldpar)
  
  
  box(lwd = 1.1)
  
  
  invisible(NULL)
}
