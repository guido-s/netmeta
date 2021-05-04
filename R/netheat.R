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
#'   for all designs or only a sensible subset (see Details).
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
#' Designs where only one treatment is involved in other designs of
#' the network or where the removal of corresponding studies would
#' lead to a splitting of the network do not contribute to the
#' inconsistency assessment. By default (\code{showall = TRUE}), these
#' designs are not incorporated into the net heat plot. If
#' \code{showall = FALSE}, additional designs with minimal
#' contribution to the inconsistency Q statistic are not incorporated
#' (i.e., designs with \code{abs(Q.inc.design)} \code{<=}
#' \code{.Machine$double.eps^0.5)}.).
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
                    showall = TRUE,
                    nchar.trts = x$nchar.trts,
                    ...) {
  
  
  meta:::chkclass(x, "netmeta")
  ##
  meta:::chklogical(random)
  missing.showall <- missing(showall)
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
    if (!is.null(attributes(decomp)$netmetabin))
      return(invisible(NULL))
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
    if (!is.null(attributes(decomp)$netmetabin))
      return(invisible(NULL))
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
    if (!is.null(attributes(decomp)$netmetabin))
      return(invisible(NULL))
    residuals <- decomp$residuals.inc.detach.random.preset
    Q.inc.design <- decomp$Q.inc.design.random.preset
  }
  
  
  if (nmak$d <= 2) {
    warning("Net heat plot not available due to small number of designs: ",
            nmak$d,
            call. = FALSE)
    return(invisible(NULL))
  }
  
  
  if (!showall) {
    sel.Q.inc.design <-
      unlist(lapply(split(Q.inc.design, names(Q.inc.design)),
                    function(x)
                      sum(abs(x)) < .Machine$double.eps^0.5))
    drop.designs <- names(sel.Q.inc.design)[sel.Q.inc.design]
  }
  
  
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
            if (any(x$narms > 2)) " along more than one design." else ".",
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
  ##
  if (length(d1) > 0)
    h1 <- hclust(d1)
  else {
    warning("Insufficient number of designs ",
            "(available or selected by the program  specification) ",
            "for a net heat plot.",
            call. = FALSE)
    return(invisible(NULL))
  }
  ##
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
  
  
  oldpar <- va.image(t(tn) + max(abs(tn)),
                     design, Hpn, h1, wi,
                     col = mycol)
  ##  
  on.exit(par(oldpar))
  
  
  if (min(round(t1, 10)) != max(round(t1, 10)))
    legend.col(rev(mycol),
               seq(from = max(-max(t1), -8), to = min(-min(t1), 8),
                   len = length(mycol)),
               oldpar)
  ##
  box(lwd = 1.1)
  
  
  invisible(NULL)
}
