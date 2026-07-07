#' Scatter plot to evaluate local inconsistency in network meta-analysis
#' 
#' @description
#' Scatter plot to show the difference between direct and indirect evidence in
#' network meta-analysis to assess local inconsistency.
#'
#' @aliases plot.netsplit
#' 
#' @param x An object of class \code{netsplit}.
#' @param type A character string indicating which figure type is to be used.
#'   Either \code{"Bland-Altman"} or \code{"Wilson"}, can be abbreviated.
#'   See Details.
#' @param pooled A character string indicating whether results for the
#'   common (\code{"common"}) or random effects model (\code{"random"})
#'   should be plotted. Can be abbreviated.
#' @param subset An optional logical vector specifying a subset of
#'   comparisons to consider (must be of same length as the total number of
#'   comparisons).
#' @param only.reference A logical indicating whether only comparisons
#'   with the reference group should be considered.
#' @param level The level used to calculate confidence intervals for local
#'   local inconsistency.
#' @param vertical A logical indicating whether confidence intervals should be
#'   plotted vertically or horizontally; only considered if argument
#'   \code{type = "Wilson"}. The interpretation is the same.
#' @param size.points A single numeric defining the point size. If \code{NULL},
#'   point sizes are proportional to the proportion of indirect evidence.
#' @param points.only.inc A logical indicating whether (larger) points should
#'   be only shown for comparisons with local inconsistency.
#' @param col Colour of points and confidence intervals.
#' @param col.inc Colour of points and confidence intervals for comparisons
#'   with local inconsistency.
#' @param labels A logical value indicating whether to show comparisons labels,
#'   or an expression or character vector specifying the labels to use.
#' @param labels.only.inc A logical indicating whether labels are only shown
#'   for comparisons with local inconsistency.
#' @param size.labels A single numeric defining the size of comparison labels.
#' @param col.labels Colour of comparison labels.
#' @param col.labels.inc Colour of labels for comparisons with local
#'   inconsistency.
#' @param max_overlap An integer specifying the maximum number of overlapping
#'   labels allowed in the scatter plot.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names in comparisons.
#' @param text.pooled A character string used to describe the statistical model.
#' @param main Main title.
#' @param xlab A label for the x-axis.
#' @param ylab A label for the y-axis.
#' @param xlim The x limits (min, max) of the plot.
#' @param ylim The y limits (min, max) of the plot.
#' @param tag The text for the tag label (see \code{\link[ggplot2]{labs}}).
#' @param \dots Additional arguments (ignored).
#' 
#' @details
#' A scatter plot is drawn in the active graphics window.
#' 
#' By default (argument \code{type = "Bland-Altman"}), a variant of the
#' Bland-Altman plot is shown (Bland & Altman, 1995). In this case, the
#' network estimates are shown on the horizontal axis and the difference between
#' direct and indirect effect estimates on the vertical axis. The confidence
#' intervals are those from the differences.
#' 
#' If argument \code{type = "Wilson"}, a scatter plot with direct effect
#' estimates on the horizontal and indirect estimates on the vertical axis
#' is produced (Wilson et al., 2026). The confidence intervals are
#' calculated from the indirect (argument \code{vertical = TRUE}) or direct
#' effect estimates and the standard error of the difference between indirect
#' and direct effect estimates.
#' 
#' For both plot types, the dashed green line corresponds to comparisons with
#' equal direct and indirect effect estimates.
#'
#' @return
#' An object of class \code{"plot.netsplit"}. The object is a list containing
#' the following components:
#' \item{plot}{A ggplot2 object.}
#' \item{data}{Data set used to create the plot.}
#' 
#' @author Federico Bonofiglio \email{bonostat@@gmx.de},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{forest.netsplit}}
#' 
#' @references
#' Bland JM, Altman DG (1995):
#' Comparing methods of measurement: Why plotting difference against standard
#' method is misleading.
#' \emph{The Lancet},
#' \bold{346}, 1085--87
#' 
#' Wilson H, Schönstein A, Robson S, Bonofiglio F (2026):
#' A novel approach for visualizing local consistency in
#' network meta-analysis.
#' \emph{Research Synthesis Methods},
#' 1--15
#' 
#' @keywords hplot
#' 
#' @examples
#' # Examples: example(netsplit)
#' 
#' @method plot netsplit
#' @export

plot.netsplit <- function(x,
                          type = "Bland-Altman",
                          pooled = ifelse(x$x$random, "random", "common"),
                          subset = NULL,
                          only.reference = x$only.reference,
                          level = x$level.ma,
                          #
                          vertical = TRUE,
                          #
                          size.points = NULL,
                          points.only.inc = FALSE,
                          col = "gray50",
                          col.inc = "red",
                          #
                          labels = TRUE,
                          labels.only.inc = FALSE,
                          size.labels = 5,
                          col.labels = col,
                          col.labels.inc = col.inc,
                          max_overlap = 10,
                          #
                          nchar.trts = x$nchar.trts,
                          #
                          text.pooled,
                          main,
                          xlab = NULL,
                          ylab = NULL,
                          xlim = NULL,
                          ylim = NULL,
                          tag = NULL,
                          ...) {
  
  #
  #
  # (1) Check and set arguments
  #
  #
  
  chkclass(x, "netsplit")
  x <- updateversion(x)
  #
  type <- setchar(type, c("Wilson", "Bland-Altman"))
  #
  pooled <- setchar(pooled, c("common", "random", "fixed"))
  pooled[pooled == "fixed"] <- "common"
  #
  missing.only.reference <- missing(only.reference)
  if (!missing.only.reference)
    chklogical(only.reference)
  #
  if (missing(text.pooled)) {
    if (pooled == "common")
      text.pooled <- gsub(" effect ", " effects ", gs("text.common"))
    else
      text.pooled <- gs("text.random")
  }
  else if (!is.null(text.pooled))
    chkchar(text.pooled, length = 1)
  #
  chklevel(level, length = 1)
  chklogical(vertical)
  #
  same.size <- !is.null(size.points)
  if (same.size)
    chknumeric(size.points, min = 0, zero = TRUE, length = 1)
  #
  chklogical(points.only.inc)
  #
  chkcolor(col, length = 1)
  chkcolor(col.inc, length = 1)
  #
  chklogical(labels)
  chknumeric(size.labels, min = 0, zero = TRUE, length = 1)
  chkcolor(col.labels, length = 1)
  chkcolor(col.labels.inc, length = 1)
  chklogical(labels.only.inc)
  chknumeric(max_overlap, min = 0, zero = TRUE, length = 1)
  #
  nchar.trts <- replaceNULL(nchar.trts, 666)
  chknumeric(nchar.trts, length = 1)
  #
  if (missing(main)) {
    if (type == "Wilson")
      main <- paste0("Direct against indirect estimates with ",
                     round(100 * level, 1),
                     "% confidence intervals for their differences")
    else if (type == "Bland-Altman")
      main <- paste0("Bland-Altman plot with ",
                     round(100 * level, 1),
                     "% confidence intervals for differences")
  }
  else if (!is.null(main))
    chkchar(main, length = 1)
  #
  if (!is.null(xlab))
    chkchar(xlab, length = 1)
  if (!is.null(ylab))
    chkchar(ylab, length = 1)
  #
  if (!is.null(xlim))
    chknumeric(xlim, length = 2)
  if (!is.null(ylim))
    chknumeric(ylim, length = 2)
  if (!is.null(tag))
    chkchar(tag, length = 1)
  #
  calculate.ci <- level != x$level.ma
  #
  # Get rid of warning 'Undefined global functions or variables'
  #
  comparison <- TE <- prop.direct<- prop.indirect <-
    direct <- indirect <- mixed <- lower <- upper <-
    pval <- label_width <- Value <- Bound <- signif_label <- NULL
  
  
  #
  #
  # (2) Extract results for common or random effects model
  #
  #
  
  if (pooled == "common") {
    dat.mixed <- x$common
    #
    dat.direct <- x$direct.common
    #
    dat.indirect <- x$indirect.common %>% mutate(prop.direct = x$prop.common)
    #
    dat.diff <- x$compare.common
  }
  else {
    dat.mixed <- x$random
    #
    dat.direct <- x$direct.random
    #
    dat.indirect <- x$indirect.random %>% mutate(prop.direct = x$prop.random)
    dat.indirect$prop.indirect <- 1 - x$prop.random
    #
    dat.diff <- x$compare.random
  }
  #
  dat.indirect %<>% mutate(prop.indirect = 1 - prop.direct)
  #
  dat.mixed %<>%
    select(comparison, TE) %>%
    rename(mixed = TE)
  #
  dat.direct %<>%
    select(comparison, TE) %>%
    rename(direct = TE)
  #
  dat.indirect %<>%
    select(comparison, TE, prop.direct, prop.indirect) %>%
    rename(indirect = TE)
  #
  if (calculate.ci) {
    ci.diff <- ci(dat.diff$TE, dat.diff$seTE, level = level)
    dat.diff %<>% mutate(lower = ci.diff$lower, upper = ci.diff$upper)
  }
  #
  dat.diff %<>%
    select(comparison, TE, lower, upper, p) %>%
    rename(diff = TE, pval = p)
  #
  dat <- inner_join(dat.mixed, dat.direct, by = "comparison")
  dat <- inner_join(dat, dat.indirect, by = "comparison")
  dat <- inner_join(dat, dat.diff, by = "comparison")
  #
  if (type == "Bland-Altman") {
    tmp.lower <- dat$lower
    #
    dat %<>% mutate(direct = mixed, indirect = diff)
    #
    vertical <- TRUE
  }
  else {
    if (vertical) {
      dat %<>%
        mutate(lower = direct - lower, upper = direct - upper)
    }
    else {
      dat %<>%
        mutate(lower = indirect + lower,
               upper = indirect + upper)
    }
  }
  #
  dat %<>%
    mutate(signif =
             factor(if_else(pval < 1 - level, "Y", "N"), levels = c("N", "Y")),
           signif_label =
             factor(if_else(pval < 1 - level, "Yl", "Nl"),
                    levels = c("Nl", "Yl")))
    
  
  #
  #
  # (3) Select treatment comparisons to show in scatter plot
  #
  #
  
  if (!missing(subset)) {
    #
    # Catch subset from data:
    #
    mc <- match.call()
    sfsp <- sys.frame(sys.parent())
    #
    error <- try(subset.x <- catch("subset", mc, dat, sfsp), silent = TRUE)
    if (!any(class(error) == "try-error"))
      subset <- subset.x
    #
    if (!is.null(subset)) {
      chklength(subset, length(dat$comparison),
                text = paste0("Argument 'subset' must be of length ",
                              length(dat$comparison), "."))
      if (!is.logical(subset))
        stop("Argument 'subset' must be a logical vector.")
    }
    #
    dat <- dat %>% filter(subset)
  }
  #
  if (only.reference) {
    if (x$reference.group == "") {
      warning("First treatment used as reference as argument ",
              "'reference.group' was unspecified in netsplit().",
              call. = FALSE)
      x$reference.group <- compsplit(x$comparison, x$sep.trts)[[1]][1]
    }
    #
    sel <- apply(!is.na(sapply(compsplit(dat$comparison, x$sep.trts),
                               match, x$reference.group)), 2, sum) >= 1
    #
    dat <- dat %>% filter(sel)
  }
  #
  # Keep only comparisons for which we have both direct and indirect.
  #
  dat %<>% filter(!is.na(direct), !is.na(indirect))
  #
  any_sign <- any(dat$pval < 1 - level)
  #
  trts <- unique(sort(unlist(compsplit(dat$comparison, split = x$sep.trts))))
  trts.abbr <- treats(trts, nchar.trts = nchar.trts)
  #
  for (i in seq_len(nrow(dat))) {
    trts.i <- unlist(compsplit(dat$comparison[i], split = x$sep.trts))
    trts.i <- factor(trts.i, levels = trts, labels = trts.abbr)
    dat$comparison[i] <- paste(trts.i, collapse = x$sep.trts)
  }
  #
  # Show (larger) points only for comparisons with local inconsistency
  #
  if (same.size)
    size.points <- rep(size.points, nrow(dat))
  #
  all.inc <- all(dat$signif == "Y", na.rm = TRUE)
  all.noninc <- all(dat$signif == "N", na.rm = TRUE)
  #
  if (all.inc | all.noninc)
    points.only.inc <- FALSE
  #
  if (points.only.inc) {
    if (same.size)
      size.points[dat$signif == "N"] <- 0
    else
      #
      dat %<>% mutate(prop.indirect = if_else(signif == "Y", prop.indirect, 0))
  }
  #
  any.inc <- any(dat$signif == "Y", na.rm = TRUE)
  
  
  dat.lines <- dat %>%
    select(comparison, direct, indirect, lower, upper, signif, signif_label) %>%
    pivot_longer(cols = c(lower, upper),
                 names_to = "Bound", values_to = "Value") %>%
    mutate(label_width = str_length(comparison))
  
  
  #
  #
  # (4) Create scatter plot
  #
  #
  
  xylim <- max(abs(dat.lines$Value))
  axis.color <- "gray40"
  #
  if (vertical) {
    maps <- aes(x = direct, y = Value, group = comparison, color = signif)
    mapslab <- aes(x = direct, y = Value, label = comparison, color = signif_label)
  }
  else {
    maps <- aes(x = Value, y = indirect, group = comparison, color = signif)
    mapslab <- aes(x = Value, y = indirect, label = comparison, color = signif_label)
  }
  #
  if (same.size)
    mapspoint <- aes(x = direct, y = indirect, color = signif,
                     size = size.points)
  else
    mapspoint <- aes(x = direct, y = indirect, color = signif,
                     size = 100 * prop.indirect)
  #
  p <- ggplot()
  #
  if (type == "Wilson") {
    if (is.null(xlab)){
      xlab <- "Direct effect estimate"
      if (x$sm != "")
        xlab <- paste0(xlab, " (", tolower(meta:::xlab_meta(x$sm, FALSE)), ")")
    }
    #
    if (is.null(ylab)){
      ylab <- "Indirect effect estimate"
      if (x$sm != "")
        ylab <- paste0(ylab, " (", tolower(meta:::xlab_meta(x$sm, FALSE)), ")")
    }
    #
    p <- p +
      geom_abline(intercept = 0, slope = 1,
                  linetype = "dashed", color = "green4") +
      #
      geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
  }
  else {
    if (is.null(xlab)){
      xlab <- "Network estimate"
      if (x$sm != "")
        xlab <- paste0(xlab, " (", tolower(meta:::xlab_meta(x$sm, FALSE)), ")")
    }
    #
    if (is.null(ylab)){
      ylab <- "Difference between direct and indirect evidence"
      if (x$sm != "")
        ylab <- paste0(ylab, "\n(", tolower(meta:::xlab_meta(x$sm, FALSE)), ")")
    }
    #
    p <- p +
      geom_vline(xintercept = 0) +
      geom_abline(intercept = 0, slope = 0,
                  linetype = "dashed", color = "green4")
  }
  #
  p <- p +
    geom_line(data = dat.lines, mapping = maps) +
    geom_point(data = dat, mapping = mapspoint) +
    # Significance colors
    # If only significant data shown then stick to 'col.inc'
    scale_color_manual(
      values =
        c("N" = if (all(dat$pval < 1 - level)) col.inc else col,
          "Y" = col.inc,
          "Nl" = if (all(dat$pval < 1 - level)) col.labels.inc else col.labels,
          "Yl" = col.labels.inc)) +
    # Define theme
    theme(
      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.margin = margin(20, 40, 20, 20),
      # panel.background = element_rect(fill = alpha("white", 1)),
      axis.ticks = element_line(color = axis.color),
      axis.line = element_line(color = axis.color),
      # panel.grid = element_line(color = alpha("white", 0)),
      axis.text = element_text(color = axis.color),
      axis.title = element_text(size  = 9, color = axis.color),
      panel.background = element_rect(fill = "gray99"),
      panel.grid.major = element_line(color = "gray95"),
      panel.grid.minor = element_blank(),
      #
      legend.position = "bottom")
  #
  if (!same.size) {
    p <- p +
      guides(colour = "none", size = guide_legend(title = "Indirect (%)"))
  }
  else {
    p <- p + guides(colour = "none", size = "none")
  }
  #
  # Add plot labels
  #
  if (any_sign) {
    p <- p +
      labs(title = main, subtitle = text.pooled, x = xlab, y = ylab,
           caption = paste0(
             "Significant results (p < ", 1 - level, ") are displayed in ",
             col.inc),
           tag = tag)
  }
  else {
    p <- p +
      labs(title = main, subtitle = text.pooled, x = xlab, y = ylab, tag = tag)
  }
  #
  if (labels) {
    dat.labels <- dat.lines %>%
      filter(Bound == "lower")
    #
    if (labels.only.inc)
      dat.labels %<>% filter(signif == "Y")
    #
    p <- p +
      geom_text_repel(
        data = dat.labels,
        mapping = mapslab,
        size = size.labels,
        max.overlaps = max_overlap)
  }
  #
  if (!is.null(xlim) & is.null(ylim))
    p <- p + coord_cartesian(xlim = xlim)
  #
  else if (is.null(xlim) & !is.null(ylim))
    p <- p + coord_cartesian(ylim = ylim)
  #
  else if (!is.null(xlim) & !is.null(ylim))
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  #
  else if (type == "Wilson") {
    p <- p +
      coord_cartesian(xlim = c(-xylim, xylim), ylim = c(-xylim, xylim))
  }
  #
  res <- list(plot = p, data = dat)
  #
  if (any(trts != trts.abbr)) {
    res$trts <- trts
    res$trts.abbr <- trts.abbr
  }
  #
  class(res) <- "plot.netsplit"
  #
  res
}

#' @method print plot.netsplit
#' @export

print.plot.netsplit <- function(x, ...) {
  print(x$plot)
  #
  invisible(NULL)
}
