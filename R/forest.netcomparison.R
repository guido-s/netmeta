#' Forest plot for complex interventions in component network
#' meta-analysis
#' 
#' @description
#' Draws a forest plot in the active graphics window (using grid
#' graphics system).
#'
#' @aliases forest.netcomparison plot.netcomparison
#' 
#' @param x An object of class \code{netcomparison}.
#' @param pooled A character string indicating whether results for the
#'   common (\code{"common"}) or random effects model (\code{"random"})
#'   should be plotted. Can be abbreviated.
#' @param leftcols A character vector specifying (additional) columns
#'   to be plotted on the left side of the forest plot or a logical
#'   value (see \code{\link[meta]{forest.meta}} help page for details).
#' @param leftlabs A character vector specifying labels for
#'   (additional) columns on left side of the forest plot (see
#'   \code{\link[meta]{forest.meta}} help page for details).
#' @param rightcols A character vector specifying (additional) columns
#'   to be plotted on the right side of the forest plot or a logical
#'   value (see \code{\link[meta]{forest.meta}} help page for details).
#' @param rightlabs A character vector specifying labels for
#'   (additional) columns on right side of the forest plot (see
#'   \code{\link[meta]{forest.meta}} help page for details).
#' @param nchar.comps A numeric defining the minimum number of
#'   characters used to create unique names for components.
#' @param digits Minimal number of significant digits for treatment
#'   effects and confidence intervals, see \code{print.default}.
#' @param digits.stat Minimal number of significant digits for tests
#'   of overall effect, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall effects, see \code{print.default}.
#' @param smlab A label printed at top of figure. By default, text
#'   indicating either common or random effects model is printed.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in forest plots. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios, for example.
#' @param lab.NA A character string to label missing values.
#' @param equal.size A logical indicating whether all squares should
#'   be of equal size. Otherwise, the square size is proportional to
#'   the precision of estimates.
#' @param \dots Additional arguments for \code{\link[meta]{forest.meta}}
#'   function.
#' 
#' @details
#' A forest plot, also called confidence interval plot, is drawn in
#' the active graphics window. For more information see help page of
#' \code{\link[meta]{forest.meta}} function.
#'
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netcomparison}}, \code{\link{netcomb}},
#'   \code{\link{discomb}}, \code{\link[meta]{forest.meta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Linde2016)
#' 
#' # Only consider studies including Face-to-face PST (to reduce
#' # runtime of example)
#' #
#' face <- subset(Linde2016, id %in% c(16, 24, 49, 118))
#' 
#' # Conduct random effects network meta-analysis
#' #
#' net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'   data = face, ref = "placebo", sm = "OR", common = FALSE)
#' 
#' # Additive model for treatment components (with placebo as inactive
#' # treatment)
#' #
#' nc1 <- netcomb(net1, inactive = "placebo")
#'
#' # Some comparisons
#' #
#' t1 <- c("F + TCA", "F + Plac", "SSRI + Plac + TCA")
#' t2 <- c("UC", "Plac", "UC")
#' #
#' netcomparison(nc1, t1, t2)
#' #
#' forest(netcomparison(nc1, t1, t2))
#' forest(netcomparison(nc1, t1, t2), nchar.comps = 4)
#' forest(netcomparison(nc1, c("F", "TCA"), "UC"), nchar.comps = 4)
#' 
#' @method forest netcomparison
#' @export

forest.netcomparison <- function(x,
                                 pooled =
                                   ifelse(x$random, "random", "common"),
                                 leftcols = c("studlab", "treat1", "treat2"),
                                 leftlabs = c("Comparison", "Trt 1", "Trt 2"),
                                 rightcols =
                                   c("effect", "ci", "statistic", "pval"),
                                 rightlabs = c(NA, NA, "z", "p-value"),
                                 nchar.comps = x$nchar.trts,
                                 digits = gs("digits.forest"),
                                 digits.stat = gs("digits.stat"),
                                 digits.pval = gs("digits.pval"),
                                 smlab = NULL,
                                 backtransf = x$backtransf,
                                 lab.NA = gs("lab.NA"),
                                 equal.size = gs("equal.size"),
                                 ...) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  chkclass(x, "netcomparison")
  ##
  pooled <- setchar(pooled, c("common", "random", "fixed"))
  pooled[pooled == "fixed"] <- "common"
  ##
  nchar.comps <- replaceNULL(nchar.comps, 666)
  chknumeric(nchar.comps, min = 1, length = 1)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 0, length = 1)
  ##
  chklogical(backtransf)
  chkchar(lab.NA)
  
  
  ##
  ##
  ## (2) Extract results for common and random effects model
  ##
  ##
  if (pooled == "common") {
    TE   <- x$TE.common
    seTE <- x$seTE.common
    stat <- x$statistic.common
    pval <- x$pval.common
    ##
    if (is.null(smlab))
      smlab <- "Common Effects Model"
  }
  ##
  else if (pooled == "random") {
    TE   <- x$TE.random
    seTE <- x$seTE.random
    stat <- x$statistic.random
    pval <- x$pval.random
    ##
    if (is.null(smlab))
      smlab <- "Random Effects Model"
  }
  ##
  statistic <- formatN(stat, digits = digits.stat, text.NA = lab.NA)
  pval <- formatPT(pval, digits = digits.pval, lab.NA = lab.NA)
  ##
  ## Abbreviated component and treatment labels
  ##
  n.comparisons <- length(x$treat1)
  treat1 <- rep("", n.comparisons)
  treat2 <- rep("", n.comparisons)
  ##
  comps <- c(x$comps, x$inactive)
  comps.abbr <- treats(comps, nchar.comps)
  ##
  for (i in seq_len(n.comparisons))
    treat1[i] <- compos(x$treat1[i], comps, comps.abbr,
                        x$x$sep.comps, x$add1[i] == " ")
  for (i in seq_len(n.comparisons))
    treat2[i] <- compos(x$treat2[i], comps, comps.abbr,
                        x$x$sep.comps, x$add2[i] == " ")
  ##
  comps1.list <- compsplit(x$treat1, x$x$sep.comps)
  comps2.list <- compsplit(x$treat2, x$x$sep.comps)
  ##
  comps1.list <- lapply(comps1.list, charfac,
                        levels = comps, labels = comps.abbr)
  comps2.list <- lapply(comps2.list, charfac,
                        levels = comps, labels = comps.abbr)
  ##
  comparison <- rep("", n.comparisons)
  ##
  for (i in seq_len(n.comparisons)) {
    sel1.i <- !comps1.list[[i]] %in% comps2.list[[i]]
    sel2.i <- !comps2.list[[i]] %in% comps1.list[[i]]
    ##
    if (any(sel1.i) | any(sel2.i))
      comparison[i] <-
        paste0(paste(comps1.list[[i]][sel1.i],
                     collapse = paste0(x$add1[i], x$x$sep.comps, x$add1[i])),
               x$x$sep.trts,
               paste(comps2.list[[i]][sel2.i],
                     collapse = paste0(x$add2[i], x$x$sep.comps, x$add2[i])))
  }
  ##
  dat <- data.frame(comparison, treat1, treat2, TE, seTE, statistic, pval)
  
  
  ##
  ##
  ## (3) Generate forest plot
  ##
  ##
  m1 <-
    suppressWarnings(metagen(TE, seTE, data = dat, sm = x$x$sm,
                             studlab = comparison, backtransf = backtransf,
                             level = x$level,
                             method.tau = "DL", method.tau.ci = "",
                             warn = FALSE))
  ##
  forest(m1,
         digits = digits,
         common = FALSE, random = FALSE,
         overall = FALSE, hetstat = FALSE, test.subgroup = FALSE,
         leftcols = leftcols,
         leftlabs = leftlabs,
         rightcols = rightcols,
         rightlabs = rightlabs,
         smlab = smlab,
         lab.NA = lab.NA,
         weight.study = if (equal.size) "same" else pooled,
         just.addcols = "right",
         ...)
  
  
  invisible(NULL)
}


#' @rdname forest.netcomparison
#' @method plot netcomparison
#' @export

plot.netcomparison <- function(x, ...)
  forest(x, ...)
