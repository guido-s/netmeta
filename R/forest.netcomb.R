#' Forest plot for additive network meta-analysis
#' 
#' @description
#' Draws a forest plot in the active graphics window (using grid
#' graphics system).
#'
#' @aliases forest.netcomb plot.netcomb
#' 
#' @param x An object of class \code{netcomb}.
#' @param reference.group Reference treatment(s).
#' @param baseline.reference A logical indicating whether results
#'   should be expressed as comparisons of other treatments versus the
#'   reference treatment (default) or vice versa.
#' @param pooled A character string indicating whether results for the
#'   common (\code{"common"}) or random effects model (\code{"random"})
#'   should be plotted. Can be abbreviated.
#' @param equal.size A logical indicating whether all squares should
#'   be of equal size. Otherwise, the square size is proportional to
#'   the precision of estimates.
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
#' @param digits Minimal number of significant digits for treatment
#'   effects and confidence intervals, see \code{print.default}.
#' @param smlab A label printed at top of figure. By default, text
#'   indicating either common or random effects model is printed.
#' @param sortvar An optional vector used to sort the individual
#'   studies (must be of same length as the total number of
#'   treatments).
#' @param overall.hetstat A logical indicating whether to print heterogeneity
#'   measures.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in forest plots. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios, for example.
#' @param lab.NA A character string to label missing values.
#' @param add.data An optional data frame with additional columns to
#'   print in forest plot (see Details).
#' @param addrows.below.overall A numeric value indicating how many
#'   empty rows are printed between meta-analysis results and
#'   heterogeneity statistics.
#' @param drop.reference.group A logical indicating whether the
#'   reference group should be printed in the forest plot.
#' @param \dots Additional arguments for \code{\link[meta]{forest.meta}}
#'   function.
#' 
#' @details
#' A forest plot, also called confidence interval plot, is drawn in
#' the active graphics window.
#' 
#' Argument \code{sortvar} can be either a numeric or character vector
#' with length of number of treatments. If \code{sortvar} is numeric
#' the \code{\link[base]{order}} function is utilised internally to
#' determine the order of values. If \code{sortvar} is character it
#' must be a permutation of the treatment names. It is also possible
#' to to sort by treatment comparisons (\code{sortvar = TE}, etc.),
#' standard error (\code{sortvar = seTE}), and number of studies with
#' direct treatment comparisons (\code{sortvar = k}).
#' 
#' Argument \code{add.data} can be used to add additional columns to
#' the forest plot. This argument must be a data frame with the same
#' row names as the treatment effects matrices in R object \code{x},
#' i.e., \code{x$TE.common} or \code{x$TE.random}.
#' 
#' For more information see help page of \code{\link[meta]{forest.meta}}
#' function.
#'
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netcomb}}, \code{\link{discomb}},
#'   \code{\link[meta]{forest.meta}}
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
#' #
#' forest(nc1)
#' 
#' \dontrun{
#' # Specify, order of treatments
#' #
#' trts <- c("TCA", "SSRI", "SNRI", "NRI", "Low-dose SARI", "NaSSa",
#'   "rMAO-A", "Ind drug", "Hypericum", "Face-to-face CBT",
#'   "Face-to-face PST", "Face-to-face interpsy", "Face-to-face psychodyn",
#'   "Other face-to-face", "Remote CBT", "Self-help CBT", "No contact CBT",
#'   "Face-to-face CBT + SSRI", "Face-to-face interpsy + SSRI",
#'   "Face-to-face PST + SSRI", "UC", "Placebo")
#' #
#' # Note, three treatments are actually combinations of 'SSRI' with
#' # other components:
#' # "Face-to-face CBT + SSRI",
#' # "Face-to-face interpsy + SSRI",
#' # "Face-to-face PST + SSRI"
#' 
#' # Conduct random effects network meta-analysis
#' #
#' net2 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'   data = Linde2016, ref = "placebo",
#'   seq = trts, sm = "OR", common = FALSE)
#' #
#' net2
#' 
#' # Additive model for treatment components (with placebo as inactive
#' # treatment)
#' #
#' nc2 <- netcomb(net2, inactive = "placebo")
#' #
#' forest(nc2)
#' }
#' 
#' @method forest netcomb
#' @export

forest.netcomb <- function(x,
                           pooled = ifelse(x$random, "random", "common"),
                           reference.group = x$reference.group,
                           baseline.reference = x$baseline.reference,
                           equal.size = gs("equal.size"),
                           leftcols = "studlab",
                           leftlabs = "Treatment",
                           rightcols = c("effect", "ci"),
                           rightlabs = NULL,
                           digits = gs("digits.forest"),
                           smlab = NULL,
                           sortvar = x$seq,
                           overall.hetstat = gs("overall.hetstat"),
                           backtransf = x$backtransf,
                           lab.NA = gs("lab.NA"),
                           add.data,
                           addrows.below.overall =
                             if (x$overall.hetstat) 2 else
                               gs("addrows.below.overall"),
                           drop.reference.group = gs("drop.reference.group"),
                           ...) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  chkclass(x, "netcomb")
  x <- updateversion(x)
  ##
  is.discomb <- inherits(x, "discomb")
  ##
  pooled <- setchar(pooled, c("common", "random", "fixed"))
  pooled[pooled == "fixed"] <- "common"
  ##
  chknumeric(digits, min = 0, length = 1)
  ##
  chklogical(baseline.reference)
  chklogical(drop.reference.group)
  ##
  overall.hetstat <- replaceNULL(overall.hetstat, FALSE)
  #
  chklogical(backtransf)
  chkchar(lab.NA)
  
  
  ##
  ##
  ## (2) Extract results for common and random effects model
  ##     and calculate P-scores
  ##
  ##
  labels <- colnames(x$TE.common)
  ##
  if (reference.group == "") {
    warning("First treatment used as reference as ",
            "argument 'reference.group' is unspecified.")
    reference.group <- labels[1]
  }
  else
    reference.group <- setref(reference.group, labels)
  ##
  if (pooled == "common") {
    TE   <- x$TE.common
    seTE <- x$seTE.common
    ##
    text.common <- "(Common Effects Model)"
    ##
    if (is.null(smlab))
      if (baseline.reference)
        smlab <- paste0("Comparison: other vs '", reference.group, "'\n",
                        text.common)
      else
        smlab <- paste0("Comparison: '", reference.group, "' vs other\n",
                       text.common)
  }
  ##
  if (pooled == "random") {
    TE   <- x$TE.random
    seTE <- x$seTE.random
    if (is.null(smlab))
      if (baseline.reference)
        smlab <- paste0("Comparison: other vs '", reference.group,
                        "'\n(Random Effects Model)")
      else
        smlab <- paste0("Comparison: '", reference.group,
                        "' vs other\n(Random Effects Model)")
  }
  
  
  ##
  ##
  ## (3) Extract comparisons with reference group
  ##
  ##
  if (baseline.reference) {
    dat <- data.frame(TE = TE[, colnames(TE) == reference.group],
                      seTE = seTE[, colnames(seTE) == reference.group],
                      trts = colnames(TE),
                      k = NA,
                      row.names = colnames(TE),
                      as.is = TRUE)
    if (!is.discomb)
      dat$k <- x$A.matrix[, colnames(TE) == reference.group]
  }
  else {
    dat <- data.frame(TE = TE[rownames(TE) == reference.group, ],
                      seTE = seTE[rownames(seTE) == reference.group, ],
                      trts = rownames(TE),
                      k = NA,
                      row.names = colnames(TE),
                      as.is = TRUE)
    if (!is.discomb)
      dat$k <- x$A.matrix[rownames(TE) == reference.group, ]
  }
  ##
  rm(TE)
  rm(seTE)
  ##
  if (!missing(add.data)) {
    if (!is.data.frame(add.data))
      stop("Argument 'add.data' must be a data frame.",
           call. = FALSE)
    if (nrow(add.data) != length(labels))
      stop("Dataset 'add.data' must have ", nrow(dat),
           " rows (corresponding to number of treatments)",
           call. = FALSE)
    if (any(rownames(add.data) != labels))
      stop("Dataset 'add.data' must have the following row names:\n",
           paste(paste0("'", labels, "'"), collapse = " - "),
           call. = FALSE)
    ##
    dat <- cbind(dat, add.data)
  }
  
  
  ##
  ##
  ## (4) Sort dataset according to argument sortvar
  ##
  ##
  sortvar.c <- deparse(substitute(sortvar))
  sortvar.c <- gsub("\"", "", sortvar.c)
  ##
  idx5 <- charmatch(tolower(sortvar.c), "te", nomatch = NA)
  sel5 <- !is.na(idx5) & idx5 == 1
  if (any(sel5))
    sortvar <- dat$TE
  ##
  idx6 <- charmatch(tolower(sortvar.c), "-te", nomatch = NA)
  sel6 <- !is.na(idx6) & idx6 == 1
  if (any(sel6))
    sortvar <- -dat$TE
  ##
  idx7 <- charmatch(tolower(sortvar.c), "sete", nomatch = NA)
  sel7 <- !is.na(idx7) & idx7 == 1
  if (any(sel7))
    sortvar <- dat$seTE
  ##
  idx8 <- charmatch(tolower(sortvar.c), "-sete", nomatch = NA)
  sel8 <- !is.na(idx8) & idx8 == 1
  if (any(sel8))
    sortvar <- -dat$seTE
  ##
  if (!is.discomb) {
    idx9 <- charmatch(tolower(sortvar.c), "k", nomatch = NA)
    sel9 <- !is.na(idx9) & idx9 == 1
    if (any(sel9))
      sortvar <- dat$k
    ##
    idx10 <- charmatch(tolower(sortvar.c), "-k", nomatch = NA)
    sel10 <- !is.na(idx10) & idx10 == 1
    if (any(sel10))
      sortvar <- -dat$k
  }
  ##  
  if (!is.null(sortvar)) {
    if (is.character(sortvar))
      sort <- setseq(sortvar, labels)
    else
      sort <- order(sortvar)
    ##
    dat <- dat[sort, ]
  }
  
  
  ##
  ##
  ## (5) Generate forest plot
  ##
  ##
  if (drop.reference.group)
    dat <- subset(dat, trts != reference.group)
  ##
  trts <- dat$trts
  m1 <-
    suppressWarnings(metagen(TE, seTE, data = dat,
                             sm = x$sm,
                             studlab = trts, backtransf = backtransf,
                             method.tau = "DL", method.tau.ci = "",
                             warn = FALSE))
  ##
  forest(m1,
         digits = digits,
         #
         overall = FALSE, common = FALSE, random = FALSE,
         overall.hetstat = overall.hetstat,
         test.subgroup = FALSE,
         #
         leftcols = leftcols,
         leftlabs = leftlabs,
         rightcols = rightcols,
         rightlabs = rightlabs,
         #
         smlab = smlab, lab.NA = lab.NA,
         #
         weight.study = if (equal.size) "same" else pooled,
         #
         addrows.below.overall = addrows.below.overall,
         #
         ...)
  
  
  invisible(NULL)
}


#' @rdname forest.netcomb
#' @method plot netcomb
#' @export

plot.netcomb <- function(x, ...)
  forest(x, ...)
