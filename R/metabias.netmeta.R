#' Test of funnel plot asymmetry in network meta-analysis
#' 
#' @description
#' Test of funnel plot asymmetry in network meta-analysis
#' 
#' @param x An object of class \code{netmeta}.
#' @param order A mandatory character or numerical vector specifying
#'   the order of treatments or list of comparators (see Details).
#' @param pooled A character string indicating whether results for the
#'   common (\code{"common"}) or random effects model (\code{"random"})
#'   should be used in test of funnel plot asymmetry. Can be
#'   abbreviated.
#' @param method.bias A character vector indicating which test(s) for
#'   funnel plot asymmatrx to use. Admissible values are
#'   \code{"Begg"}, \code{"Egger"}, and \code{"Thompson"}, can be
#'   abbreviated. See function \code{\link[meta]{metabias.meta}}.
#' @param lump.comparator A logical indicating whether comparators
#'   should be lumped, e.g., to specify inactive treatments.
#'   information on direct comparisons should be added to the plot.
#' @param \dots Additional arguments (passed on to
#'   \code{\link[meta]{metabias.meta}}).
#' 
#' @details
#' Test of funnel plot asymmetry in network meta-analysis
#' 
#' Argument \code{order} is mandatory to determine the order of
#' treatments (Chaimani et al., 2013):
#' 
#' \emph{\dQuote{[...] investigators should order
#' the treatments in a meaningful way and make assumptions about how
#' small studies differ from large ones. For example, if they
#' anticipate that newer treatments are favored in small trials, then
#' they could name the treatments from oldest to newest so that all
#' comparisons refer to \sQuote{old versus new intervention}. Other
#' possibilities include defining the comparisons so that all refer to
#' an active treatment versus placebo or sponsored versus
#' non-sponsored intervention.}}
#' 
#' Alternatively, it is possible to only provide a single or few
#' treatment name(s) in argument \code{order} to define the
#' comparator(s). In this case only comparisons with this / these
#' treatment(s) will be considered. If argument \code{lump.comparator}
#' is \code{TRUE}, all comparators will be lumped into a single
#' group.
#' 
#' @return
#' A list with class \code{metabias} containing the following
#' components if a test for funnel plot asymmetry is conducted:
#' \item{statistic}{Test statistic.}
#' \item{df}{The degrees of freedom of the test statistic in
#'   the case that it follows a t distribution.}
#' \item{pval}{The p-value for the test.}
#' \item{estimate}{Estimates used to calculate test statisic.}
#' \item{method}{A character string indicating what type of test was
#'   used.}
#' \item{title}{Title of Cochrane review.}
#' \item{complab}{Comparison label.}
#' \item{outclab}{Outcome label.}
#' \item{var.model}{A character string indicating whether none,
#'   multiplicative, or additive residual heterogeneity variance was
#'   assumed.}
#' \item{method.bias}{As defined above.}
#' \item{x}{Network meta-analysis object.}
#' \item{version}{Version of R package \bold{meta} used to create
#'   object.}
#' \item{version.netmeta}{Version of R package \bold{netmeta} used to
#'   create object.}
#'
#' Or a list with the following elements if test is not conducted due
#' to the number of studies:
#' \item{k}{Number of comparisons.}
#' \item{k.min}{Minimum number of comparisons to perform test for
#'   funnel plot asymmetry.}
#' \item{version}{Version of R package \bold{meta} used to create
#'   object.}
#' \item{version.netmeta}{Version of R package \bold{netmeta} used to
#'   create object.}
#' 
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{funnel.netmeta}},
#'   \code{\link[meta]{metabias}}
#' 
#' @references
#' Chaimani A & Salanti G (2012):
#' Using network meta-analysis to evaluate the existence of
#' small-study effects in a network of interventions.
#' \emph{Research Synthesis Methods},
#' \bold{3}, 161--76
#' 
#' Chaimani A, Higgins JP, Mavridis D, Spyridonos P, Salanti G (2013):
#' Graphical tools for network meta-analysis in STATA.
#' PLOS ONE,
#' \bold{8}, e76654
#' 
#' @examples
#' \dontrun{
#' data(Senn2013)
#' 
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013, sm = "MD")
#' 
#' # Test for asymmetry in 'comparison-adjusted' funnel plot not
#' # conducted as argument 'order' is missing
#' #
#' try(metabias(net1))
#' 
#' # Test for funnel plot asymmetry comparing active treatments with
#' # placebo
#' metabias(net1, order = "pl")
#'
#' # Rank test
#' #
#' metabias(net1, order = "pl", method.bias = "Begg")
#'
#' 
#' # Test for funnel plot asymmetry based on (non-sensical) alphabetic
#' # order of treatments with placebo as last treatment
#' #
#' ord <- c("a", "b", "me", "mi", "pi", "r", "si", "su", "v", "pl")
#' metabias(net1, order = ord)
#' }
#'
#' @method metabias netmeta
#' @export

metabias.netmeta <- function(x,
                             order,
                             pooled = ifelse(x$random, "random", "common"),
                             method.bias = "Egger",
                             lump.comparator = gs("lump.comparator"),
                             ...) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  chkclass(x, "netmeta")
  x <- updateversion(x)
  sep.trts <- x$sep.trts
  ##
  pooled <- setchar(pooled, c("common", "random", "fixed"))
  pooled[pooled == "fixed"] <- "common"
  chklogical(lump.comparator)
  ##
  if (missing(order))
    stop("Argument 'order' with a meaningful order of treatments ",
         "must be provided.\n  ",
         "(see help page of funnel.netmeta for some examples).")
  else {
    if (length(order) != length(x$trts)) {
      order <- setchar(order, x$trts)
      order.all <- c(x$trts[!(x$trts %in% order)], order)
    }
    else
      order.all <- order <- setseq(order, x$trts)
  }
  ##
  if ((length(order) == length(order.all)) & lump.comparator) {
    warning("Argument 'lump.comparator' ignored as full treatment order is ",
            "provided.",
            call. = FALSE)
    lump.comparator <- FALSE
  }
  ##
  method.bias <- setmethodbias(method.bias, 1:3)
  
  
  ##
  ##
  ## (2) Get data
  ##
  ##
  TE <- x$TE
  seTE <- x$seTE
  studlab <- x$studlab
  ##
  treat1 <- x$treat1
  treat2 <- x$treat2
  ##
  comparison <- paste(treat1, treat2, sep = sep.trts)
  comparison21 <- paste(treat2, treat1, sep = sep.trts)
  ##
  treat1.pos <- as.numeric(factor(treat1, levels = order.all))
  treat2.pos <- as.numeric(factor(treat2, levels = order.all))
  ##
  wo <- treat1.pos > treat2.pos
  ##
  if (any(wo)) {
    ttreat1.pos <- treat1.pos
    treat1.pos[wo] <- treat2.pos[wo]
    treat2.pos[wo] <- ttreat1.pos[wo]
    ##
    TE[wo] <- -TE[wo]
    ##
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
    ##
    comparison[wo] <- comparison21[wo]
  }
  ##
  o <- order(treat1.pos, treat2.pos)
  ##
  TE <- TE[o]
  seTE <- seTE[o]
  ##
  studlab <- studlab[o]
  ##
  treat1 <- treat1[o]
  treat2 <- treat2[o]
  comparison <- comparison[o]
  ##
  res <- data.frame(studlab,
                    treat1, treat2, comparison,
                    TE, TE.direct = NA, TE.adj = NA, seTE)
  
  
  ##
  ##
  ## (3) Calculate 'comparison-adjusted' treatment effects
  ##
  ##
  if (is.numeric(treat1))
    treat1 <- as.character(treat1)
  if (is.numeric(treat2))
    treat2 <- as.character(treat2)
  ##
  if (pooled == "common")
    for (i in seq_along(res$TE))
      res$TE.direct[i] <- x$TE.direct.common[treat1[i], treat2[i]]
  else
    for (i in seq_along(res$TE))
      res$TE.direct[i] <- x$TE.direct.random[treat1[i], treat2[i]]
  ##
  res$TE.adj <- res$TE - res$TE.direct
  ##
  if (length(order) != length(order.all))
    res <- subset(res, treat1 %in% order | treat2 %in% order)
  ##
  if (lump.comparator) {
    res <- subset(res, !(treat1 %in% order & treat2 %in% order))
    ##
    res$treat2[res$treat2 %in% order] <- "...inactive..."
    ##
    res$comparison <- paste(res$treat1, res$treat2, sep = sep.trts)
  }
  
  
  ##
  ##
  ## (4) Conduct meta-analysis
  ##
  ##
  m.adj <-
    suppressWarnings(metagen(res$TE.adj, res$seTE,
                             studlab = res$studlab, sm = x$sm,
                             method.tau = "DL", method.tau.ci = ""))
  
  
  ##
  ##
  ## (5) Test of funnel plot asymmetry
  ##
  ##
  res <- metabias(m.adj, ...)
  ##
  res$version.netmeta <- packageDescription("netmeta")$Version
  
  
  res
}
