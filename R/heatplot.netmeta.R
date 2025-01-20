#' Heat Plot
#'
#' @description
#' Produces a heat plot containing treatment estimates with confidence
#' intervals for all possible pairwise comparisons.
#'
#' @param x An object of class \code{netmeta}.
#' @param pooled A character string indicating whether results for the
#'   common (\code{"common"}) or random effects model
#'   (\code{"random"}) should be plotted. Can be abbreviated.
#' @param seq A character or numerical vector specifying the sequence
#'   of treatments in rows and columns of the heat plot.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param low.colour A string indicating the colour of low relative
#'   treatment effects for the heat plot (e.g odds ratio of ~0.5)
#' @param mid.colour A string indicating the colour of null relative
#'   treatment effects for the heat plot (e.g odds ratio of ~1.0).
#' @param high.colour A string indicating the colour of high relative
#'   treatment effects for the heat plot (e.g odds ratio of ~2.0).
#' @param size The size of cell entries with the relative treatment
#'   effect and confidence intervals.
#' @param size.trt The size of treatment names placed on the top and
#'   left of the plot.
#' @param size.axis The size of labels on the top and left of the plot
#' @param digits Minimal number of significant digits for treatment
#'   effects and confidence intervals, see \code{print.default}.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in forest plots. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios, for example.
#' @param \dots Additional arguments.
#'
#' @return
#' League heat plot, where a color scale is used to represent the
#' values of relative treatment effects.
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Senn2013)
#' 
#' # Only consider first five studies (to reduce runtime of example)
#' #
#' studies <- unique(Senn2013$studlab)
#' Senn2013.5 <- subset(Senn2013, studlab %in% studies[1:5])
#' 
#' # Conduct random effects network meta-analysis with
#' # placebo as reference treatment
#' #
#' net1 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'   data = Senn2013.5, sm = "MD", common = FALSE, reference = "plac")
#'       
#' # Generate a heat plot (with abbreviated treatment labels)
#' #
#' heatplot(net1, nchar.trts = 4) 
#'
#' @method heatplot netmeta
#' @export

heatplot.netmeta <- function(x,
                             pooled = ifelse(x$random, "random", "common"),
                             seq = x$seq,
                             nchar.trts = x$nchar.trts,
                             ##
                             low.colour = "red",
                             mid.colour = "white",
                             high.colour = "springgreen4",
                             ##
                             size = 6,
                             size.trt = 16,
                             size.axis = 12,
                             ##
                             digits = gs("digits.forest"),
                             backtransf = x$backtransf,
                             ...) {
  
  chkclass(x, "netmeta")
  x <- updateversion(x)
  sm <- x$sm
  ##
  pooled <- setchar(pooled, c("common", "random", "fixed"))
  pooled[pooled == "fixed"] <- "common"
  ##
  if (is.null(x$nchar.trts))
    nchar.trts <- 666
  chknumeric(nchar.trts, min = 1, length = 1)
  ##
  chknumeric(digits, min = 0, length = 1)
  ##
  if (is_untransformed(sm))
    backtransf <- TRUE
  backtransf <- replaceNULL(backtransf, TRUE)
  chklogical(backtransf)
  
  
  trts <- x$trts
  ##
  if (!missing(seq) & is.null(seq))
    stop("Argument 'seq' must be not NULL.")
  ##
  if (is.null(seq) | (length(seq) == 1 & x$d == 1))
    seq1 <- 1:length(trts)
  else
    seq1 <- charmatch(setseq(seq, x$trts), x$trts)
  
  
  if (pooled == "common") {
    TE.nma <- x$TE.common[seq1, seq1]
    lower.nma <- x$lower.common[seq1, seq1]
    upper.nma <- x$upper.common[seq1, seq1]
  }
  else {
    TE.nma <- x$TE.random[seq1, seq1]
    lower.nma <- x$lower.random[seq1, seq1]
    upper.nma <- x$upper.random[seq1, seq1]
  }
  ##
  noeffect <- 1L * (backtransf & is_relative_effect(sm))
  #
  if (backtransf) {
    TE.nma    <- backtransf(TE.nma, sm)
    lower.nma <- backtransf(lower.nma, sm)
    upper.nma <- backtransf(upper.nma, sm)
    #
    # Switch lower and upper limit for VE if results have been
    # backtransformed
    #
    if (sm == "VE") {
      tmp.l <- lower.nma
      lower.nma <- upper.nma
      upper.nma <- tmp.l
    }
  }
  ##
  TE.nma    <- round(TE.nma, digits)
  lower.nma <- round(lower.nma, digits)
  upper.nma <- round(upper.nma, digits)
  ##
  trts <- trts[seq1]
  ct <- heattrts(TE.nma)
  ##
  hdata <- data.frame(Treatment = c(ct$treat2, ct$treat1),
                      Comparator = c(ct$treat1, ct$treat2))
  ##
  hdata$TE.nma <- c(lowertri(TE.nma), uppertri(TE.nma))
  hdata$lower.nma <- c(lowertri(lower.nma), uppertri(lower.nma))
  hdata$upper.nma <- c(lowertri(upper.nma), uppertri(upper.nma))
  ##
  trts.abbr <- treats(trts, nchar.trts)
  ##
  hdata$Treatment <-
    as.character(factor(hdata$Treatment,
                        levels = trts, labels = trts.abbr))
  ##
  hdata$Comparator <-
    as.character(factor(hdata$Comparator,
                        levels = trts, labels = trts.abbr))
  
  
  ##
  ##
  ## Create heat plot
  ##
  ##
  Treatment <- Comparator <- TE.nma <- lower.nma <- upper.nma <- NULL
  ##
  hplot <-
    ggplot(data = hdata,
           aes(x = Treatment, y = Comparator, fill = TE.nma)) +
    geom_tile() +
    geom_text(aes(label = paste0(formatN(TE.nma, digits), "\n",
                                 formatCI(round(lower.nma, digits),
                                          round(upper.nma, digits)))),
              size = size)
  ##
  hplot <- hplot +
    scale_fill_gradient2(low = low.colour,
                         mid = mid.colour,
                         high = high.colour,
                         midpoint = noeffect) +
    theme_dark() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none", panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(size = size.trt),
          axis.title = element_text(size = size.axis)) +
    scale_x_discrete(limits = trts.abbr, expand = c(0, 0),
                     position = "top") +
    scale_y_discrete(limits = rev(trts.abbr), expand = c(0, 0))
  
  return(hplot)
}


heattrts <- function(x) {
  trts <- rownames(x)
  n <- length(trts)
  ##
  treat1 <- treat2 <- vector("character", 0)
  k <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (i != j) {
        k <- k + 1
        treat1[k] <- trts[j]
        treat2[k] <- trts[i]
      }
    }
  }
  ##
  data.frame(treat1, treat2)
}
