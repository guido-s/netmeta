#' Plot treatment ranking(s) of network meta-analyses
#' 
#' @description
#' Produce an image plot of treatment ranking(s) generated with R
#' function \code{netrank}.
#' 
#' @param ... A single netrank object or a list of netrank objects.
#' @param name An optional character vector providing descriptive
#'   names for the network meta-analysis objects.
#' @param common A logical indicating whether results for the common
#'   effects model should be plotted.
#' @param random A logical indicating whether results for the random
#'   effects model should be plotted.
#' @param seq A character or numerical vector specifying the sequence
#'   of treatments on the x-axis.
#' @param low A character string defining the colour for a P-score of
#'   0, see \code{\link[ggplot2]{scale_fill_gradient2}}.
#' @param mid A character string defining the colour for a P-score of
#'   0.5, see \code{\link[ggplot2]{scale_fill_gradient2}}.
#' @param high A character string defining the colour for a P-score of
#'   1, see \code{\link[ggplot2]{scale_fill_gradient2}}.
#' @param col Colour of text.
#' @param main Title.
#' @param main.size Font size of title, see
#'   \code{\link[ggplot2]{element_text}}.
#' @param main.col Colour of title, see
#'   \code{\link[ggplot2]{element_text}}.
#' @param main.face Font face of title, see
#'   \code{\link[ggplot2]{element_text}}.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param axis.size Font size of axis text, see
#'   \code{\link[ggplot2]{element_text}}.
#' @param axis.col Colour of axis text, see
#'   \code{\link[ggplot2]{element_text}}.
#' @param axis.face Font face of axis text, see
#'   \code{\link[ggplot2]{element_text}}.
#' @param na.value Colour for missing values, see
#'   \code{\link[ggplot2]{scale_fill_gradient2}}.
#' @param angle Angle for text on x-axis, see
#'   \code{\link[ggplot2]{element_text}}.
#' @param hjust.x A numeric between 0 and 1 with horizontal
#'   justification of text on x-axis, see
#'   \code{\link[ggplot2]{element_text}}.
#' @param vjust.x A numeric between 0 and 1 with vertical
#'   justification of text on x-axis, see
#'   \code{\link[ggplot2]{element_text}}.
#' @param hjust.y A numeric between 0 and 1 with horizontal
#'   justification of text on y-axis, see
#'   \code{\link[ggplot2]{element_text}}.
#' @param vjust.y A numeric between 0 and 1 with vertical
#'   justification of text on y-axis, see
#'   \code{\link[ggplot2]{element_text}}.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param fixed Deprecated argument (replaced by 'common').
#' @param comb.fixed Deprecated argument (replaced by 'common').
#' @param comb.random Deprecated argument (replaced by 'random').
#'
#' @details
#' This function produces an image plot of network rankings (Palpacuer
#' et al., 2018, Figure 4). Note, a scatter plot of two network
#' rankings can be generated with \code{\link{plot.netposet}}.
#'
#' By default, treatments are ordered by decreasing P-scores of the
#' first network meta-analysis object. Argument \code{seq} can be used
#' to specify a differenct treatment order.
#' 
#' @return
#' A ggplot2 object or NULL if no ranking was conducted.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}, Clément
#'   Palpacuer \email{clementpalpacuer@@gmail.com}
#' 
#' @seealso \code{\link{netrank}}, \code{\link{netmeta}},
#'   \code{\link{netposet}}, \code{\link{hasse}}
#' 
#' @references
#' Palpacuer C, Duprez R, Huneau A, Locher C, Boussageon R, Laviolle
#' B, et al. (2018):
#' Pharmacologically controlled drinking in the treatment of alcohol
#' dependence or alcohol use disorders: a systematic review with
#' direct and network meta-analyses on nalmefene, naltrexone,
#' acamprosate, baclofen and topiramate.
#' \emph{Addiction},
#' \bold{113}, 220--37
#' 
#' @keywords hplot
#' 
#' @examples
#' \dontrun{
#' # Use depression dataset
#' #
#' data(Linde2015)
#' 
#' # Define order of treatments
#' #
#' trts <- c("TCA", "SSRI", "SNRI", "NRI",
#'   "Low-dose SARI", "NaSSa", "rMAO-A", "Hypericum", "Placebo")
#'
#' # Outcome labels
#' #
#' outcomes <- c("Early response", "Early remission")
#' 
#' # (1) Early response
#' #
#' p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#' #
#' net1 <- netmeta(p1, common = FALSE,
#'   seq = trts, ref = "Placebo")
#' 
#' # (2) Early remission
#' #
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#' #
#' net2 <- netmeta(p2, common = FALSE,
#'   seq = trts, ref = "Placebo")
#' 
#' # Image plot of treatment rankings (two outcomes)
#' #
#' plot(netrank(net1, small.values = "bad"),
#'   netrank(net2, small.values = "bad"),
#'   name = outcomes, digits = 2)
#'
#' 
#' # Outcome labels
#' #
#' outcomes <- c("Early response", "Early remission",
#'   "Lost to follow-up", "Lost to follow-up due to AEs",
#'   "Adverse events (AEs)")
#' 
#' # (3) Loss to follow-up
#' #
#' p3 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss1, loss2, loss3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#' #
#' net3 <- netmeta(p3, common = FALSE, seq = trts, ref = "Placebo")
#' 
#' # (4) Loss to follow-up due to adverse events
#' #
#' p4 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss.ae1, loss.ae2, loss.ae3), n = list(n1, n2, n3),
#'   studlab = id, data = subset(Linde2015, id != 55), sm = "OR")
#' #
#' net4 <- netmeta(p4, common = FALSE, seq = trts, ref = "Placebo")
#' 
#' # (5) Adverse events
#' #
#' p5 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(ae1, ae2, ae3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#' #
#' net5 <- netmeta(p5, common = FALSE, seq = trts, ref = "Placebo")
#' 
#' # Image plot of treatment rankings (two outcomes)
#' #
#' plot(netrank(net1, small.values = "bad"),
#'   netrank(net2, small.values = "bad"),
#'   netrank(net3, small.values = "good"),
#'   netrank(net4, small.values = "good"),
#'   netrank(net5, small.values = "good"),
#'   name = outcomes, digits = 2)
#' }
#' 
#' @method plot netrank
#' @export


plot.netrank <- function(...,
                         name,
                         common, random,
                         ##
                         seq,
                         ##
                         low = "red", mid = "yellow", high = "green",
                         ##
                         col = "black",
                         ##
                         main,
                         main.size = 14, main.col = col,
                         main.face = "bold",
                         ##
                         legend = TRUE,
                         ##
                         axis.size = 12, axis.col = col,
                         axis.face = "plain",
                         ##
                         na.value = "grey50",
                         ##
                         angle = 45,
                         ##
                         hjust.x = 1, vjust.x = 1,
                         hjust.y = 1, vjust.y = 0,
                         ##
                         nchar.trts,
                         ##
                         digits = 3,
                         ##
                         fixed,
                         comb.fixed, comb.random,
                         ##
                         warn.deprecated = gs("warn.deprecated")) {
  
  ##
  ##
  ## (1) Identify list elements of class netrank
  ##
  ##
  args <- list(...)
  names.args <- names(args)
  ##
  n.args <- length(args)
  n.i <- seq_len(n.args)
  ##
  is.netrank <- rep_len(FALSE, n.args)
  ##
  for (i in n.i)
    is.netrank[i] <- inherits(args[[i]], "netrank")
  ##
  n.netrank <- sum(is.netrank)
  ##
  if (n.netrank == 0)
    stop("At least one R object of class 'netrank' must be provided.",
         call. = FALSE)
  ##
  for (i in n.i)
    if (is.netrank[i])
      args[[i]] <- updateversion(args[[i]])
  ##
  missing.name <- missing(name)
  missing.common <- missing(common)
  missing.random <- missing(random)
  missing.seq <- missing(seq)
  missing.main <- missing(main)
  missing.nchar.trts <- missing(nchar.trts)
  missing.fixed <- missing(fixed)
  missing.comb.fixed <- missing(comb.fixed)
  missing.comb.random <- missing(comb.random)
  
  
  ##
  ##
  ## (2) Identify additional arguments
  ##
  ##
  formal.args <- c("name", "common", "random",
                   "seq", "low", "mid", "high", "col",
                   "main", "main.size", "main.col", "main.face",
                   "legend", "axis.size", "axis.col", "axis.face",
                   "na.value", "angle",
                   "hjust.x", "vjust.x", "hjust.y", "vjust.y",
                   "nchar.trts", "digits",
                   "fixed", "comb.fixed", "comb.random", "warn.deprecated")
  ##
  for (i in n.i) {
    if (!is.netrank[i]) {
      cm <- charmatch(names.args[i], formal.args)
      ##
      if (is.na(cm))
        stop("Argument '", names.args[i],
             "' must be an R object of class 'netrank'.",
             call. = FALSE)
      else if (cm == 0)
        stop("Argument '", names.args[i],
             "' matches multiple formal arguments.",
             call. = FALSE)
      else {
        if (cm == 1) {
          missing.name <- FALSE
          name <- args[[i]]
        }
        else if (cm == 2) {
          missing.common <- FALSE
          common <- args[[i]]
        }
        else if (cm == 3) {
          missing.random <- FALSE
          random <- args[[i]]
        }
        else if (cm == 4) {
          missing.seq <- FALSE
          seq <- args[[i]]
        }
        else if (cm == 5)
          low <- args[[i]]
        else if (cm == 6)
          mid <- args[[i]]
        else if (cm == 7)
          high <- args[[i]]
        else if (cm == 8)
          col <- args[[i]]
        else if (cm == 9) {
          missing.main <- FALSE
          main <- args[[i]]
        }
        else if (cm == 10)
          main.size <- args[[i]]
        else if (cm == 11)
          main.col <- args[[i]]
        else if (cm == 12)
          main.face <- args[[i]]
        else if (cm == 13)
          legend <- args[[i]]
        else if (cm == 14)
          axis.size <- args[[i]]
        else if (cm == 15)
          axis.col <- args[[i]]
        else if (cm == 16)
          axis.face <- args[[i]]
        else if (cm == 17)
          na.value <- args[[i]]
        else if (cm == 18)
          angle <- args[[i]]
        else if (cm == 19)
          hjust.x <- args[[i]]
        else if (cm == 20)
          vjust.x <- args[[i]]
        else if (cm == 21)
          hjust.y <- args[[i]]
        else if (cm == 22)
          vjust.y <- args[[i]]
        else if (cm == 23) {
          missing.nchar.trts <- FALSE
          nchar.trts <- args[[i]]
        }
        else if (cm == 24)
          digits <- args[[i]]
        else if (cm == 25) {
          missing.comb.common <- FALSE
          comb.common <- args[[i]]
        }
        else if (cm == 26) {
          missing.comb.random <- FALSE
          comb.random <- args[[i]]
        }
        else if (cm == 27)
          warn.deprecated <- args[[i]]
      }
    }
  }
  
  
  ##
  ##
  ## (3) Check arguments
  ##
  ##
  chkchar(low)
  chkchar(mid)
  chkchar(high)
  ##
  chklogical(legend)
  chknumeric(main.size, min = 0, length = 1)
  chkchar(main.face)
  chknumeric(axis.size, min = 0, length = 1)
  chkchar(axis.face)
  chknumeric(angle, min = -360, max = 360, length = 1)
  chknumeric(hjust.x, min = 0, max = 1, length = 1)
  chknumeric(vjust.x, min = 0, max = 1, length = 1)
  chknumeric(hjust.y, min = 0, max = 1, length = 1)
  chknumeric(vjust.y, min = 0, max = 1, length = 1)
  chknumeric(digits, min = 0, length = 1)
  ##
  print.warning1 <- FALSE
  print.warning2 <- FALSE
  print.warning3 <- FALSE
  print.warning4 <- FALSE
  print.warning5 <- FALSE
  
  
  ##
  ##
  ## (4) Set names of network meta-analysis objects
  ##
  ##
  if (missing.name)
    name <- paste0("netmeta", seq_len(n.netrank))
  else {
    if (length(name) != n.netrank) 
      stop("Number of network meta-analyses and ",
           "names provided in argument 'name' differ.",
           call. = FALSE)
    ##
    if (length(unique(name)) != length(name)) {
      warning1 <-
        paste0("Network meta-analyses are labelled 'netmeta1' to 'netmeta",
               n.netrank,
               "' as values of argument 'name' are not all disparate.")
      print.warning1 <- TRUE
      ##
      name <- paste0("netmeta", seq_len(n.netrank))
    }
  }
  
  
  ##
  ##
  ## (5) Determine common
  ##
  ##
  if (missing.common & (!missing.fixed | !missing.comb.fixed)) {
    chklogical(warn.deprecated)
    common <-
      deprecated2(common, missing.common, comb.fixed, missing.comb.fixed,
                  warn.deprecated)
    common <-
      deprecated2(common, missing.common, fixed, missing.fixed,
                  warn.deprecated)
    missing.common <- FALSE
  }
  ##
  if (missing.common & missing.fixed & missing.comb.fixed) {
    cfs <- logical(0)
    ##
    for (i in n.i)
      if (is.netrank[i])
        cfs[i] <- args[[i]]$x$common
      else
        cfs[i] <- FALSE
    ##
    cfs <- unique(cfs[is.netrank])
    ##
    if (length(cfs) != 1) {
      common <- TRUE
      warning2 <- paste0("Argument 'common' set to TRUE ",
                         "(as it is not unique in network meta-analyses).")
      print.warning2 <- TRUE
    }
    else
      common <- cfs
  }
  else {
    chklogical(warn.deprecated)
    common <-
      deprecated2(common, missing.common, comb.fixed, missing.comb.fixed,
                  warn.deprecated)
    common <-
      deprecated2(common, missing.common, fixed, missing.fixed,
                  warn.deprecated)
    missing.common <- FALSE
  }
  ##
  chklogical(common)
  
  
  ##
  ##
  ## (6) Determine random
  ##
  ##
  if (missing.random & missing.comb.random) {
    crs <- logical(0)
    ##
    for (i in n.i)
      if (is.netrank[i])
        crs[i] <- args[[i]]$x$random
      else
        crs[i] <- FALSE
    ##
    crs <- unique(crs[is.netrank])
    ##
    if (length(crs) != 1) {
      random <- TRUE
      warning3 <- paste0("Argument 'random' set to TRUE ",
                         "(as it is not unique in network meta-analyses).")
      print.warning3 <- TRUE
    }
    else
      random <- crs
  }
  else {
    chklogical(warn.deprecated)
    random <-
      deprecated2(random, missing.random, comb.random, missing.comb.random,
                  warn.deprecated)
    missing.random <- FALSE
  }
  ##
  chklogical(random)
  ##
  if (common & random) {
    warning4 <- paste0("P-scores for random effects model displayed ",
                       "as both common and random effects ",
                       "network meta-analysis was conducted.")
    print.warning4 <- TRUE
    common <- FALSE
  }
  else if (!common & !random) {
    warning("No plot generated as ranking was neither for common nor ",
            "random effects network meta-analysis conducted.")
    return(invisible(NULL))
  }
  
  
  ##
  ##
  ## (7) Determine nchar.trts
  ##
  ##
  if (missing.nchar.trts) {
    cns <- logical(0)
    ##
    for (i in n.i)
      if (is.netrank[i])
        cns[i] <- args[[i]]$x$nchar.trts
      else
        cns[i] <- 666
    ##
    cns <- unique(cns[is.netrank])
    ##
    if (length(cns) != 1) {
      nchar.trts <- 666
      warning5 <- paste0("Argument 'nchar.trts' set to 666 ",
                         "(as it is not unique in network meta-analyses).")
      print.warning5 <- TRUE
    }
    else
      nchar.trts <- cns
  }
  else
    chknumeric(nchar.trts, min = 1, length = 1)
  
  
  ##
  ##
  ## (8) Determine list of all treatments
  ##
  ##
  trts <- character(0)
  ##
  for (i in n.i)
    if (is.netrank[i])
      trts <- c(trts, names(args[[i]]$ranking.common))
  ##
  trts <- unique(trts)
  ##  
  if (!missing.seq)
    seq <- setseq(seq, trts)
  else {
    trts1 <- data.frame(treat = trts, ranking = NA,
                        row.names = trts,
                        stringsAsFactors = FALSE)
    first <- min(seq_along(is.netrank)[is.netrank])
    trts.first <- names(args[[first]]$ranking.common)
    ##
    if (random)
      trts1[trts.first, "ranking"] <- args[[first]]$ranking.random
    else
      trts1[trts.first, "ranking"] <- args[[first]]$ranking.common
    ##
    trts1 <- trts1[rev(order(trts1$ranking, na.last = FALSE)), ]
    seq <- trts1$treat
  }
  ##
  n.trts <- length(trts)
  
  
  ##
  ##
  ## (9) Generate dataset for plot
  ##
  ##
  dat <- data.frame(name = character(0), treat = character(0),
                    ranking = numeric(0),
                    stringsAsFactors = FALSE)
  ##  
  for (i in n.i) {
    if (is.netrank[i]) {
      dat.i <- data.frame(name = rep_len(name[i], n.trts),
                          treat = trts, ranking = NA,
                          row.names = trts,
                          stringsAsFactors = FALSE)
      ##
      trts.i <- names(args[[i]]$ranking.common)
      ##
      if (random)
        dat.i[trts.i, "ranking"] <- args[[i]]$ranking.random
      else 
        dat.i[trts.i, "ranking"] <- args[[i]]$ranking.common
      ##
      dat <- rbind(dat, dat.i)
    }
  }
  ##
  row.names(dat) <- seq_len(nrow(dat))
  ##
  dat$ranking <- as.character(formatN(round(dat$ranking, digits),
                                      digits = digits, text.NA = ""))
  ##
  dat$name <- factor(dat$name, levels = rev(name))
  ##
  dat$treat <- factor(dat$treat, levels = seq,
                      labels = treats(seq, nchar.trts))
  ## Get rid of warnings 'Undefined global functions or variables'
  treat <- dat$treat
  ranking <- dat$ranking
  
  
  ##
  ##
  ## (10) Create ggplot2 object
  ##
  ##
  plt <- ggplot(dat,
                aes(x = treat, y = name, fill = as.numeric(ranking)))
  
  plt <- plt + theme_classic()
  
  plt <- plt + geom_tile() + ggplot2::xlab("") + ggplot2::ylab("") +
    theme(line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(size = axis.size, face = axis.face,
                                     colour = axis.col,
                                     angle = angle,
                                     hjust = hjust.x, vjust = vjust.x),
          axis.text.y = element_text(size = axis.size, face = axis.face,
                                     colour = axis.col,
                                     angle = 0,
                                     hjust = hjust.y, vjust = vjust.y),
          legend.title = element_text(colour = col,
                                      size = 7, face = "bold"),
          legend.text = element_text(colour = col, size = 10)) +
    scale_fill_gradient2(low = low, mid = mid, high = high,
                         midpoint = 0.5, limit = c(0, 1),
                         space = "Lab", name = "P-scores",
                         na.value = na.value)
  
  plt <- plt + geom_text(aes(x = treat, y = name, label = ranking),
                         colour = col)
  
  if (!legend)
    plt <- plt + theme(legend.position = "none")
  
  if (!missing.main)
    plt <- plt + ggtitle(main) +
      theme(plot.title = element_text(size = main.size, face = main.face,
                                      colour = main.col))
  
  
  ##
  ## Print warnings
  ##
  if (print.warning1)
    warning(warning1, call. = FALSE)
  if (print.warning2)
    warning(warning2, call. = FALSE)
  if (print.warning3)
    warning(warning3, call. = FALSE)
  if (print.warning4)
    warning(warning4, call. = FALSE)
  if (print.warning5)
    warning(warning5, call. = FALSE)
  
  
  plt
}
