#' Plot network ranking(s)
#' 
#' @description
#' Plot network ranking(s)
#' 
#' @param ... A single netrank object or a list of netrank objects.
#' @param name name
#' @param comb.fixed comb.fixed
#' @param comb.random comb.random
#' @param seq seq
#' @param digits digits
#' 
#' @return
#' A ggplot2 object:
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netrank}}, \code{\link{netmeta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' # Use depression dataset
#' #
#' data(Linde2015)
#' 
#' # Define order of treatments
#' #
#' trts <- c("TCA", "SSRI", "SNRI", "NRI",
#'           "Low-dose SARI", "NaSSa", "rMAO-A", "Hypericum",
#'           "Placebo")
#'
#' # Outcome labels
#' #
#' outcomes <- c("Early response", "Early remission")
#' 
#' # (1) Early response
#' #
#' p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'                event = list(resp1, resp2, resp3),
#'                n = list(n1, n2, n3),
#'                studlab = id, data = Linde2015, sm = "OR")
#' #
#' net1 <- netmeta(p1, comb.fixed = FALSE,
#'                 seq = trts, ref = "Placebo")
#' 
#' # (2) Early remission
#' #
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'                event = list(remi1, remi2, remi3),
#'                n = list(n1, n2, n3),
#'                studlab = id, data = Linde2015, sm = "OR")
#' #
#' net2 <- netmeta(p2, comb.fixed = FALSE,
#'                 seq = trts, ref = "Placebo")
#' 
#' # Image plot of treatment rankings (two outcomes)
#' #
#' plot(netrank(net1, small.values = "bad"),
#'      netrank(net2, small.values = "bad"),
#'      name = outcomes, digits = 2)
#'
#' \dontrun{
#' #
#' # Outcome labels
#' #
#' outcomes <- c("Early response", "Early remission",
#'               "Lost to follow-up", "Lost to follow-up due to AEs",
#'               "Adverse events (AEs)")
#' 
#' # (3) Loss to follow-up
#' #
#' p3 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'                event = list(loss1, loss2, loss3),
#'                n = list(n1, n2, n3),
#'                studlab = id, data = Linde2015, sm = "OR")
#' #
#' net3 <- netmeta(p3, comb.fixed = FALSE,
#'                 seq = trts, ref = "Placebo")
#' 
#' # (4) Loss to follow-up due to adverse events
#' #
#' p4 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'                event = list(loss.ae1, loss.ae2, loss.ae3),
#'                n = list(n1, n2, n3),
#'                studlab = id, data = subset(Linde2015, id != 55),
#'                sm = "OR")
#' #
#' net4 <- netmeta(p4, comb.fixed = FALSE,
#'                 seq = trts, ref = "Placebo")
#' 
#' # (5) Adverse events
#' #
#' p5 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'                event = list(ae1, ae2, ae3),
#'                n = list(n1, n2, n3),
#'                studlab = id, data = Linde2015, sm = "OR")
#' #
#' net5 <- netmeta(p5, comb.fixed = FALSE,
#'                 seq = trts, ref = "Placebo")
#' 
#' # Image plot of treatment rankings (two outcomes)
#' #
#' plot(netrank(net1, small.values = "bad"),
#'      netrank(net2, small.values = "bad"),
#'      netrank(net3, small.values = "good"),
#'      netrank(net4, small.values = "good"),
#'      netrank(net5, small.values = "good"),
#'      name = outcomes, digits = 2)
#' }
#' 
#' @method plot netrank
#' @export
#' @export plot.netrank
#'
#' @importFrom ggplot2 ggplot aes theme_classic geom_tile xlab ylab
#'   theme element_blank element_text scale_fill_gradient2 geom_text


plot.netrank <- function(..., name,
                         comb.fixed, comb.random,
                         ##
                         seq,
                         digits = gs("digits.prop")) {
  
  
  ##
  ##
  ## (1) Extract list elements and basic checks
  ##
  ##
  chklogical <- meta:::chklogical
  ##
  args <- list(...)
  ##
  n.netmeta <- length(args)
  n.i <- seq_len(n.netmeta)
  ##  
  for (i in n.i) {
    if (!inherits(args[[i]], "netrank"))
      stop("All elements of argument '...' must be of class 'netrank'.",
           call. = FALSE)
  }


  ##
  ##
  ## (2) Set names of network meta-analysis objects
  ##
  ##
  if (missing(name))
    name <- paste0("netmeta", n.i)
  else {
    if (length(name) != n.netmeta) 
      stop("Number of network meta-analyses and ",
           "names provided in argument 'name' differ.",
           call. = FALSE)
    ##
    if (length(unique(name)) != length(name)) {
      warning1 <- paste0("Network meta-analyses are labelled 'netmeta1' to 'netmeta",
                         n.netmeta,
                         "' as values of argument 'name' are not all disparate.")
      print.warning1 <- TRUE
      ##
      name <- paste0("netmeta", n.i)
    }
  }
  
  
  ##
  ## (3) Determine comb.fixed
  ##
  if (missing(comb.fixed)) {
    cfs <- logical(0)
    ##
    for (i in n.i)
      cfs[i] <- args[[i]]$x$comb.fixed
    ##
    cfs <- unique(cfs)
    ##
    if (length(cfs) != 1) {
      comb.fixed <- TRUE
      warning2 <- paste0("Argument 'comb.fixed' set to TRUE ",
                         "(as it is not unique in network meta-analyses).")
      print.warning2 <- TRUE
    }
    else
      comb.fixed <- cfs
  }
  else
    comb.fixed <- chklogical(comb.fixed)
  
  
  ##
  ## (4) Determine comb.random
  ##
  if (missing(comb.random)) {
    crs <- logical(0)
    ##
    for (i in n.i)
      crs[i] <- args[[i]]$x$comb.random
    ##
    crs <- unique(crs)
    ##
    if (length(crs) != 1) {
      comb.random <- TRUE
      warning3 <- paste0("Argument 'comb.random' set to TRUE ",
                         "(as it is not unique in network meta-analyses).")
      print.warning3 <- TRUE
    }
    else
      comb.random <- crs
  }
  else
    comb.random <- chklogical(comb.random)
  

  if (comb.fixed & comb.random) {
    warning("P-scores for random effects model displayed ",
            "as both fixed and random effects network meta-analysis ",
            "was conducted.")
    comb.fixed <- FALSE
  }
  else if (!comb.fixed & !comb.random) {
    warning("No plot generated as either fixed nor ",
            "random effects network meta-analysis was conducted.")
    return(invisible(NULL))
  }
  
  
  trts <- character(0)
  ##
  for (i in n.i)
    trts <- c(trts, names(args[[i]]$Pscore.fixed))
  ##
  trts <- unique(trts)

  
  if (!missing(seq))
    seq <- setseq(seq, trts)
  else
    seq <- trts
  ##
  n.trts <- length(trts)
  
  
  dat <- data.frame(name = character(0), treat = character(0),
                    pscore = numeric(0),
                    stringsAsFactors = FALSE)
  ##  
  for (i in n.i) {
    dat.i <- data.frame(name = rep_len(name[i], n.trts),
                        treat = trts, pscore = NA,
                        row.names = trts,
                        stringsAsFactors = FALSE)
    ##
    trts.i <- names(args[[i]]$Pscore.fixed)
    ##
    if (comb.random)
      dat.i[trts.i, "pscore"] <- args[[i]]$Pscore.random
    else 
      dat.i[trts.i, "pscore"] <- args[[i]]$Pscore.fixed
    ##
    dat <- rbind(dat, dat.i)
  }
  ##
  row.names(dat) <- seq_len(nrow(dat))
  ##
  dat$pscore <- as.character(meta:::formatN(round(dat$pscore, digits),
                                            digits = digits, text.NA = ""))
  ##
  dat$name <- factor(dat$name, levels = rev(name))
  dat$treat <- factor(dat$treat, levels = seq)
  
  
  plt <- ggplot(dat,
                aes(x = dat$treat, y = dat$name,
                    fill = as.numeric(as.character(dat$pscore))))
  ##
  ## OPTIONS
  ##
  plt <- plt + theme_classic()
  
  plt <- plt + geom_tile() + xlab("") + ylab("") +
    theme(line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(colour = "grey20", size = 12, angle = 80,
                                     hjust = .5,vjust = .5,face = "plain"),
          axis.text.y= element_text(colour = "grey20",size = 12,angle = 0,
                                    hjust = 1,vjust = 0,face = "plain"),
          legend.title = element_text(size = 7, face = "bold"),
          legend.text = element_text(size = 10)) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                         midpoint = 0.5, limit = c(0,1),
                         space = "Lab", name = "P-scores")
  
  plt <- plt + geom_text(aes(x = dat$treat, y = dat$name, label = dat$pscore))
  
  plt
}
