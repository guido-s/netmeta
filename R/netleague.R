#' Create and print league table for network meta-analysis results
#' 
#' @description
#' A league table is a square matrix showing all pairwise comparisons
#' in a network meta-analysis. Typically, both treatment estimates and
#' confidence intervals are shown.
#' 
#' @aliases netleague print.netleague
#' 
#' @param x An object of class \code{netmeta} or \code{netleague}
#'   (mandatory).
#' @param y An object of class \code{netmeta} (optional).
#' @param comb.fixed A logical indicating whether a league table
#'   should be printed for the fixed effects (common effects) network
#'   meta-analysis.
#' @param comb.random A logical indicating whether a league table
#'   should be printed for the random effects network meta-analysis.
#' @param seq A character or numerical vector specifying the sequence
#'   of treatments in rows and columns of a league table.
#' @param ci A logical indicating whether confidence intervals should
#'   be shown.
#' @param backtransf A logical indicating whether printed results
#'   should be back transformed. If \code{backtransf = TRUE}, results
#'   for \code{sm = "OR"} are printed as odds ratios rather than log
#'   odds ratios, for example.
#' @param direct A logical indicating whether league table with
#'   network estimates (default) or estimates from direct comparisons
#'   should be generated if argument \code{y} is not missing.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param bracket A character with bracket symbol to print lower
#'   confidence interval: "[", "(", "\{", "".
#' @param separator A character string with information on separator
#'   between lower and upper confidence interval.
#' @param text.NA A character string to label missing values.
#' @param big.mark A character used as thousands separator.
#' @param \dots Additional arguments (ignored at the moment).
#' 
#' @details
#' A league table is a square matrix showing all pairwise comparisons
#' in a network meta-analysis. Typically, both treatment estimates and
#' confidence intervals are shown.
#' 
#' If argument \code{y} is not provided, the league table contains the
#' network estimates from network meta-analysis object \code{x} in the
#' lower triangle and the direct treatment estimates from pairwise
#' comparisons in the upper triangle.
#' 
#' If argument \code{y} is provided, the league table contains
#' information on treatment comparisons from network meta-analysis
#' object \code{x} in the lower triangle and from network
#' meta-analysis object \code{y} in the upper triangle. This is, for
#' example, useful to print information on efficacy and safety in the
#' same league table.
#' 
#' This implementation reports pairwise comparisons of the treatment
#' in the row versus the treatment in the column in the lower triangle
#' and column versus row in the upper triangle. This is a common
#' presentation for network meta-analyses which allows to easily
#' compare direction and magnitude of treatment effects. For example,
#' given treatments A, B, and C, the results reported in the first row
#' and second column as well as second row and first column are from
#' the pairwise comparison A versus B. Note, this presentation is
#' different from the printout of a network meta-analysis object which
#' reports opposite pairwise comparisons in the lower and upper
#' triangle, e.g., A versus B in the first row and second column and B
#' versus A in the second row and first column.
#' 
#' If the same network meta-analysis object is used for arguments
#' \code{x} and \code{y}, reciprocal treatment estimates will be shown
#' in the upper triangle (see examples), e.g., the comparison B versus
#' A.
#' 
#' R function \code{\link{netrank}} can be used to change the order of
#' rows and columns in the league table (see examples).
#'
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}, Gerta
#'   Rücker \email{ruecker@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netposet}},
#'   \code{\link{netrank}}
#' 
#' @keywords print
#' 
#' @examples
#' # Network meta-analysis of count mortality statistics
#' #
#' data(Woods2010)
#' 
#' p0 <- pairwise(treatment, event = r, n = N,
#'                studlab = author, data = Woods2010, sm = "OR")
#' net0 <- netmeta(p0)
#' 
#' oldopts <- options(width = 100)
#' 
#' # League table for fixed and random effects model with
#' # - network estimates in lower triangle
#' # - direct estimates in upper triangle
#' #
#' netleague(net0, digits = 2, bracket = "(", separator = " - ")
#' 
#' # League table for fixed effects model
#' #
#' netleague(net0, comb.random = FALSE, digits = 2)
#' 
#' # Change order of treatments according to treatment ranking (random
#' # effects model)
#' #
#' netleague(net0, comb.fixed = FALSE, digits = 2,
#'           seq = netrank(net0))
#' #
#' print(netrank(net0), comb.fixed = FALSE)
#' 
#' \dontrun{
#' # Create a CSV file with league table for random effects model
#' #
#' league0 <- netleague(net0, digits = 2, bracket = "(", separator = " to ")
#' #
#' write.table(league0$random, file = "league0-random.csv",
#'             row.names = FALSE, col.names = FALSE,
#'             sep = ",")
#' #
#' # Create Excel files with league tables (using R package WriteXLS
#' # which requires Perl https://www.perl.org/)
#' #
#' library(WriteXLS)
#' #
#' # League table from random effects model
#' #
#' WriteXLS(league0$random, ExcelFileName = "league0-random.xls",
#'          SheetNames = "leaguetable (random)", col.names = FALSE)
#' #
#' # League tables from fixed and random effects models
#' #
#' WriteXLS(list(league0$fixed, league0$random),
#'          ExcelFileName = "league0-both.xls",
#'          SheetNames = c("leaguetable (fixed)", "leaguetable (random)"),
#'          col.names = FALSE)
#' }
#' 
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
#' options(width = 200)
#' netleague(net1, digits = 2)
#' 
#' netleague(net1, digits = 2, ci = FALSE)
#' netleague(net2, digits = 2, ci = FALSE)
#' 
#' # League table for two outcomes with
#' # - network estimates of first outcome in lower triangle
#' # - network estimates of second outcome in upper triangle
#' #
#' netleague(net1, net2, digits = 2, ci = FALSE)
#' 
#' netleague(net1, net2, seq = netrank(net1, small = "bad"), ci = FALSE)
#' netleague(net1, net2, seq = netrank(net2, small = "bad"), ci = FALSE)
#' 
#' print(netrank(net1, small = "bad"))
#' print(netrank(net2, small = "bad"))
#' 
#' 
#' # Report results for network meta-analysis twice
#' #
#' netleague(net1, net1, seq = netrank(net1, small = "bad"), ci = FALSE,
#'           backtransf = FALSE)
#' netleague(net1, net1, seq = netrank(net1, small = "bad"), ci = FALSE,
#'           backtransf = FALSE, direct = TRUE)
#' 
#' options(oldopts)
#' 
#' \dontrun{
#' # Generate a partial order of treatment rankings 
#' #
#' np <- netposet(net1, net2, outcomes = outcomes, small.values = rep("bad",2))
#' hasse(np)
#' plot(np)
#' }
#' 
#' @rdname netleague
#' @export netleague


netleague <- function(x, y,
                      comb.fixed = x$comb.fixed, comb.random = x$comb.random,
                      seq = x$seq, ci = TRUE, backtransf = TRUE,
                      direct = FALSE,
                      digits = gs("digits"),
                      bracket = gs("CIbracket"),
                      separator = gs("CIseparator"),
                      text.NA = ".",
                      big.mark = gs("big.mark")) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  meta:::chkclass(x, "netmeta")
  ##
  chkchar <- meta:::chkchar
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  formatCI <- meta:::formatCI
  formatN <- meta:::formatN
  is.relative.effect <- meta:::is.relative.effect
  ##
  if (!missing(y)) {
    meta:::chkclass(y, "netmeta")
    ##
    if (length(x$seq) != length(y$seq))
      stop("Arguments 'x' and 'y' must have the same number of treatments.")
    if (any(sort(x$seq) != sort(y$seq)))
      stop("Arguments 'x' and 'y' must have the same treatments.")
    ##
    x.is.y <- all.equal(x, y)
    x.is.y <- is.logical(x.is.y) && x.is.y
  }
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  ##
  if (missing(seq)) {
    seq.f <- seq.r <- seq
  }
  else {
    if (is.null(seq))
      stop("Argument 'seq' must be not NULL.")
    else if (inherits(seq, "netrank")) {
      pscore.f <- seq$Pscore.fixed
      pscore.r <- seq$Pscore.random
      ##
      seq.f <- setseq(names(pscore.f)[rev(order(pscore.f))], x$seq)
      seq.r <- setseq(names(pscore.r)[rev(order(pscore.r))], x$seq)
    }
    else
      seq.f <- seq.r <- setseq(seq, x$seq)
  }
  ##
  chklogical(ci)
  chklogical(backtransf)
  chklogical(direct)
  chknumeric(digits, min = 0, single = TRUE)
  ##
  bracket.old <- gs("CIbracket")
  separator.old <- gs("CIseparator")
  cilayout(bracket, separator)
  on.exit(cilayout(bracket.old, separator.old))
  ##
  chkchar(text.NA)
  chkchar(big.mark)
  
  
  ##
  ##
  ## (2) Back-transform log odds ratios & Co
  ##
  ##
  if (!missing(y) & direct) {
    TE.fixed.x    <- x$TE.direct.fixed
    lower.fixed.x <- x$lower.direct.fixed
    upper.fixed.x <- x$upper.direct.fixed
    ##
    if (comb.random) {
      TE.random.x    <- x$TE.direct.random
      lower.random.x <- x$lower.direct.random
      upper.random.x <- x$upper.direct.random
    }
  }
  else {
    TE.fixed.x    <- x$TE.fixed
    lower.fixed.x <- x$lower.fixed
    upper.fixed.x <- x$upper.fixed
    ##
    if (comb.random) {
      TE.random.x    <- x$TE.random
      lower.random.x <- x$lower.random
      upper.random.x <- x$upper.random
    }
  }
  ##
  if (backtransf & is.relative.effect(x$sm)) {
    TE.fixed.x    <- exp(TE.fixed.x)
    lower.fixed.x <- exp(lower.fixed.x)
    upper.fixed.x <- exp(upper.fixed.x)
    ##
    if (comb.random) {
      TE.random.x    <- exp(TE.random.x)
      lower.random.x <- exp(lower.random.x)
      upper.random.x <- exp(upper.random.x)
    }
  }
  ##
  if (!missing(y)) {
    if (direct) {
      TE.fixed.y    <- y$TE.direct.fixed
      lower.fixed.y <- y$lower.direct.fixed
      upper.fixed.y <- y$upper.direct.fixed
      ##
      if (comb.random) {
        TE.random.y    <- y$TE.direct.random
        lower.random.y <- y$lower.direct.random
        upper.random.y <- y$upper.direct.random
      }
    }
    else {
      TE.fixed.y    <- y$TE.fixed
      lower.fixed.y <- y$lower.fixed
      upper.fixed.y <- y$upper.fixed
      ##
      if (comb.random) {
        TE.random.y    <- y$TE.random
        lower.random.y <- y$lower.random
        upper.random.y <- y$upper.random
      }
    }
    ##
    if (backtransf & is.relative.effect(y$sm)) {
      TE.fixed.y    <- exp(TE.fixed.y)
      lower.fixed.y <- exp(lower.fixed.y)
      upper.fixed.y <- exp(upper.fixed.y)
      ##
      if (comb.random) {
        TE.random.y    <- exp(TE.random.y)
        lower.random.y <- exp(lower.random.y)
        upper.random.y <- exp(upper.random.y)
      }
    }
    ##
    if (x.is.y) {
      TE.fixed.y    <- t(TE.fixed.y)
      lower.fixed.y <- t(lower.fixed.y)
      upper.fixed.y <- t(upper.fixed.y)
      ##
      if (comb.random) {
        TE.random.y    <- t(TE.random.y)
        lower.random.y <- t(lower.random.y)
        upper.random.y <- t(upper.random.y)
      }
    }
  }
  else {
    if (backtransf & is.relative.effect(x$sm)) {
      TE.fixed.y    <- exp(x$TE.direct.fixed)
      lower.fixed.y <- exp(x$lower.direct.fixed)
      upper.fixed.y <- exp(x$upper.direct.fixed)
      ##
      if (comb.random) {
        TE.random.y    <- exp(x$TE.direct.random)
        lower.random.y <- exp(x$lower.direct.random)
        upper.random.y <- exp(x$upper.direct.random)
      }
    }
    else {
      TE.fixed.y    <- x$TE.direct.fixed
      lower.fixed.y <- x$lower.direct.fixed
      upper.fixed.y <- x$upper.direct.fixed
      ##
      if (comb.random) {
        TE.random.y    <- x$TE.direct.random
        lower.random.y <- x$lower.direct.random
        upper.random.y <- x$upper.direct.random
      }
    }
  }
  ##
  ## Comparisons are column versus row
  ##
  TE.fixed.x <- t(TE.fixed.x)
  lower.fixed.x <- t(lower.fixed.x)
  upper.fixed.x <- t(upper.fixed.x)
  ##
  if (comb.random) {
    TE.random.x <- t(TE.random.x)
    lower.random.x <- t(lower.random.x)
    upper.random.x <- t(upper.random.x)
  }
  ##
  TE.fixed.y <- t(TE.fixed.y)
  lower.fixed.y <- t(lower.fixed.y)
  upper.fixed.y <- t(upper.fixed.y)
  ##
  if (comb.random) {
    TE.random.y <- t(TE.random.y)
    lower.random.y <- t(lower.random.y)
    upper.random.y <- t(upper.random.y)
  }
  
  
  ##
  ##
  ## (3) Print league table for fixed effect model
  ##
  ##
  TE.fixed.x    <- round(   TE.fixed.x[seq.f, seq.f], digits)
  lower.fixed.x <- round(lower.fixed.x[seq.f, seq.f], digits)
  upper.fixed.x <- round(upper.fixed.x[seq.f, seq.f], digits)
  ##
  if (ci) {
    nl.NA <- is.na(TE.fixed.x)
    nl.f <- paste(formatN(TE.fixed.x, digits = digits,
                          text.NA = text.NA, big.mark = big.mark),
                  formatCI(lower.fixed.x, upper.fixed.x, lab.NA = text.NA,
                           big.mark = big.mark))
    nl.f[nl.NA] <- text.NA
  }
  else
    nl.f <- formatN(TE.fixed.x, digits = digits,
                    text.NA = text.NA, big.mark = big.mark)
  ##
  nl.f <- matrix(nl.f, nrow = nrow(TE.fixed.x), ncol = ncol(TE.fixed.x))
  diag(nl.f) <- rownames(TE.fixed.x)
  ##
  TE.fixed.y    <- round(   TE.fixed.y[seq.f, seq.f], digits)
  lower.fixed.y <- round(lower.fixed.y[seq.f, seq.f], digits)
  upper.fixed.y <- round(upper.fixed.y[seq.f, seq.f], digits)
  ##
  if (ci) {
    nl.NA <- is.na(TE.fixed.y)
    nl.f.y <- paste(formatN(TE.fixed.y,
                            digits = digits, text.NA = text.NA, big.mark = big.mark),
                    formatCI(lower.fixed.y, upper.fixed.y, lab.NA = text.NA,
                             big.mark = big.mark))
    nl.f.y[nl.NA] <- text.NA
  }
  else
    nl.f.y <- formatN(TE.fixed.y, digits = digits,
                      text.NA = text.NA, big.mark = big.mark)
  ##
  nl.f.y <- matrix(nl.f.y, nrow = nrow(TE.fixed.y), ncol = ncol(TE.fixed.y))
  ##
  nl.f[upper.tri(nl.f)] <- t(nl.f.y)[upper.tri(nl.f)]
  ##
  nl.f <- as.data.frame(nl.f, stringsAsFactors = FALSE)
  
  
  ##
  ##
  ## (4) Print league table for random effects model
  ##
  ##
  if (comb.random) {
    TE.random.x    <- round(   TE.random.x[seq.r, seq.r], digits)
    lower.random.x <- round(lower.random.x[seq.r, seq.r], digits)
    upper.random.x <- round(upper.random.x[seq.r, seq.r], digits)
    ##
    if (ci) {
      nl.NA <- is.na(TE.random.x)
      nl.r <- paste(formatN(TE.random.x, digits = digits,
                            text.NA = text.NA, big.mark = big.mark),
                    formatCI(lower.random.x, upper.random.x, lab.NA = text.NA,
                             big.mark = big.mark))
      nl.r[nl.NA] <- text.NA
    }
    else
      nl.r <- formatN(TE.random.x, digits = digits,
                      text.NA = text.NA, big.mark = big.mark)
    ##
    nl.r <- matrix(nl.r, nrow = nrow(TE.random.x), ncol = ncol(TE.random.x))
    diag(nl.r) <- rownames(TE.random.x)
    ##
    TE.random.y    <- round(   TE.random.y[seq.r, seq.r], digits)
    lower.random.y <- round(lower.random.y[seq.r, seq.r], digits)
    upper.random.y <- round(upper.random.y[seq.r, seq.r], digits)
    ##
    if (ci) {
      nl.NA <- is.na(TE.random.y)
      nl.r.y <- paste(formatN(TE.random.y, digits = digits,
                              text.NA = text.NA),
                      formatCI(lower.random.y, upper.random.y, lab.NA = text.NA))
      nl.r.y[nl.NA] <- text.NA
    }
    else
      nl.r.y <- formatN(TE.random.y, digits = digits,
                        text.NA = text.NA, big.mark = big.mark)
    ##
    nl.r.y <- matrix(nl.r.y, nrow = nrow(TE.random.y), ncol = ncol(TE.random.y))
    ##
    nl.r[upper.tri(nl.r)] <- t(nl.r.y)[upper.tri(nl.f)]
    ##
    nl.r <- as.data.frame(nl.r, stringsAsFactors = FALSE)
  }
  
  
  res <- list(fixed = nl.f,
              random = if (comb.random) nl.r else NA,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              seq = seq, ci = ci, backtransf = backtransf,
              digits = digits)
  ##
  class(res) <- "netleague"
  
  
  res
}





#' @rdname netleague
#' @method print netleague
#' @export
#' @export print.netleague


print.netleague <- function(x,
                            comb.fixed = x$comb.fixed,
                            comb.random = x$comb.random,
                            ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  meta:::chkclass(x, "netleague")
  ##
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  
  
  ##
  ##
  ## (2) Print league table for fixed effect model
  ##
  ##
  if (comb.fixed) {
    cat("League table (fixed effect model):\n")
    ##
    prmatrix(x$fixed, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(x$fixed)),
             collab = rep("", ncol(x$fixed)))
    if (comb.random)
      cat("\n")
  }
  
  
  ##
  ##
  ## (3) Print league table for random effects model
  ##
  ##
  if (comb.random) {
    cat("League table (random effects model):\n")
    ##
    prmatrix(x$rando, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(x$random)),
             collab = rep("", ncol(x$random)))
  }
  
  
  invisible(NULL)
}
