#' Create league table with network meta-analysis results
#' 
#' @description
#' A league table is a square matrix showing all pairwise comparisons
#' in a network meta-analysis. Typically, both treatment estimates and
#' confidence intervals are shown.
#' 
#' @aliases netleague print.netleague
#' 
#' @param x An object of created with \code{netmeta}, \code{netcomb},
#'   \code{discomb}, or \code{netleague} (mandatory).
#' @param y An object of class \code{netmeta}, \code{netcomb}, or \code{discomb}
#'   (optional).
#' @param common A logical indicating whether a league table should be
#'   printed for the common effects network meta-analysis.
#' @param random A logical indicating whether a league table should be
#'   printed for the random effects network meta-analysis.
#' @param seq A character or numerical vector specifying the sequence
#'   of treatments in rows and columns of a league table.
#' @param ci A logical indicating whether confidence intervals should
#'   be shown.
#' @param backtransf A logical indicating whether printed results
#'   should be back transformed. If \code{backtransf = TRUE}, results
#'   for \code{sm = "OR"} are printed as odds ratios rather than log
#'   odds ratios, for example.
#' @param direct A logical indicating whether league tables with
#'   network estimates or estimates from direct comparisons
#'   should be shown if (i) argument \code{y} is provided and the input for
#'   arguments \code{x} and \code{y} was created with \code{netmeta} or (ii)
#'   input for argument \code{x} was created with \code{netmeta} or
#'   \code{netcomb}.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param big.mark A character used as thousands separator.
#' @param text.NA A character string to label missing values.
#' @param bracket A character with bracket symbol to print lower
#'   confidence interval: "[", "(", "\{", "".
#' @param separator A character string with information on separator
#'   between lower and upper confidence interval.
#' @param lower.blank A logical indicating whether blanks between left
#'   bracket and lower confidence limit should be printed.
#' @param upper.blank A logical indicating whether blanks between
#'   separator and upper confidence limit should be printed.
#' @param writexl A logical indicating whether an Excel file should be
#'   created (R package \bold{writexl} must be available).
#' @param path A character string specifying the filename of the Excel
#'   file.
#' @param overwrite A logical indicating whether an existing Excel
#'   file should be overwritten.
#' @param details A logical specifying whether details on rows and columns
#'   should be printed.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (passed on to \code{write_xlsx}
#'   to create Excel file).
#' 
#' @details
#' A league table is a square matrix showing all pairwise comparisons
#' in a network meta-analysis (Hutton et al., 2015). Typically, both
#' treatment estimates and confidence intervals are shown.
#' 
#' If argument \code{y} is not provided, the league table contains
#' \itemize{
#' \item the network estimates in the lower triangle and the direct treatment
#'   estimates from pairwise comparisons in the upper triangle, if the input
#'   for argument \code{x} was created with \code{netmeta};
#' \item the network estimates from the component network meta-analysis in the
#'   lower triangle, if the input for argument \code{x} was created with
#'   \code{netcomb}, and, in the upper triangle, (i) the network estimates from
#'   the network meta-analysis, if argument \code{direct = FALSE}, or
#'   (ii)) the direct treatment estimates from pairwise comparisons, if
#'   argument \code{direct = TRUE};
#' \item the network estimates from the component network meta-analysis in the
#'   lower and upper triangle, if the input for argument \code{x} was created
#'   with \code{discomb}.
#' }
#' Note, for the random-effects model, the direct treatment estimates are based
#' on the common between-study variance \eqn{\tau^2} from the network
#' meta-analysis, i.e. the square of list element \code{x$tau}.
#' 
#' If argument \code{y} is provided, the league table contains
#' information on treatment comparisons from (component) network meta-analysis
#' object \code{x} in the lower triangle and from (component) network
#' meta-analysis object \code{y} in the upper triangle. This is, for
#' example, useful to print information on efficacy and safety in the
#' same league table. If argument \code{direct = TRUE}, direct estimates are
#' shown both in the lower and upper triangle if the input for arguments
#' \code{x} and \code{y} was created with \code{netmeta}.
#'
#' By default, an R object with the league tables is generated. Alternatively,
#' an Excel file is created if argument \code{writexl = TRUE}.
#' 
#' This implementation reports pairwise comparisons of the treatment
#' in the column versus the treatment in the row in the lower triangle
#' and row versus column in the upper triangle. This is a common
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
#' If the same (component) network meta-analysis object is used for arguments
#' \code{x} and \code{y}, reciprocal treatment estimates will be shown
#' in the upper triangle (see examples), e.g., the comparison B versus
#' A.
#' 
#' R function \code{\link{netrank}} can be used to change the order of
#' rows and columns in the league table (see examples).
#'
#' @return
#' An object of class \code{netleague} with corresponding \code{print}
#' function if \code{writexl = FALSE}. The object is a list containing
#' the league tables in list elements 'common' and 'random'. An Excel
#' file is created if \code{writexl = TRUE}. In this case, \code{NULL}
#' is returned in R.
#'
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}, Gerta
#'   Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netcomb}},
#'   \code{\link{discomb}}, \code{\link{netposet}}, \code{\link{netrank}},
#'   \code{\link[metadat]{dat.woods2010}},
#'   \code{\link[metadat]{dat.linde2015}}
#'
#' 
#' @references
#' Hutton B, Salanti G, Caldwell DM, et al. (2015):
#' The PRISMA Extension Statement for Reporting of Systematic Reviews
#' Incorporating Network Meta-analyses of Health Care Interventions:
#' Checklist and Explanations.
#' \emph{Annals of Internal Medicine},
#' \bold{162}, 777
#' 
#' @keywords print
#' 
#' @examples
#' # Network meta-analysis of count mortality statistics
#' #
#' p0 <- pairwise(treatment, event = r, n = N,
#'   studlab = author, data = dat.woods2010, sm = "OR")
#' net0 <- netmeta(p0)
#' 
#' oldopts <- options(width = 100)
#' 
#' # League table for common and random effects model with
#' # - network estimates in lower triangle
#' # - direct estimates in upper triangle
#' #
#' netleague(net0, digits = 2, bracket = "(", separator = " - ")
#' 
#' # League table for common effects model
#' #
#' netleague(net0, random = FALSE, digits = 2)
#' 
#' # Change order of treatments according to treatment ranking (random
#' # effects model)
#' #
#' netleague(net0, common = FALSE, digits = 2, seq = netrank(net0))
#' #
#' print(netrank(net0), common = FALSE)
#' 
#' \dontrun{
#' # Create a CSV file with league table for random effects model
#' #
#' league0 <- netleague(net0, digits = 2, bracket = "(", separator = " to ")
#' #
#' write.table(league0$random, file = "league0-random.csv",
#'   row.names = FALSE, col.names = FALSE, sep = ",")
#' #
#' # Create Excel files with league tables
#' # (if R package writexl is available)
#' #
#' netleague(net0, digits = 2, bracket = "(", separator = " to ",
#'           path = tempfile(fileext = ".xlsx"))
#' }
#' 
#' \donttest{
#' # Define order of treatments in depression dataset dat.linde2015
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
#'   studlab = id, data = dat.linde2015, sm = "OR")
#' #
#' net1 <- netmeta(p1, common = FALSE,
#'                 seq = trts, ref = "Placebo", small = "undesirable")
#' 
#' # (2) Early remission
#' #
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#' #
#' net2 <- netmeta(p2, common = FALSE,
#'                 seq = trts, ref = "Placebo", small = "undesirable")
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
#' netleague(net1, net2, seq = netrank(net1), ci = FALSE)
#' netleague(net1, net2, seq = netrank(net2), ci = FALSE)
#' 
#' netrank(net1)
#' netrank(net2)
#' 
#' 
#' # Report results for network meta-analysis twice
#' #
#' netleague(net1, net1, seq = netrank(net1), ci = FALSE,
#'   backtransf = FALSE)
#' netleague(net1, net1, seq = netrank(net1), ci = FALSE,
#'   backtransf = FALSE, direct = TRUE)
#' }
#' 
#' options(oldopts)
#' 
#' \dontrun{
#' # Generate a partial order of treatment rankings 
#' #
#' np <- netposet(net1, net2, outcomes = outcomes)
#' plot(np)
#'
#' hasse(np)
#' }
#' 
#' @rdname netleague
#' @export netleague

netleague <- function(x, y,
                      common = x$common, random = x$random,
                      seq = x$seq, ci = TRUE, backtransf = x$backtransf,
                      direct,
                      ##
                      digits = gs("digits"),
                      ##
                      big.mark = gs("big.mark"),
                      text.NA = ".",
                      ##
                      bracket = gs("CIbracket"),
                      separator = gs("CIseparator"),
                      lower.blank = gs("CIlower.blank"),
                      upper.blank = gs("CIupper.blank"),
                      ##
                      writexl = !missing(path),
                      path = "leaguetable.xlsx",
                      overwrite = FALSE,
                      #
                      details = gs("details"),
                      #
                      warn.deprecated = gs("warn.deprecated"),
                      ...) {
  
  
  ##
  ##
  ## (1) Check class and arguments
  ##
  ##
  
  chkclass(x, c("netmeta", "netcomb"))
  x <- updateversion(x)
  #
  x.netmeta <- inherits(x, "netmeta")
  x.discomb <- inherits(x, "discomb")
  x.netcomb <- !(x.netmeta | x.discomb)
  #
  x.is.y <- FALSE
  #
  avail.y <- !missing(y)
  #
  if (avail.y) {
    chkclass(y, c("netmeta", "netcomb"))
    ##
    if (length(x$seq) != length(y$seq))
      stop("Arguments 'x' and 'y' must have the same number of treatments.")
    if (any(sort(x$seq) != sort(y$seq)))
      stop("Arguments 'x' and 'y' must have the same treatments.")
    #
    y.netmeta <- inherits(y, "netmeta")
    y.discomb <- inherits(y, "discomb")
    y.netcomb <- !(y.netmeta | y.discomb)
    #
    x.is.y <- all.equal(x, y)
    x.is.y <- is.logical(x.is.y) && x.is.y
  }
  #
  if (missing(seq)) {
    seq.c <- seq.r <- seq
  }
  else {
    if (is.null(seq))
      stop("Argument 'seq' must be not NULL.")
    else if (inherits(seq, "netrank")) {
      ranking.c <- seq$ranking.common
      ranking.r <- seq$ranking.random
      ##
      seq.c <- setseq(names(ranking.c)[rev(order(ranking.c))], x$seq)
      seq.r <- setseq(names(ranking.r)[rev(order(ranking.r))], x$seq)
    }
    else
      seq.c <- seq.r <- setseq(seq, x$seq)
  }
  ##
  chklogical(ci)
  #
  chklogical(backtransf)
  #
  if (missing(direct)) {
    if (!avail.y)
      direct <- !x.discomb
    else
      direct <- FALSE
  }
  chklogical(direct)
  #
  if (direct &&
      (inherits(x, "discomb") || (avail.y && inherits(y, "discomb")))) {
    warning("Argument 'direct' set to FALSE for object created with discomb().")
    direct <- FALSE
  }
  #
  if (!direct && (inherits(x, "netmeta")) && !avail.y) {
    warning(paste("Argument 'direct' set to TRUE for single object",
                  "created with netmeta()."))
    direct <- TRUE
  }
  #
  chknumeric(digits, min = 0, length = 1)
  ##
  chkchar(text.NA)
  chkchar(big.mark)
  ##
  bracket.old <- gs("CIbracket")
  separator.old <- gs("CIseparator")
  lower.blank.old <- gs("CIlower.blank")
  upper.blank.old <- gs("CIupper.blank")
  ##
  cilayout(bracket, separator, lower.blank, upper.blank)
  on.exit(cilayout(bracket.old, separator.old,
                   lower.blank.old, upper.blank.old))
  ##
  chklogical(writexl)
  chkchar(path, length = 1)
  chklogical(overwrite)
  #
  chklogical(details)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  
  
  #
  #
  # (2) Set details
  #
  #
  
  if (details) {
    if (x.netmeta) {
      if (avail.y & direct)
        txt.x <- "\nLower triangle: results from direct comparisons"
      else
        txt.x <- "\nLower triangle: results from network meta-analysis"
    }
    else
      txt.x <- "\nLower triangle: results from component network meta-analysis"
    #
    if (avail.y) {
      if (y.netmeta) {
        if (direct)
          txt.y <- "\nUpper triangle: results from direct comparisons"
        else
          txt.y <- "\nUpper triangle: results from network meta-analysis"
      }
      else
        txt.y <-
          "\nUpper triangle: results from component network meta-analysis"
    }
    else {
      if (x.netmeta)
        txt.y <- "\nUpper triangle: results from direct comparisons"
      else if (x.netcomb) {
        if (direct)
          txt.y <- "\nUpper triangle: results from direct comparisons"
        else
          txt.y <- "\nUpper triangle: results from network meta-analysis"
      }
      else {
          txt.y <-
            "\nUpper triangle: results from component network meta-analysis"
      }
    }    
    #
    txt.x <- paste(txt.x, "(column vs row)")
    txt.y <- paste(txt.y,
                   if (x.discomb & !avail.y)
                     "(column vs row)"
                   else if (!avail.y | (avail.y & !x.is.y)) "(row vs column)"
                   else "(column vs row)")
    #
    text.details <- paste0(txt.x, txt.y, "\n")
  }
  else
    text.details <- NULL
  
  
  ##
  ##
  ## (3) Values in lower triangle
  ##
  ##
  
  if (x.netmeta & avail.y & direct) {
    TE.common.x    <- x$TE.direct.common
    lower.common.x <- x$lower.direct.common
    upper.common.x <- x$upper.direct.common
    #
    TE.random.x    <- x$TE.direct.random
    lower.random.x <- x$lower.direct.random
    upper.random.x <- x$upper.direct.random
  }
  else if (x.netcomb & avail.y & direct) {
    TE.common.x    <- x$x$TE.direct.common
    lower.common.x <- x$x$lower.direct.common
    upper.common.x <- x$x$upper.direct.common
    #
    TE.random.x    <- x$x$TE.direct.random
    lower.random.x <- x$x$lower.direct.random
    upper.random.x <- x$x$upper.direct.random
  }
  else {
    TE.common.x    <- x$TE.common
    lower.common.x <- x$lower.common
    upper.common.x <- x$upper.common
    #
    TE.random.x    <- x$TE.random
    lower.random.x <- x$lower.random
    upper.random.x <- x$upper.random
  }
  
  
  #
  #
  # (4) Values in upper triangle
  #
  #
  
  if (!avail.y) {
    if (x.netmeta) {
      TE.common.y    <- x$TE.direct.common
      lower.common.y <- x$lower.direct.common
      upper.common.y <- x$upper.direct.common
      #
      TE.random.y    <- x$TE.direct.random
      lower.random.y <- x$lower.direct.random
      upper.random.y <- x$upper.direct.random
    }
    else if (inherits(x, "discomb")) {
      TE.common.y    <- x$TE.common
      lower.common.y <- x$lower.common
      upper.common.y <- x$upper.common
      #
      TE.random.y    <- x$TE.random
      lower.random.y <- x$lower.random
      upper.random.y <- x$upper.random
    }
    else {
      if (direct) {
        TE.common.y    <- x$x$TE.direct.common
        lower.common.y <- x$x$lower.direct.common
        upper.common.y <- x$x$upper.direct.common
        #
        TE.random.y    <- x$x$TE.direct.random
        lower.random.y <- x$x$lower.direct.random
        upper.random.y <- x$x$upper.direct.random
      }
      else {
        TE.common.y    <- x$x$TE.common
        lower.common.y <- x$x$lower.common
        upper.common.y <- x$x$upper.common
        #
        TE.random.y    <- x$x$TE.random
        lower.random.y <- x$x$lower.random
        upper.random.y <- x$x$upper.random
      }
    }
  }
  else {
    if (y.netmeta & direct) {
      TE.common.y    <- y$TE.direct.common
      lower.common.y <- y$lower.direct.common
      upper.common.y <- y$upper.direct.common
      #
      TE.random.y    <- y$TE.direct.random
      lower.random.y <- y$lower.direct.random
      upper.random.y <- y$upper.direct.random
    }
    else if (y.netcomb & direct) {
      TE.common.y    <- y$x$TE.direct.common
      lower.common.y <- y$x$lower.direct.common
      upper.common.y <- y$x$upper.direct.common
      #
      TE.random.y    <- y$x$TE.direct.random
      lower.random.y <- y$x$lower.direct.random
      upper.random.y <- y$x$upper.direct.random
    }
    else {
      TE.common.y    <- y$TE.common
      lower.common.y <- y$lower.common
      upper.common.y <- y$upper.common
      #
      TE.random.y    <- y$TE.random
      lower.random.y <- y$lower.random
      upper.random.y <- y$upper.random
    }
  }
  #
  if (x.is.y | (x.discomb & !avail.y)) {
    TE.common.y    <- t(TE.common.y)
    lower.common.y <- t(lower.common.y)
    upper.common.y <- t(upper.common.y)
    #
    TE.random.y    <- t(TE.random.y)
    lower.random.y <- t(lower.random.y)
    upper.random.y <- t(upper.random.y)
  }
  
  
  #
  #
  # (5) Back-transform results
  #
  #
  
  if (backtransf) {
    TE.common.x    <- backtransf(TE.common.x, x$sm)
    lower.common.x <- backtransf(lower.common.x, x$sm)
    upper.common.x <- backtransf(upper.common.x, x$sm)
    #
    # Switch lower and upper limit for VE if results have been
    # backtransformed
    #
    if (x$sm == "VE") {
      tmp.l <- lower.common.x
      lower.common.x <- upper.common.x
      upper.common.x <- tmp.l
    }
    #
    if (random) {
      TE.random.x    <- backtransf(TE.random.x, x$sm)
      lower.random.x <- backtransf(lower.random.x, x$sm)
      upper.random.x <- backtransf(upper.random.x, x$sm)
      #
      if (x$sm == "VE") {
        tmp.l <- lower.random.x
        lower.random.x <- upper.random.x
        upper.random.x <- tmp.l
      }
    }
    #
    TE.common.y    <- backtransf(TE.common.y, x$sm)
    lower.common.y <- backtransf(lower.common.y, x$sm)
    upper.common.y <- backtransf(upper.common.y, x$sm)
    #
    # Switch lower and upper limit for VE if results have been
    # backtransformed
    #
    if (x$sm == "VE") {
      tmp.l <- lower.common.y
      lower.common.y <- upper.common.y
      upper.common.y <- tmp.l
    }
    #
    if (random) {
      TE.random.y    <- backtransf(TE.random.y, x$sm)
      lower.random.y <- backtransf(lower.random.y, x$sm)
      upper.random.y <- backtransf(upper.random.y, x$sm)
      #
      if (x$sm == "VE") {
        tmp.l <- lower.random.y
        lower.random.y <- upper.random.y
        upper.random.y <- tmp.l
      }
    }
  }
  
  
  #
  #
  # (6) Comparisons are column versus row
  #
  #
  
  TE.common.x <- t(TE.common.x)
  lower.common.x <- t(lower.common.x)
  upper.common.x <- t(upper.common.x)
  #
  TE.random.x <- t(TE.random.x)
  lower.random.x <- t(lower.random.x)
  upper.random.x <- t(upper.random.x)
  #
  TE.common.y <- t(TE.common.y)
  lower.common.y <- t(lower.common.y)
  upper.common.y <- t(upper.common.y)
  #
  TE.random.y <- t(TE.random.y)
  lower.random.y <- t(lower.random.y)
  upper.random.y <- t(upper.random.y)
  
  
  #
  #
  # (7) Create league table for common effects model
  #
  #
  
  TE.common.x    <- round(   TE.common.x[seq.c, seq.c], digits)
  lower.common.x <- round(lower.common.x[seq.c, seq.c], digits)
  upper.common.x <- round(upper.common.x[seq.c, seq.c], digits)
  #
  if (ci) {
    nl.NA <- is.na(TE.common.x)
    nl.c <- paste(formatN(TE.common.x, digits = digits,
                          text.NA = text.NA, big.mark = big.mark),
                  formatCI(lower.common.x, upper.common.x, lab.NA = text.NA,
                           big.mark = big.mark))
    nl.c[nl.NA] <- text.NA
  }
  else
    nl.c <- formatN(TE.common.x, digits = digits,
                    text.NA = text.NA, big.mark = big.mark)
  ##
  nl.c <- matrix(nl.c, nrow = nrow(TE.common.x), ncol = ncol(TE.common.x))
  diag(nl.c) <- rownames(TE.common.x)
  ##
  TE.common.y    <- round(   TE.common.y[seq.c, seq.c], digits)
  lower.common.y <- round(lower.common.y[seq.c, seq.c], digits)
  upper.common.y <- round(upper.common.y[seq.c, seq.c], digits)
  ##
  if (ci) {
    nl.NA <- is.na(TE.common.y)
    nl.c.y <- paste(formatN(TE.common.y,
                            digits = digits, text.NA = text.NA,
                            big.mark = big.mark),
                    formatCI(lower.common.y, upper.common.y, lab.NA = text.NA,
                             big.mark = big.mark))
    nl.c.y[nl.NA] <- text.NA
  }
  else
    nl.c.y <- formatN(TE.common.y, digits = digits,
                      text.NA = text.NA, big.mark = big.mark)
  ##
  nl.c.y <- matrix(nl.c.y, nrow = nrow(TE.common.y), ncol = ncol(TE.common.y))
  ##
  nl.c[upper.tri(nl.c)] <- t(nl.c.y)[upper.tri(nl.c)]
  ##
  nl.c <- as.data.frame(nl.c, stringsAsFactors = FALSE)
  #
  names(nl.c) <- rownames(nl.c) <- seq.c
  
  
  #
  #
  # (8) Create league table for random effects model
  #
  #
  
  if (random) {
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
    nl.r[upper.tri(nl.r)] <- t(nl.r.y)[upper.tri(nl.c)]
    ##
    nl.r <- as.data.frame(nl.r, stringsAsFactors = FALSE)
    #
    names(nl.r) <- rownames(nl.r) <- seq.r
  }
  
  
  #
  #
  # (9) Save Excel file
  #
  #
  
  if (writexl) {
    if (!(common | random)) {
      warning("Excel file not generated as neither ",
              "argument 'common' nor 'random' is TRUE.")
      return(invisible(NULL))
    }
    ##
    if (!is_installed_package("writexl", stop = FALSE))
      stop(paste0("Package 'writexl' missing.",
                  "\n  ",
                  "Please use the following R command for installation:",
                  "\n  install.packages(\"writexl\")"),
           call. = FALSE)
    ##
    if (file.exists(path) & !overwrite)
      warning("File '", path, "' exists. ",
              "Use argument 'overwrite = TRUE' to overwrite file.",
              call. = FALSE)
    else {
      if (common & random)
        xlsx <- list(common = nl.c, random = nl.r)
      else if (common)
        xlsx <- list(common = nl.c)
      else
        xlsx <- list(random = nl.r)
      ##
      writexl::write_xlsx(xlsx, path = path, col_names = FALSE, ...)
      message(paste0("League table", if (common & random) "s",
                     " saved in file '", path, "'."))
    }
    ##
    return(invisible(NULL))
  }
  
  
  #
  #
  # (10) Return league tables
  #
  #
  
  res <- list(common = nl.c, random = if (random) nl.r else NA,
              seq = seq, ci = ci, backtransf = backtransf,
              x = list(common = common, random = random),
              details = text.details,
              direct = direct,
              digits = digits,
              version = packageDescription("netmeta")$Version)
  #
  # Backward compatibility
  #
  res$fixed <- res$common
  res$x$fixed <- res$x$common
  ##
  class(res) <- "netleague"
  
  res
}


#' @rdname netleague
#' @method print netleague
#' @export

print.netleague <- function(x,
                            common = x$x$common,
                            random = x$x$random,
                            warn.deprecated = gs("warn.deprecated"),
                            ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkclass(x, "netleague")
  x <- updateversion(x)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  
  
  ##
  ##
  ## (2) Print league table for common effects model
  ##
  ##
  if (common) {
    cat("League table (common effects model):\n")
    ##
    prmatrix(x$common, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(x$common)),
             collab = rep("", ncol(x$common)))
    if (random)
      cat("\n")
  }
  
  
  ##
  ##
  ## (3) Print league table for random effects model
  ##
  ##
  if (random) {
    cat("League table (random effects model):\n")
    ##
    prmatrix(x$random, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(x$random)),
             collab = rep("", ncol(x$random)))
  }
  
  if (!is.null(x$details))
    cat(x$details)
  
  invisible(NULL)
}
