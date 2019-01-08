#' Partial order of treatments in network meta-analysis
#' 
#' @description
#' Partial order of treatments in network meta-analysis. The set of
#' treatments in a network is called a partially ordered set (in
#' short, a \emph{poset}), if different outcomes provide different
#' treatment ranking lists.
#' 
#' @aliases netposet print.netposet
#' 
#' @param \dots See details.
#' @param outcomes A character vector with outcome names.
#' @param treatments A character vector with treatment names.
#' @param small.values See details.
#' @param comb.fixed A logical indicating whether to show results for
#'   the fixed effects (common effects) model.
#' @param comb.random A logical indicating whether to show results for
#'   the random effects model.
#' @param x An object of class \code{netposet}.
#' @param pooled A character string indicating whether Hasse diagram
#'   should be drawn for fixed effect (\code{"fixed"}) or random
#'   effects model (\code{"random"}). Can be abbreviated.
#' 
#' @details
#' In network meta-analysis, frequently different outcomes are
#' considered which may each provide a different ordering of
#' treatments. The concept of a partially ordered set (in short, a
#' \emph{poset}, Carlsen & Bruggemann, 2014) of treatments can be used
#' to gain further insights in situations with apparently conflicting
#' orderings. This implementation for rankings in network meta-analyis
#' is described in Rücker & Schwarzer (2017).
#' 
#' In function \code{netposet}, argument \code{\dots{}} can be any of the
#' following:
#' \itemize{
#' \item arbitrary number of \code{netrank} objects providing
#'   P-scores;
#' \item arbitrary number of \code{netmeta} objects;
#' \item single ranking matrix with each column providing P-scores
#'   (Rücker & Schwarzer 2015) or SUCRA values (Salanti et al. 2011)
#'   for an outcome and rows corresponding to treatments.
#' }
#' Note, albeit in general a ranking matrix is not constrained to have
#' values between 0 and 1, \code{netposet} stops with an error in this
#' case as this function expects a matrix with P-scores or SUCRA
#' values.
#' 
#' Argument \code{outcomes} can be used to label outcomes. If argument
#' \code{outcomes} is missing,
#' \itemize{
#' \item column names of the ranking matrix are used as outcome labels
#'   (if first argument is a ranking matrix and column names are
#'   available);
#' \item capital letters 'A', 'B', \dots{} are used as outcome labels
#'   and a corresponding warning is printed.
#' }
#'
#' Argument \code{treatments} can be used to provide treatment labels
#' if the first argument is a ranking matrix. If argument
#' \code{treatment} is missing,
#' \itemize{
#' \item row names of the ranking matrix are used as treatment labels
#'   (if available);
#' \item letters 'a', 'b', \dots{} are used as treatment labels and a
#'   corresponding warning is printed.
#' }
#' If argument \code{\dots{}} consists of \code{netmeta} objects,
#' \code{netrank} is called internally to calculate P-scores. In this
#' case, argument \code{small.values} can be used to specify for each
#' outcome whether small values are good or bad; see
#' \code{\link{netrank}}. This argument is ignored for a ranking
#' matrix and \code{netrank} objects.
#' 
#' Arguments \code{comb.fixed} and \code{comb.random} can be used to
#' define whether results should be printed and plotted for fixed
#' effect and / or random effects model. If netmeta and netrank
#' objects are provided in argument \code{\dots{}}, values for
#' \code{comb.fixed} and \code{comb.random} within these objects are
#' considered; if these values are not unique, argument
#' \code{comb.fixed} and / or \code{comb.random} are set to
#' \code{TRUE}.
#' 
#' In function \code{print.netposet}, argument \code{\dots{}} is
#' passed on to the printing function.
#'
#' @return
#' An object of class \code{netposet} with corresponding \code{print},
#' \code{plot}, and \code{hasse} functions. The object is a list
#' containing the following components:
#' \item{P.fixed}{Ranking matrix with rows corresponding to treatments
#'   and columns corresponding to outcomes (fixed effect model).}
#' \item{M0.fixed}{Hasse matrix skipping unnecessary paths (fixed
#'   effect model).}
#' \item{M.fixed}{"Full" Hasse matrix (fixed effect model).}
#' \item{O.fixed}{Matrix with information about partial ordering
#'   (fixed effect model).}
#' \item{P.random}{Ranking matrix with rows corresponding to
#'   treatments and columns corresponding to outcomes (random effects
#'   model).}
#' \item{M0.random}{Hasse matrix skipping unnecessary paths (random
#'   effects model).}
#' \item{M.random}{"Full" Hasse matrix (random effects model).}
#' \item{O.random}{Matrix with information about partial ordering
#'   (random effects model).}
#' \item{small.values, comb.fixed, comb.random}{As.defined above.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netmeta}}, \code{\link{netrank}},
#'   \code{\link{hasse}}, \code{\link{plot.netposet}}
#' 
#' @references
#' Carlsen L, Bruggemann R (2014):
#' Partial order methodology: a valuable tool in chemometrics.
#' \emph{Journal of Chemometrics},
#' \bold{28}, 226--34
#' 
#' Rücker G, Schwarzer G (2015):
#' Ranking treatments in frequentist network meta-analysis works
#' without resampling methods.
#' \emph{BMC Medical Research Methodology},
#' \bold{15}, 58
#' 
#' Rücker G, Schwarzer G (2017):
#' Resolve conflicting rankings of outcomes in network meta-analysis:
#' Partial ordering of treatments.
#' \emph{Research Synthesis Methods},
#' \bold{8}, 526--36
#' 
#' Salanti G, Ades AE, Ioannidis JP (2011):
#' Graphical methods and numerical summaries for presenting results
#' from multiple-treatment meta-analysis: an overview and tutorial.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{64}, 163--71
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
#' # Partial order of treatment rankings (two outcomes)
#' #
#' po <- netposet(netrank(net1, small.values = "bad"),
#'                netrank(net2, small.values = "bad"),
#'                outcomes = outcomes)
#' 
#' # Hasse diagram
#' #
#' hasse(po)
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
#' # Partial order of treatment rankings (all five outcomes)
#' #
#' po.ranks <- netposet(netrank(net1, small.values = "bad"),
#'                      netrank(net2, small.values = "bad"),
#'                      netrank(net3, small.values = "good"),
#'                      netrank(net4, small.values = "good"),
#'                      netrank(net5, small.values = "good"),
#'                      outcomes = outcomes)
#' 
#' # Same result
#' #
#' po.nets <- netposet(net1, net2, net3, net4, net5,
#'                     small.values = c("bad", "bad", "good", "good", "good"),
#'                     outcomes = outcomes)
#' #
#' all.equal(po.ranks, po.nets)
#' 
#' # Print matrix with P-scores (random effects model)
#' #
#' po.nets$P.random
#' 
#' # Hasse diagram for all outcomes (random effects model)
#' #
#' hasse(po.ranks)
#' 
#' # Hasse diagram for outcomes early response and early remission
#' #
#' po12 <- netposet(netrank(net1, small.values = "bad"),
#'                  netrank(net2, small.values = "bad"),
#'                  outcomes = outcomes[1:2])
#' hasse(po12)
#' 
#' # Scatter plot
#' #
#' oldpar <- par(pty = "s")
#' plot(po12)
#' par(oldpar)
#' }
#' 
#' # Example using ranking matrix with P-scores
#' #
#' # Ribassin-Majed L, Marguet S, Lee A.W., et al. (2017):
#' # What is the best treatment of locally advanced nasopharyngeal
#' # carcinoma? An individual patient data network meta-analysis.
#' # Journal of Clinical Oncology, 35, 498-505
#' #
#' outcomes <- c("OS", "PFS", "LC", "DC")
#' treatments <- c("RT", "IC-RT", "IC-CRT", "CRT",
#'                 "CRT-AC", "RT-AC", "IC-RT-AC")
#' #
#' # P-scores (from Table 1)
#' #
#' pscore.os  <- c(15, 33, 63, 70, 96, 28, 45) / 100
#' pscore.pfs <- c( 4, 46, 79, 52, 94, 36, 39) / 100
#' pscore.lc  <- c( 9, 27, 47, 37, 82, 58, 90) / 100
#' pscore.dc  <- c(16, 76, 95, 48, 72, 32, 10) / 100
#' #
#' pscore.matrix <- data.frame(pscore.os, pscore.pfs, pscore.lc, pscore.dc)
#' rownames(pscore.matrix) <- treatments
#' colnames(pscore.matrix) <- outcomes
#' pscore.matrix
#' #
#' po <- netposet(pscore.matrix)
#' po12 <- netposet(pscore.matrix[, 1:2])
#' po
#' po12
#' #
#' hasse(po)
#' hasse(po12)
#' #
#' oldpar <- par(pty = "s")
#' plot(po12)
#' par(oldpar)
#' 
#' @rdname netposet
#' @export netposet


netposet <- function(..., outcomes, treatments, small.values,
                     comb.fixed, comb.random) {
  
  
  args <- list(...)
  
  
  if (!missing(comb.fixed))
    meta:::chklogical(comb.fixed)
  if (!missing(comb.random))
    meta:::chklogical(comb.random)
  
  
  any.netmeta <- any.netrank <- FALSE
  
  
  ##
  ## First argument is expected to be a ranking matrix if only a
  ## single argument is provided
  ##
  if (length(args) == 1) {
    ##
    ## (1) Check arguments
    ##
    if (!(is.matrix(args[[1]]) | is.data.frame(args[[1]])))
      stop("First argument must be a matrix or ",
           "data frame if only a single argument is provided.",
           call. = FALSE)
    ##
    if (!missing(small.values))
      warning("Argument 'small.values' is ignored as first argument ",
              "is a ranking matrix.",
              call. = FALSE)
    ##
    pscore.matrix <- args[[1]]
    ##
    if (any(pscore.matrix[!is.na(pscore.matrix)] > 1) |
        any(pscore.matrix[!is.na(pscore.matrix)] < 0))
      stop("All elements of ranking matrix must be between 0 and 1.",
           call. = FALSE)
    ##
    n.outcomes <- ncol(pscore.matrix)
    n.treatments <- nrow(pscore.matrix)
    ##
    if (n.outcomes == 1)
      stop("Minimum number of two outcomes ",
           "(equal to number of columns of ranking matrix).",
           call. = FALSE)
    ##
    if (!missing(outcomes)) {
      if (length(outcomes) != n.outcomes)
        stop("Number of outcomes provided by argument 'outcomes' ",
             "is different from number of columns of ranking matrix.",
             call. = FALSE)
    }
    else
      outcomes <- colnames(pscore.matrix)
    ##
    if (!missing(treatments)) {
      if (length(treatments) != n.treatments)
        stop("Number of treatments provided by argument 'treatments' ",
             "is different from number of rows of ranking matrix.",
             call. = FALSE)
    }
    else
      treatments <- rownames(pscore.matrix)
    ##
    ## (2) P-Score matrices
    ##
    if (is.null(outcomes)) {
      warning("Outcomes are labelled 'A' to '", LETTERS[n.outcomes],
              "' as argument 'outcomes' is missing and column names are NULL.",
              call. = FALSE)
      outcomes <- LETTERS[1:n.outcomes]
    }
    ##
    if (is.null(treatments)) {
      warning("Treatments are labelled 'a' to '", letters[n.treatments],
              "' as row names are NULL.",
              call. = FALSE)
      treatments <- letters[1:n.treatments]
    }
    ##
    rownames(pscore.matrix) <- treatments
    colnames(pscore.matrix) <- outcomes
    ##
    pscore.matrix.fixed <- pscore.matrix
    pscore.matrix.random <- pscore.matrix
    ##
    if (missing(comb.fixed))
      comb.fixed <- FALSE
    if (missing(comb.random))
      comb.random <- FALSE
  }
  ##
  ## More than one element provided in '...' (netmeta or netrank objects)
  ##
  else {
    ##
    ## (1) Check arguments
    ##
    n.outcomes <- length(args)
    ##
    if (!missing(treatments))
      warning("Argument 'treatments' ignored (only needed if ",
              "first argument is a ranking matrix).",
              call. = FALSE)
    ##
    if (missing(outcomes)) {
      warning("Outcomes are labelled 'A' to '", LETTERS[n.outcomes],
              "' as argument 'outcomes' is missing.",
              call. = FALSE)
      outcomes <- LETTERS[1:n.outcomes]  
    }
    ##
    if (length(outcomes) != n.outcomes)
      stop("Number of outcomes provided by argument 'outcomes' ",
           "and number of rankings differ.",
           call. = FALSE)
    ##
    for (i in seq_along(args)) {
      if (inherits(args[[i]], "netmeta"))
        any.netmeta <- TRUE
      else if (inherits(args[[i]], "netrank"))
        any.netrank <- TRUE
      else
        stop("All elements of argument '...' must be of class ",
             "'netmeta' or 'netrank'.",
             call. = FALSE)
      ##
      if (!missing(small.values))
        small.values <- meta:::setchar(small.values, c("good", "bad"))
      else {
        if (any.netmeta)
          warning("R function netrank() called internally with argument ",
                  "small.values = 'good' as argument 'small.values' ",
                  "is missing.",
                  call. = FALSE)
        small.values <- rep("good", n.outcomes)
      }
      ##
      if (length(small.values) != n.outcomes)
        stop("Length of argument 'small.values' must be equal to ",
             "number of outcomes.",
             call. = FALSE)
    }
    ##
    ## (2) Extract P-Scores
    ##
    pscore.list.fixed <- pscore.list.random <- list()
    comb.fixeds <- comb.randoms <- rep_len(NA, length(args))
    ##
    for (i in seq_along(args)) {
      ##
      args.i <- args[[i]]
      ##
      if (inherits(args.i, "netmeta")) {
        pscore.list.fixed[[i]] <- netrank(args.i,
                                          small.values =
                                            small.values[i])$Pscore.fixed
        pscore.list.random[[i]] <- netrank(args.i,
                                           small.values =
                                             small.values[i])$Pscore.random
        comb.fixeds[i]  <- args.i$comb.fixed
        comb.randoms[i] <- args.i$comb.random
      }
      else if (inherits(args.i, "netrank")) {
        pscore.list.fixed[[i]] <- args.i$Pscore.fixed
        pscore.list.random[[i]] <- args.i$Pscore.random
        comb.fixeds[i]  <- args.i$x$comb.fixed
        comb.randoms[i] <- args.i$x$comb.random
      }
    }
    
    
    pscore.treatments <- lapply(pscore.list.fixed, names)
    n.treatments <- unlist(lapply(pscore.treatments, length))
    ##
    if (length(unique(n.treatments)) != 1) {
      sel.max <- seq_along(n.treatments)[n.treatments == max(n.treatments)][1]
      ##
      treatments <- pscore.treatments[[sel.max]]
      ##
      ## Act on rankings with missing treatments
      ##
      for (j in seq_along(n.treatments)[n.treatments < max(n.treatments)]) {
        treatments.j <- pscore.treatments[[j]]
        missing.j <- !(treatments %in% treatments.j)
        ##
        if (any(treatments.j != treatments[!missing.j]))
          stop("Treatment names of all rankings must be in same order.",
               call. = FALSE)
        ##
        pscore.j <- pscore.list.fixed[[j]]
        pscore.list.fixed[[j]] <- pscore.list.fixed[[sel.max]]
        pscore.list.fixed[[j]][treatments[!missing.j]] <- pscore.j
        pscore.list.fixed[[j]][treatments[missing.j]]  <- NA
        ##
        pscore.j <- pscore.list.random[[j]]
        pscore.list.random[[j]] <- pscore.list.random[[sel.max]]
        pscore.list.random[[j]][treatments[!missing.j]] <- pscore.j
        pscore.list.random[[j]][treatments[missing.j]]  <- NA
        ##
        pscore.treatments[[j]] <- treatments
        }
    }
    else
      treatments <- pscore.treatments[[1]]
    ##
    for (i in seq_along(pscore.treatments))
      if (any(pscore.treatments[[i]] != treatments)) {
        if (all(sort(pscore.treatments[[i]]) == sort(treatments)))
          stop("Different order of treatments provided:\n ",
               paste(paste("'", treatments, "'", sep = ""),
                     collapse = " - "),
               "\n ",
               paste(paste("'", pscore.treatments[[i]], "'", sep = ""),
                     collapse = " - "),
               call. = FALSE)
        else
          stop("Different treatments provided:\n ",
               paste(paste("'", treatments, "'", sep = ""),
                     collapse = " - "),
               "\n ",
               paste(paste("'", pscore.treatments[[i]], "'", sep = ""),
                     collapse = " - "),
               call. = FALSE)
      }
    ##
    pscore.matrix.fixed <- matrix(unlist(pscore.list.fixed,
                                         use.names = FALSE),
                                  ncol = length(outcomes), byrow = FALSE)
    pscore.matrix.random <- matrix(unlist(pscore.list.random,
                                          use.names = FALSE),
                                   ncol = length(outcomes), byrow = FALSE)
    rownames(pscore.matrix.fixed) <- treatments
    rownames(pscore.matrix.random) <- treatments
    colnames(pscore.matrix.fixed) <- outcomes
    colnames(pscore.matrix.random) <- outcomes
    ##
    text.netmeta.netrank <- if (any.netmeta) "netmeta" else ""
    text.netmeta.netrank <- if (any.netmeta & any.netrank) "netmeta and netrank"
    text.netmeta.netrank <- if (!any.netmeta & any.netrank) "netrank"
    text.netmeta.netrank <- paste(text.netmeta.netrank, "objects.")
    ##
    if (missing(comb.fixed)) {
      comb.fixed <- unique(comb.fixeds)
      if (length(comb.fixed) != 1) {
        warning("Argument 'comb.fixed' set to TRUE as different values are ",
                "available in ", text.netmeta.netrank,
                call. = FALSE)
        comb.fixed <- TRUE
      }
    }
    ##
    if (missing(comb.random)) {
      comb.random <- unique(comb.randoms)
      if (length(comb.random) != 1) {
        warning("Argument 'comb.random' set to TRUE as different values are ",
                "available in ", text.netmeta.netrank,
                call. = FALSE)
        comb.random <- TRUE
      }
    }
  } 
  
  
  n <- nrow(pscore.matrix.fixed)
  o <- ncol(pscore.matrix.fixed)
  
  
  Pos.fixed <- M.fixed <- matrix(0, nrow = n, ncol = n)
  ##
  rownames(Pos.fixed) <- colnames(Pos.fixed) <- rownames(pscore.matrix.fixed)
  rownames(M.fixed) <- colnames(M.fixed) <- rownames(pscore.matrix.fixed)
  ##
  Pos.random <- M.random <- matrix(0, nrow = n, ncol = n)
  ##
  rownames(Pos.random) <- colnames(Pos.random) <- rownames(pscore.matrix.random)
  rownames(M.random) <- colnames(M.random) <- rownames(pscore.matrix.random)
  
  
  ## Pos[i, j] counts how many rankings judge treatment i superior to
  ## treatment j (number between 0 and o)
  ##
  for (i in 1:n)
    for (j in 1:n)
      if (i != j)
        for (k in 1:o)
          if (!is.na(pscore.matrix.fixed[i, k]) &
              !is.na(pscore.matrix.fixed[j, k]) &
              pscore.matrix.fixed[i, k] >= pscore.matrix.fixed[j, k])
            Pos.fixed[i, j] <- Pos.fixed[i, j] + 1
  ##
  for (i in 1:n)
    for (j in 1:n)
      if (i != j)
        for (k in 1:o)
          if (!is.na(pscore.matrix.random[i, k]) &
              !is.na(pscore.matrix.random[j, k]) &
              pscore.matrix.random[i, k] >= pscore.matrix.random[j, k])
            Pos.random[i, j] <- Pos.random[i, j] + 1
  
  
  ## Calculate "full" Hasse matrix M
  ##
  M.fixed[Pos.fixed == o]   <- 1
  M.random[Pos.random == o] <- 1
  
  
  ## Matrix with information about partial ordering
  ##
  PO.fixed <- M.fixed[1:2, ]
  PO.fixed[1, ] <- rowSums(M.fixed)
  PO.fixed[2, ] <- colSums(M.fixed)
  rownames(PO.fixed) <- c("inferior", "superior")
  ##
  PO.random <- M.random[1:2, ]
  PO.random[1, ] <- rowSums(M.random)
  PO.random[2, ] <- colSums(M.random)
  rownames(PO.random) <- c("inferior", "superior")
  
  
  ## Skipping each direct path where a path of length 2 is found
  ##
  M0.fixed <- M.fixed - sign(M.fixed %*% M.fixed)
  ##
  M0.random <- M.random - sign(M.random %*% M.random)
  
  
  if (missing(small.values) | !any.netmeta | length(args) == 1)
    small.values <- NULL
  
  
  res <- list(P.fixed = pscore.matrix.fixed,
              M0.fixed = M0.fixed,
              M.fixed = M.fixed,
              O.fixed = PO.fixed,
              P.random = pscore.matrix.random,
              M0.random = M0.random,
              M.random = M.random,
              O.random = PO.random,
              small.values = small.values,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              call = match.call(),
              version = packageDescription("netmeta")$Version)
  ##
  class(res) <- "netposet"
  
  res
}





#' @rdname netposet
#' @method print netposet
#' @export
#' @export print.netposet


print.netposet <- function(x,
                           pooled = ifelse(x$comb.random, "random", "fixed"),
                           ...) {
  
  pooled <- meta:::setchar(pooled, c("fixed", "random"))
  
  if (pooled == "fixed") {
    if (!(!x$comb.fixed & !x$comb.random))
      cat("Fixed effect model\n")
    print(x$M0.fixed, ...)
    if (x$comb.random)
      cat("\n")
  }
  ##
  else {
    cat("Random effects model\n")
    print(x$M0.random, ...)
  }
  
  
  ## diag(M0) <- NA
  ## M0[is.na(M0)] <- "."
  ## prmatrix(M0, quote = FALSE, right = TRUE)
  
  invisible(NULL)
}
