#' Additive network meta-analysis for combinations of treatments
#' (disconnected networks)
#' 
#' @description
#' Some treatments in a network meta-analysis may be combinations of
#' other treatments or have common components. The influence of
#' individual components can be evaluated in an additive network
#' meta-analysis model assuming that the effect of treatment
#' combinations is the sum of the effects of its components. This
#' function implements this additive model in a frequentist way and is
#' particularly intended for disconnected networks.
#' 
#' @param TE Estimate of treatment effect, i.e. difference between
#'   first and second treatment (e.g. log odds ratio, mean difference,
#'   or log hazard ratio).
#' @param seTE Standard error of treatment estimate.
#' @param treat1 Label/Number for first treatment.
#' @param treat2 Label/Number for second treatment.
#' @param studlab An optional - but important! - vector with study
#'   labels (see \code{\link{netmeta}}).
#' @param data An optional data frame containing the study
#'   information.
#' @param subset An optional vector specifying a subset of studies to
#'   be used.
#' @param inactive A character string defining the inactive treatment
#'   (see Details).
#' @param sep.comps A single character to define separator between
#'   treatment components.
#' @param C.matrix C matrix (see Details).
#' @param sm A character string indicating underlying summary measure,
#'   e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
#'   \code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param level The level used to calculate confidence intervals for
#'   individual comparisons.
#' @param level.comb The level used to calculate confidence intervals
#'   for pooled estimates.
#' @param comb.fixed A logical indicating whether a fixed effects
#'   (common effects) network meta-analysis should be conducted.
#' @param comb.random A logical indicating whether a random effects
#'   network meta-analysis should be conducted.
#' @param reference.group Reference treatment.
#' @param baseline.reference A logical indicating whether results
#'   should be expressed as comparisons of other treatments versus the
#'   reference treatment (default) or vice versa. This argument is
#'   only considered if \code{reference.group} has been specified.
#' @param seq A character or numerical vector specifying the sequence
#'   of treatments in printouts.
#' @param tau.preset An optional value for the square-root of the
#'   between-study variance \eqn{\tau^2}.
#' @param tol.multiarm A numeric for the tolerance for consistency of
#'   treatment estimates and corresponding variances in multi-arm
#'   studies which are consistent by design.
#' @param details.chkmultiarm A logical indicating whether treatment
#'   estimates and / or variances of multi-arm studies with
#'   inconsistent results or negative multi-arm variances should be
#'   printed.
#' @param sep.trts A character used in comparison names as separator
#'   between treatment labels.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names (see Details).
#' @param title Title of meta-analysis / systematic review.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if studies are excluded from meta-analysis due to zero
#'   standard errors).
#' 
#' @details
#' Treatments in network meta-analysis (NMA) can be complex
#' interventions. Some treatments may be combinations of others or
#' have common components. The standard analysis provided by
#' \code{\link{netmeta}} is a NMA where all existing (single or
#' combined) treatments are considered as different nodes in the
#' network. Exploiting the fact that some treatments are combinations
#' of common components, an additive component network meta-analysis
#' (CNMA) model can be used to evaluate the influence of individual
#' components. This model assumes that the effect of a treatment
#' combination is the sum of the effects of its components which
#' implies that common components cancel out in comparisons.
#' 
#' This R function can be used for disconnected networks. Use
#' \code{\link{netmeta}} and \code{\link{netcomb}} for connected
#' networks.
#' 
#' The additive CNMA model has been implemented using Bayesian methods
#' (Mills et al., 2012; Welton et al., 2013). This function implements
#' the additive model in a frequentist way (Rücker et al., 2018).
#' 
#' The underlying multivariate model is given by
#' 
#' \deqn{\bold{\delta} = \bold{B} \bold{\theta}, \bold{\theta} =
#' \bold{C} \bold{\beta}}
#' 
#' with
#' \itemize{
#' \item[\eqn{\bold{\delta}}] vector of true treatment effects
#'   (differences) from individual studies,
#' \item[\eqn{\bold{B}}] is a design matrix describing the structure
#'   of the network,
#' \item[\eqn{\bold{\theta}}] parameter vector that represents the
#'   existing combined treatments,
#' \item[\eqn{\bold{C}}] matrix describing how the treatments are
#'   composed,
#' \item[\eqn{\bold{\beta}}] is a parameter vector representing the
#'   treatment components.
#' }
#' All parameters are estimated using weighted least squares
#' regression.
#' 
#' Argument \code{inactive} can be used to specify a single component
#' that does not have any therapeutic value. Accordingly, it is
#' assumed that the treatment effect of the combination of this
#' component with an additional treatment component is equal to the
#' treatment effect of the additional component alone.
#' 
#' Argument \code{sep.comps} can be used to specify the separator
#' between individual components. By default, the matrix \strong{C} is
#' calculated internally from treatment names. However, it is possible
#' to specify a different matrix using argument \code{C.matrix}.
#' 
#' @return
#' An object of classes \code{discomb} and \code{netcomb} with
#' corresponding \code{print}, \code{summary}, and \code{forest}
#' functions. The object is a list containing the following
#' components:
#' \item{inactive, sep.comps, C.matrix, sm}{As defined above.} 
#' \item{level.comb, comb.fixed, comb.random}{As defined above.}
#' \item{reference.group, baseline.reference, seq}{As defined above.}
#' \item{tau.preset, sep.trts, nchar.trts, backtransf}{As defined
#'   above.}
#' \item{k}{Total number of studies.}
#' \item{m}{Total number of pairwise comparisons.}
#' \item{n}{Total number of treatments.}
#' \item{c}{Total number of components.}
#' \item{comparisons.fixed, comparisons.random}{Lists with components
#'   studlab, treat1, treat2, TE, seTE, lower, upper, z, p level, df
#'   (fixed and random effects model).}
#' \item{components.fixed, components.random}{Lists with components
#'   TE, seTE, lower, upper, z, p level, df (fixed and random effects
#'   model).}
#' \item{combinations.fixed, combinations.random}{Lists with
#'   components TE, seTE, lower, upper, z, p level, df (fixed and
#'   random effects model).}
#' \item{Q.additive}{Overall heterogeneity / inconsistency statistic
#'   (additive model).}
#' \item{df.Q.additive}{Degrees of freedom for test of heterogeneity /
#'   inconsistency (additive model).}
#' \item{pval.Q.additive}{P-value for test of heterogeneity /
#'   inconsistency (additive model).}
#' \item{Q.standard}{Overall heterogeneity / inconsistency statistic
#'   (standard model).}
#' \item{df.Q.standard}{Degrees of freedom for test of heterogeneity /
#'   inconsistency (standard model).}
#' \item{pval.Q.standard}{P-value for test of heterogeneity /
#'   inconsistency (standard model).}
#' \item{Q.diff}{Test statistic for difference in goodness of fit
#'   between standard and additive model.}
#' \item{df.Q.diff}{Degrees of freedom for difference in goodness of
#'   fit between standard and additive model.}
#' \item{pval.Q.diff}{P-value for difference in goodness of fit
#'   between standard and additive model.}
#' \item{B.matrix}{Edge-vertex incidence matrix (\emph{m}x\emph{n}).}
#' \item{X}{Full design matrix (\emph{m}x\emph{n}).}
#' \item{trts}{Treatments included in network meta-analysis.}
#' \item{title}{Title of meta-analysis / systematic review.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \link{netcomb}, \link{forest.netcomb},
#'   \link{summary.netcomb}, \link{netmeta}, \link{netconnection}
#' 
#' @references
#' König J, Krahn U, Binder H (2013):
#' Visualizing the flow of evidence in network meta-analysis and
#' characterizing mixed treatment comparisons.
#' \emph{Statistics in Medicine},
#' \bold{32}, 5414--29
#' 
#' Mills EJ, Thorlund K, Ioannidis JP (2012):
#' Calculating additive treatment effects from multiple randomized
#' trials provides useful estimates of combination therapies.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{65}, 1282--8
#' 
#' Rücker G, Petropoulou M, Schwarzer G (2018):
#' Network meta-analysis of multicomponent interventions.
#' \emph{Manuscript submitted for publication}.
#' 
#' Welton NJ, Caldwell DM, Adamopoulos E, Vedhara K (2009):
#' Mixed treatment comparison meta-analysis of complex interventions:
#' psychological interventions in coronary heart disease.
#' \emph{American Journal of Epidemiology},
#' \bold{169}: 1158--65
#' 
#' @examples
#' # Artificial dataset
#' #
#' t1 <- c("A + B", "A + C", "A"    , "A"    , "D", "D", "E")
#' t2 <- c("C"    , "B"    , "B + C", "A + D", "E", "F", "F")
#' #
#' mean    <- c(4.1, 2.05, 0, 0, 0.1, 0.1, 0.05)
#' se.mean <- rep(0.1, 7)
#' #
#' study <- paste("study", c(1:4, 5, 5, 5))
#' #
#' dat <- data.frame(mean, se.mean, t1, t2, study,
#'                   stringsAsFactors = FALSE)
#' #
#' trts <- c("A", "A + B", "A + C", "A + D",
#'           "B", "B + C", "C", "D", "E", "F")
#' #
#' comps <- LETTERS[1:6]
#' 
#' # Use netconnection() to display network information
#' #
#' netconnection(t1, t2, study)
#' 
#' dc1 <- discomb(mean, se.mean, t1, t2, study, seq = trts)
#' dc1
#' 
#' forest(dc1, ref = "F")
#' 
#' # Define C matrix manually (which will produce the same results)
#' #
#' C <- rbind(c(1, 0, 0, 0, 0, 0),  # A
#'            c(1, 1, 0, 0, 0, 0),  # A + B
#'            c(1, 0, 1, 0, 0, 0),  # A + C
#'            c(1, 0, 0, 1, 0, 0),  # A + D
#'            c(0, 1, 0, 0, 0, 0),  # B
#'            c(0, 1, 1, 0, 0, 0),  # B + C
#'            c(0, 0, 1, 0, 0, 0),  # C
#'            c(0, 0, 0, 1, 0, 0),  # D
#'            c(0, 0, 0, 0, 1, 0),  # E
#'            c(0, 0, 0, 0, 0, 1))  # F
#' #                  
#' colnames(C) <- comps
#' rownames(C) <- trts
#' #
#' dc2 <- discomb(mean, se.mean, t1, t2, study, seq = trts,
#'                C.matrix = C)
#' #
#' # Compare C matrices
#' #
#' all.equal(dc1$C.matrix, dc2$C.matrix)
#' 
#' @export discomb


discomb <- function(TE, seTE,
                    treat1, treat2,
                    studlab, data = NULL, subset = NULL,
                    ##
                    inactive = NULL,
                    sep.comps = "+", 
                    C.matrix,
                    ##
                    sm,
                    level = gs("level"),
                    level.comb = gs("level.comb"),
                    comb.fixed = gs("comb.fixed"),
                    comb.random = gs("comb.random") | !is.null(tau.preset),
                    ##
                    reference.group = "",
                    baseline.reference = TRUE,
                    seq = NULL,
                    ##
                    tau.preset = NULL,
                    ##
                    tol.multiarm = 0.0005,
                    details.chkmultiarm = FALSE,
                    ##
                    sep.trts = ":",
                    nchar.trts = 666,
                    ##
                    backtransf = gs("backtransf"),
                    ##
                    title = "",
                    warn = TRUE) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkchar <- meta:::chkchar
  chklevel <- meta:::chklevel
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  ##
  chkchar(sep.comps, nchar = 1)
  ##
  chklevel(level)
  chklevel(level.comb)
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  ##
  chklogical(baseline.reference)
  ##
  if (!is.null(tau.preset))
    chknumeric(tau.preset, min = 0, single = TRUE)
  ##
  chknumeric(tol.multiarm, min = 0, single = TRUE)
  chklogical(details.chkmultiarm)
  ##
  missing.sep.trts <- missing(sep.trts)
  chkchar(sep.trts)
  chknumeric(nchar.trts, min = 1, single = TRUE)
  ##
  chklogical(backtransf)
  ##
  chkchar(title)
  chklogical(warn)
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  nulldata <- is.null(data)
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch TE, treat1, treat2, seTE, studlab from data:
  ##
  TE <- eval(mf[[match("TE", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  if (missing(sm))
    if (!is.null(data) && !is.null(attr(data, "sm")))
      sm <- attr(data, "sm")
    else
      sm <- ""
  ##
  seTE <- eval(mf[[match("seTE", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  ##
  treat1 <- eval(mf[[match("treat1", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  treat2 <- eval(mf[[match("treat2", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  ##
  k.Comp <- length(TE)
  ##
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  if (is.factor(studlab))
    studlab <- as.character(studlab)
  
  
  ##
  ##
  ## (3) Use subset for analysis
  ##
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  if (!is.null(subset)) {
    if ((is.logical(subset) & (sum(subset) > k.Comp)) ||
        (length(subset) > k.Comp))
      stop("Length of subset is larger than number of studies.")
    ##
    TE <- TE[subset]
    seTE <- seTE[subset]
    treat1 <- treat1[subset]
    treat2 <- treat2[subset]
    studlab <- studlab[subset]
  }
  ##
  labels <- sort(unique(c(treat1, treat2)))
  ##
  if (compmatch(labels, sep.trts)) {
    if (!missing.sep.trts)
      warning("Separator '", sep.trts, "' used in at least one treatment label. Try to use predefined separators: ':', '-', '_', '/', '+', '.', '|', '*'.")
    ##
    if (!compmatch(labels, ":"))
      sep.trts <- ":"
    else if (!compmatch(labels, "-"))
      sep.trts <- "-"
    else if (!compmatch(labels, "_"))
      sep.trts <- "_"
    else if (!compmatch(labels, "/"))
      sep.trts <- "/"
    else if (!compmatch(labels, "+"))
      sep.trts <- "+"
    else if (!compmatch(labels, "."))
      sep.trts <- "-"
    else if (!compmatch(labels, "|"))
      sep.trts <- "|"
    else if (!compmatch(labels, "*"))
      sep.trts <- "*"
    else
      stop("All predefined separators (':', '-', '_', '/', '+', '.', '|', '*') are used in at least one treatment label.",
           "\n   Please specify a different character that should be used as separator (argument 'sep.trts').",
           call. = FALSE)
  }
  ##
  if (!is.null(seq))
    seq <- setseq(seq, labels)
  else {
    seq <- labels
    if (is.numeric(seq))
      seq <- as.character(seq)
  }
  ##
  if (reference.group != "")
    reference.group <- setref(reference.group, labels)
  
  
  ##
  ##
  ## (4) Additional checks
  ##
  ##
  if (any(treat1 == treat2))
    stop("Treatments must be different (arguments 'treat1' and 'treat2').")
  ##
  if (length(studlab) != 0)
    studlab <- as.character(studlab)
  else {
    if (warn)
      warning("No information given for argument 'studlab'. Assuming that comparisons are from independent studies.")
    studlab <- seq(along = TE)
  }
  ##
  ## Check for correct number of comparisons
  ##
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)
      abs(x - round(x)) < tol
  ##
  tabnarms <- table(studlab)
  sel.narms <- !is.wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  ##
  if (sum(sel.narms) == 1)
    stop(paste("Study '", names(tabnarms)[sel.narms],
               "' has a wrong number of comparisons.",
               "\n  Please provide data for all treatment comparisons",
               " (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
               sep = ""))
  if (sum(sel.narms) > 1)
    stop(paste("The following studies have a wrong number of comparisons: ",
               paste(paste("'", names(tabnarms)[sel.narms], "'", sep = ""),
                     collapse = ", "),
               "\n  Please provide data for all treatment comparisons",
               " (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
               sep = ""))
  ##
  ## Check NAs and zero standard errors
  ##
  excl <- is.na(TE) | is.na(seTE) | seTE <= 0
  ##
  if (any(excl)) {
    dat.NAs <- data.frame(studlab = studlab[excl],
                          treat1 = treat1[excl],
                          treat2 = treat2[excl],
                          TE = format(round(TE[excl], 4)),
                          seTE = format(round(seTE[excl], 4))
                          )
    warning("Comparison",
            if (sum(excl) > 1) "s",
            " with missing TE / seTE or zero seTE not considered in network meta-analysis.",
            call. = FALSE)
    cat(paste("Comparison",
              if (sum(excl) > 1) "s",
              " not considered in network meta-analysis:\n", sep = ""))
    prmatrix(dat.NAs, quote = FALSE, right = TRUE,
             rowlab = rep("", sum(excl)))
    cat("\n")
    ##
    studlab <- studlab[!(excl)]
    treat1  <- treat1[!(excl)]
    treat2  <- treat2[!(excl)]
    TE      <- TE[!(excl)]
    seTE    <- seTE[!(excl)]
    ##
    seq <- seq[seq %in% unique(c(treat1, treat2))]
    labels <- labels[labels %in% unique(c(treat1, treat2))]
  }
  ##
  ## Check for correct number of comparisons (after removing
  ## comparisons with missing data)
  ##
  tabnarms <- table(studlab)
  sel.narms <- !is.wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  ##
  if (sum(sel.narms) == 1)
    stop(paste("After removing comparisons with missing treatment effects",
               " or standard errors,\n  study '",
               names(tabnarms)[sel.narms],
               "' has a wrong number of comparisons.",
               " Please check data and\n  consider to remove study",
               " from network meta-analysis.",
               sep = ""))
  if (sum(sel.narms) > 1)
    stop(paste("After removing comparisons with missing treatment effects",
               " or standard errors,\n  the following studies have",
               " a wrong number of comparisons: ",
               paste(paste("'", names(tabnarms)[sel.narms], "'", sep = ""),
                     collapse = ", "),
               "\n  Please check data and consider to remove studies",
               " from network meta-analysis.",
               sep = ""))
  ##
  ## Check for correct treatment order within comparison
  ##
  wo <- treat1 > treat2
  ##
  if (any(wo)) {
    if (warn)
      warning("Note, treatments within a comparison have been re-sorted in increasing order.", call. = FALSE)
    TE[wo] <- -TE[wo]
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
  }
  
  
  ##
  ##
  ## (5) Create C.matrix
  ##
  ##
  netc <- netconnection(treat1, treat2, studlab)
  ##
  if (missing(C.matrix)) {
    C.matrix <- createC(netc, sep.comps, inactive)
    C.matrix <- C.matrix[labels, , drop = FALSE]
  }
  else {
    ##
    if (!(is.matrix(C.matrix) | is.data.frame(C.matrix))) 
      stop("Argument 'C.matrix' must be a matrix or data frame.", 
           call. = FALSE)
    ##
    wrong.labels <- FALSE
    if (is.null(rownames(C.matrix)))
      wrong.labels <- TRUE
    else {
      if (length(unique(labels)) == length(unique(tolower(labels)))) 
        idx <- charmatch(tolower(rownames(C.matrix)), 
                         tolower(labels), nomatch = NA)
      else idx <- charmatch(rownames(C.matrix), labels, nomatch = NA)
      if (any(is.na(idx)) || any(idx == 0)) 
        wrong.labels <- TRUE
    }
    if (wrong.labels) 
      stop(paste("Row names of argument 'C.matrix' must be a ", 
                 "permutation of treatment names:\n  ",
                 paste(paste("'", labels, "'", sep = ""), collapse = " - "),
                 sep = ""),
           call. = FALSE)
    ##
    C.matrix <- C.matrix[labels, , drop = FALSE]
  }
  if (is.data.frame(C.matrix))
    C.matrix <- as.matrix(C.matrix)
  ##
  c <- ncol(C.matrix) # number of components
  
  
  p0 <- prepare(TE, seTE, treat1, treat2, studlab)
  ##
  o <- order(p0$order)
  ##
  B.matrix <- createB(p0$treat1.pos[o], p0$treat2.pos[o])
  ##
  colnames(B.matrix) <- labels
  rownames(B.matrix) <- studlab
  
  
  ##
  ## Design matrix based on treatment components
  ##
  X <- B.matrix %*% C.matrix
  ##
  colnames(X) <- colnames(C.matrix)
  rownames(X) <- studlab
  
  
  tdata <- data.frame(studies = p0$studlab[o], narms = p0$narms[o])
  tdata <- unique(tdata[order(tdata$studies, tdata$narms), ])
  ##
  studies <- tdata$studies
  narms <- tdata$narms
  n.a <- sum(narms)  
  ##
  comps <- colnames(C.matrix) # treatment components
  ##
  n <- length(labels)
  m <- length(TE)
  k <- length(unique(studlab))
  
  
  ##
  ## Fixed effects models
  ##
  df.Q.additive <- n.a - k - qr(X)$rank
  ##
  if (netc$n.subnets == 1) {
    net <- netmeta(TE, seTE, treat1, treat2, studlab)
    ##
    Q <- net$Q
    df.Q <- net$df.Q
    pval.Q <- net$pval.Q
    ##
    df.Q.diff <- n - 1 - qr(X)$rank
  }
  else {
    Q <- df.Q <- pval.Q <- NA
    df.Q.diff <- NA
  }
  
  
  res.f <- nma.additive(p0$TE[o], p0$weights[o], p0$studlab[o],
                        p0$treat1[o], p0$treat2[o], level.comb,
                        X, C.matrix, B.matrix,
                        Q, df.Q.additive, df.Q.diff,
                        n, sep.trts)
  
  
  ##
  ## Calculate heterogeneity statistics (additive model)
  ##
  Q.additive <- res.f$Q.additive
  ##
  if (!is.null(tau.preset))
    tau <- tau.preset
  else
    tau <- res.f$tau
  ##
  I2 <- res.f$I2
  
  
  ##
  ## Random effects models
  ##
  p1 <- prepare(TE, seTE, treat1, treat2, studlab, tau)
  ##
  res.r <- nma.additive(p1$TE[o], p1$weights[o], p1$studlab[o],
                        p1$treat1[o], p1$treat2[o], level.comb,
                        X, C.matrix, B.matrix,
                        Q, df.Q.additive, df.Q.diff,
                        n, sep.trts)
  
  
  NAs <- rep(NA, length(res.f$comparisons$TE))
  
  
  res <- list(studlab = p0$studlab[o],
              treat1 = p0$treat1[o],
              treat2 = p0$treat2[o],
              ##
              TE = p0$TE[o],
              seTE = seTE[o],
              seTE.adj = sqrt(1 / p0$weights[o]),
              ##
              event1 = NA,
              event2 = NA,
              n1 = NA,
              n2 = NA,
              ##
              k = k,
              m = m,
              n = n,
              d = NA,
              c = c,
              ##
              trts = labels,
              k.trts = NA,
              n.trts = NA,
              events.trts = NA,
              ##
              studies = NA,
              narms = NA,
              ##
              designs = NA,
              ##
              comps = names(res.f$components$TE),
              k.comps = NA,
              n.comps = NA,
              events.comps = NA,
              ##
              TE.nma.fixed = NAs,
              seTE.nma.fixed = NAs,
              lower.nma.fixed = NAs,
              upper.nma.fixed = NAs,
              zval.nma.fixed = NAs,
              pval.nma.fixed = NAs,
              ##
              TE.cnma.fixed = res.f$comparisons$TE,
              seTE.cnma.fixed = res.f$comparisons$seTE,
              lower.cnma.fixed = res.f$comparisons$lower,
              upper.cnma.fixed = res.f$comparisons$upper,
              zval.cnma.fixed = res.f$comparisons$z,
              pval.cnma.fixed = res.f$comparisons$p,
              ##
              TE.fixed = res.f$all.comparisons$TE,
              seTE.fixed = res.f$all.comparisons$seTE,
              lower.fixed = res.f$all.comparisons$lower,
              upper.fixed = res.f$all.comparisons$upper,
              zval.fixed = res.f$all.comparisons$z,
              pval.fixed = res.f$all.comparisons$p,
              ##
              TE.nma.random = NAs,
              seTE.nma.random = NAs,
              lower.nma.random = NAs,
              upper.nma.random = NAs,
              zval.nma.random = NAs,
              pval.nma.random = NAs,
              ##
              TE.cnma.random = res.r$comparisons$TE,
              seTE.cnma.random = res.r$comparisons$seTE,
              lower.cnma.random = res.r$comparisons$lower,
              upper.cnma.random = res.r$comparisons$upper,
              zval.cnma.random = res.r$comparisons$z,
              pval.cnma.random = res.r$comparisons$p,
              ##
              TE.random = res.r$all.comparisons$TE,
              seTE.random = res.r$all.comparisons$seTE,
              lower.random = res.r$all.comparisons$lower,
              upper.random = res.r$all.comparisons$upper,
              zval.random = res.r$all.comparisons$z,
              pval.random = res.r$all.comparisons$p,
              ##
              Comp.fixed = unname(res.f$components$TE),
              seComp.fixed = unname(res.f$components$seTE),
              lower.Comp.fixed = unname(res.f$components$lower),
              upper.Comp.fixed = unname(res.f$components$upper),
              zval.Comp.fixed = unname(res.f$components$z),
              pval.Comp.fixed = unname(res.f$components$p),
              ##
              Comp.random = unname(res.r$components$TE),
              seComp.random = unname(res.r$components$seTE),
              lower.Comp.random = unname(res.r$components$lower),
              upper.Comp.random = unname(res.r$components$upper),
              zval.Comp.random = unname(res.r$components$z),
              pval.Comp.random = unname(res.r$components$p),
              ##
              Comb.fixed = unname(res.f$combinations$TE),
              seComb.fixed = unname(res.f$combinations$seTE),
              lower.Comb.fixed = unname(res.f$combinations$lower),
              upper.Comb.fixed = unname(res.f$combinations$upper),
              zval.Comb.fixed = unname(res.f$combinations$z),
              pval.Comb.fixed = unname(res.f$combinations$p),
              ##
              Comb.random = unname(res.r$combinations$TE),
              seComb.random = unname(res.r$combinations$seTE),
              lower.Comb.random = unname(res.r$combinations$lower),
              upper.Comb.random = unname(res.r$combinations$upper),
              zval.Comb.random = unname(res.r$combinations$z),
              pval.Comb.random = unname(res.r$combinations$p),
              ##
              Q.additive = Q.additive, 
              df.Q.additive = df.Q.additive, 
              pval.Q.additive = res.f$pval.Q.additive,
              tau = tau,
              I2 = I2,
              ##
              Q.standard = Q,
              df.Q.standard = df.Q,
              pval.Q.standard = pval.Q, 
              ##
              Q.diff = res.f$Q.diff,
              df.Q.diff = df.Q.diff,
              pval.Q.diff = res.f$pval.Q.diff, 
              ##
              B.matrix = B.matrix,
              C.matrix = C.matrix,
              ##
              n.matrix = NA,
              events.matrix = NA,
              ##
              sm = sm,
              method = "Inverse",
              level = level,
              level.comb = level.comb,
              comb.fixed = comb.fixed,
              comb.random = comb.random, 
              ##
              reference.group = reference.group,
              baseline.reference = baseline.reference,
              all.treatments = NULL,
              seq = seq,
              ##
              tau.preset = tau.preset,
              ##
              sep.trts = sep.trts,
              nchar.trts = nchar.trts,
              ##
              inactive = inactive,
              sep.comps = sep.comps,
              ##
              backtransf = backtransf, 
              ##
              title = title,
              ##
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- c("discomb", "netcomb")
  
  
  res
}
