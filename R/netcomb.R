#' Additive network meta-analysis for combinations of treatments
#' 
#' @description
#' Some treatments in a network meta-analysis may be combinations of
#' other treatments or have common components. The influence of
#' individual components can be evaluated in an additive network
#' meta-analysis model assuming that the effect of treatment
#' combinations is the sum of the effects of its components. This
#' function implements this additive model in a frequentist way.
#' 
#' @param x An object of class \code{netmeta}.
#' @param inactive A character string defining the inactive treatment
#'   component (see Details).
#' @param sep.comps A single character to define separator between
#'   treatment components.
#' @param C.matrix C matrix (see Details).
#' @param common A logical indicating whether a common effects network
#'   meta-analysis should be conducted.
#' @param random A logical indicating whether a random effects network
#'   meta-analysis should be conducted.
#' @param tau.preset An optional value for the square-root of the
#'   between-study variance \eqn{\tau^2}.
#' @param overall.hetstat A logical indicating whether to print heterogeneity
#'   measures.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots.
#' @param details.chkident A logical indicating whether details on
#'   unidentifiable components should be printed.
#' @param nchar.comps A numeric defining the minimum number of
#'   characters used to create unique names for components (see
#'   Details).
#' @param func.inverse R function used to calculate the pseudoinverse
#'   of the Laplacian matrix L (see \code{\link{netmeta}}).
#' @param na.unident A logical indicating whether unidentifiable
#'   components and combinations should be set to missing values.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
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
#' The additive CNMA model has been implemented using Bayesian methods
#' (Mills et al., 2012; Welton et al., 2013). This function implements
#' the additive model in a frequentist way (Rücker et al., 2020).
#' 
#' The underlying multivariate model is given by
#' 
#' \deqn{\bold{\delta} = \bold{B} \bold{\theta}, \bold{\theta} =
#' \bold{C} \bold{\beta}}
#' 
#' with \describe{
#' \item{\eqn{\bold{\delta}}}{vector of true treatment effects
#'   (differences) from individual studies,}
#' \item{\eqn{\bold{B}}}{design matrix describing the structure of the
#'   network,}
#' \item{\eqn{\bold{\theta}}}{parameter vector that represents the
#'   existing combined treatments,}
#' \item{\eqn{\bold{C}}}{matrix describing how the treatments are
#'   composed,}
#' \item{\eqn{\bold{\beta}}}{parameter vector representing the
#'   treatment components.}
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
#' By default, component names are not abbreviated in
#' printouts. However, in order to get more concise printouts,
#' argument \code{nchar.comps} can be used to define the minimum
#' number of characters for abbreviated component names (see
#' \code{\link{abbreviate}}, argument \code{minlength}). R function
#' \code{\link{treats}} is utilised internally to create abbreviated
#' component names.
#' 
#' @note
#' This function calculates effects for individual components and
#' complex interventions present in the network.
#'
#' R function \code{\link{netcomplex}} can be used to calculate the
#' effect for arbitrary complex interventions in a component network
#' meta-analysis. Furthermore, R function \code{\link{netcomparison}}
#' can be used to calculate the effect for comparisons of two
#' arbitrary complex intervention in a component network
#' meta-analysis.
#'
#' @return
#' An object of class \code{netcomb} with corresponding \code{print},
#' \code{summary}, and \code{forest} functions. The object is a list
#' containing the following components:
#' \item{studlab}{Study labels.}
#' \item{treat1}{Label/Number for first treatment.}
#' \item{treat2}{Label/Number for second treatment.}
#' \item{TE}{Estimate of treatment effect, i.e. difference between
#'   first and second treatment.}
#' \item{seTE}{Standard error of treatment estimate.}
#' \item{seTE.adj.common, seTE.adj.random}{Standard error of treatment
#'   estimate, adjusted for multi-arm studies.}
#' \item{design}{Design of study providing pairwise comparison.}
#' \item{event1}{Number of events in first treatment group.}
#' \item{event2}{Number of events in second treatment group.}
#' \item{n1}{Number of observations in first treatment group.}
#' \item{n2}{Number of observations in second treatment group.}
#' \item{k}{Total number of studies.}
#' \item{m}{Total number of pairwise comparisons.}
#' \item{n}{Total number of treatments.}
#' \item{d}{Total number of designs (corresponding to the unique set
#'   of treatments compared within studies).}
#' \item{c}{Total number of components.}
#' \item{trts}{Treatments included in network meta-analysis.}
#' \item{k.trts}{Number of studies evaluating a treatment.}
#' \item{n.trts}{Number of observations receiving a treatment (if
#'   arguments \code{n1} and \code{n2} are provided).}
#' \item{events.trts}{Number of events observed for a treatment (if
#'   arguments \code{event1} and \code{event2} are provided).}
#' \item{studies}{Study labels coerced into a factor with its levels
#'   sorted alphabetically.}
#' \item{narms}{Number of arms for each study.}
#' \item{designs}{Unique list of designs present in the network. A
#'   design corresponds to the set of treatments compared within a
#'   study.}
#' \item{comps}{Unique list of components present in the network.}
#' \item{TE.nma.common, TE.nma.random}{A vector of length \emph{m} of
#'   consistent treatment effects estimated by network meta-analysis
#'   (nma) (common and random effects model).}
#' \item{seTE.nma.common, seTE.nma.random}{A vector of length \emph{m}
#'   of effective standard errors estimated by network meta-analysis
#'   (common and random effects model).}
#' \item{lower.nma.common, lower.nma.random}{A vector of length
#'   \emph{m} of lower confidence interval limits for consistent
#'   treatment effects estimated by network meta-analysis (common
#'   and random effects model).}
#' \item{upper.nma.common, upper.nma.random}{A vector of length
#'   \emph{m} of upper confidence interval limits for the consistent
#'   treatment effects estimated by network meta-analysis (common
#'   and random effects model).}
#' \item{statistic.nma.common, statistic.nma.random}{A vector of length \emph{m}
#'   of z-values for test of treatment effect for individual
#'   comparisons (common and random effects model).}
#' \item{pval.nma.common, pval.nma.random}{A vector of length \emph{m}
#'   of p-values for test of treatment effect for individual
#'   comparisons (common and random effects model).}
#' \item{TE.cnma.common, TE.cnma.random}{A vector of length \emph{m} of
#'   consistent treatment effects estimated by the additive (common and
#'   random effects) model.}
#' \item{seTE.cnma.common, seTE.cnma.random}{A vector of length
#'   \emph{m} with standard errors estimated by the additive (common
#'   and random effects) model.}
#' \item{lower.cnma.common, lower.cnma.random}{A vector of length
#'   \emph{m} of lower confidence interval limits for consistent
#'   treatment effects estimated by the additive (common and random
#'   effects) model.}
#' \item{upper.cnma.common, upper.cnma.random}{A vector of length
#'   \emph{m} of upper confidence interval limits for consistent
#'   treatment effects estimated by the additive (common and random
#'   effects) model.}
#' \item{statistic.cnma.common, statistic.cnma.random}{A vector of length
#'   \emph{m} of z-values for the test of an overall effect estimated
#'   by the additive (common and random effects) model.}
#' \item{pval.cnma.common, pval.cnma.random}{A vector of length
#'   \emph{m} of p-values for the test of an overall effect estimated
#'   by the additive (common and random effects) model.}
#' \item{TE.common, TE.random}{\emph{n}x\emph{n} matrix with overall
#'   treatment effects estimated by the additive (common and random
#'   effects) model.}
#' \item{seTE.common, seTE.random}{\emph{n}x\emph{n} matrix with
#'   standard errors estimated by the additive (common and random
#'   effects) model.}
#' \item{lower.common, upper.common, lower.random,
#'   upper.random}{\emph{n}x\emph{n} matrices with lower and upper
#'   confidence interval limits estimated by the additive (common and
#'   random effects) model.}
#' \item{statistic.common, pval.common, statistic.random,
#'   pval.random}{\emph{n}x\emph{n} matrices with z-values and
#'   p-values for test of overall effect estimated by the additive
#'   (common and random effects) model.}
#' \item{Comp.common, Comp.random}{A vector of component effects (common
#'   and random effects model).}
#' \item{seComp.common, seComp.random}{A vector with corresponding
#'   standard errors (common and random effects model).}
#' \item{lower.Comp.common, lower.Comp.random}{A vector with lower
#'   confidence limits for components (common and random effects
#'   model).}
#' \item{upper.Comp.common, upper.Comp.random}{A vector with upper
#'   confidence limits for components (common and random effects
#'   model).}
#' \item{statistic.Comp.common, statistic.Comp.random}{A vector with z-values for
#'   the overall effect of components (common and random effects
#'   model).}
#' \item{pval.Comp.common, pval.Comp.random}{A vector with p-values for
#'   the overall effect of components (common and random effects
#'   model).}
#' \item{Comb.common, Comb.random}{A vector of combination effects
#'   (common and random effects model).}
#' \item{seComb.common, seComb.random}{A vector with corresponding
#'   standard errors (common and random effects model).}
#' \item{lower.Comb.common, lower.Comb.random}{A vector with lower
#'   confidence limits for combinations (common and random effects
#'   model).}
#' \item{upper.Comb.common, upper.Comb.random}{A vector with upper
#'   confidence limits for combinations (common and random effects
#'   model).}
#' \item{statistic.Comb.common, statistic.Comb.random}{A vector with
#'   z-values for the overall effect of combinations (common and random
#'   effects model).}
#' \item{pval.Comb.common, pval.Comb.random}{A vector with p-values for
#'   the overall effect of combinations (common and random effects
#'   model).}
#' \item{Q.additive}{Overall heterogeneity / inconsistency statistic
#'   (additive model).}
#' \item{df.Q.additive}{Degrees of freedom for test of heterogeneity /
#'   inconsistency (additive model).}
#' \item{pval.Q.additive}{P-value for test of heterogeneity /
#'   inconsistency (additive model).}
#' \item{tau}{Square-root of between-study variance (additive model).}
#' \item{I2, lower.I2, upper.I2}{I-squared, lower and upper confidence
#'   limits.}
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
#' \item{A.matrix}{Adjacency matrix (\emph{n}x\emph{n}).}
#' \item{B.matrix}{Edge-vertex incidence matrix (\emph{m}x\emph{n}).}
#' \item{C.matrix}{As defined above.}
#' \item{sm}{Summary measure.}
#' \item{level.ma}{Level for confidence intervals.}
#' \item{common, random, tau.preset}{As defined above.}
#' \item{sep.trts}{A character used in comparison names as separator
#'   between treatment labels.}
#' \item{nchar.comps}{A numeric defining the minimum number of
#'   characters used to create unique component names.}
#' \item{inactive, sep.comps}{As defined above.}
#' \item{backtransf}{A logical indicating whether results should be
#'   back transformed in printouts and forest plots.}
#' \item{title}{Title of meta-analysis / systematic review.}
#' \item{x}{As defined above.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create
#'   object.}
#' 
#' @author Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}, Guido
#'   Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{discomb}}, \code{\link{netmeta}},
#'   \code{\link{forest.netcomb}}, \code{\link{print.netcomb}},
#'   \code{\link{netcomplex}}, \code{\link{netcomparison}}
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
#' Rücker G, Petropoulou M, Schwarzer G (2020):
#' Network meta-analysis of multicomponent interventions.
#' \emph{Biometrical Journal},
#' \bold{62}, 808--21
#' 
#' Welton NJ, Caldwell DM, Adamopoulos E, Vedhara K (2009):
#' Mixed treatment comparison meta-analysis of complex interventions:
#' psychological interventions in coronary heart disease.
#' \emph{American Journal of Epidemiology},
#' \bold{169}: 1158--65
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
#' net1
#' forest(net1, xlim = c(0.2, 50))
#' 
#' # Additive model for treatment components (with placebo as inactive
#' # treatment)
#' #
#' nc1 <- netcomb(net1, inactive = "placebo")
#' nc1
#' forest(nc1, xlim = c(0.2, 50))
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
#' net2
#' forest(net2, xlim = c(0.2, 50))
#' 
#' # Additive model for treatment components (with placebo as inactive
#' # treatment)
#' #
#' nc2 <- netcomb(net2, inactive = "placebo")
#' nc2
#' forest(nc2, xlim = c(0.2, 50))
#' }
#' 
#' @export netcomb


netcomb <- function(x,
                    inactive = NULL,
                    sep.comps = "+",
                    C.matrix,
                    common = x$common,
                    random = x$random | !is.null(tau.preset),
                    tau.preset = NULL,
                    details.chkident = FALSE,
                    nchar.comps = x$nchar.trts,
                    ##
                    func.inverse = invmat,
                    overall.hetstat = x$overall.hetstat,
                    backtransf = x$backtransf,
                    ##
                    na.unident = TRUE,
                    warn.deprecated = gs("warn.deprecated"),
                    ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkclass(x, "netmeta")
  ##
  x <- updateversion(x)
  ##
  chkchar(sep.comps, nchar = 0:1, length = 1)
  ##
  if (!is.null(tau.preset))
    chknumeric(tau.preset, min = 0, length = 1)
  ##
  chklogical(details.chkident)
  nchar.comps <- replaceNULL(nchar.comps, 666)
  chknumeric(nchar.comps, min = 1, length = 1)
  chklogical(overall.hetstat)
  chklogical(backtransf)
  chklogical(na.unident)
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
  random <-
    deprecated(random, missing(random), args, "comb.random", warn.deprecated)
  chklogical(random)
  
  
  ##
  ##
  ## (2) Create C.matrix
  ##
  ##
  if (missing(C.matrix)) {
    #
    # Create C-matrix from netmeta object
    #
    if (sep.comps == "")
      C.matrix <-
        createC(x, inactive = inactive,
                sep.comps = "...this_is_not_a_separator...")
    else
      C.matrix <-
        createC(x, inactive = inactive, sep.comps = sep.comps)
    #
    inactive <- attr(C.matrix, "inactive")
    C.matrix <- as.matrix(C.matrix)[x$trts, , drop = FALSE]
  }
  else {
    ##
    ## Check if appropriate C.matrix is provided
    ##
    if (!(is.matrix(C.matrix) | is.data.frame(C.matrix)))
      stop("Argument 'C.matrix' must be a matrix or data frame.",
           call. = FALSE)
    ##
    if (nrow(C.matrix) != x$n)
      stop("Argument 'C.matrix' has wrong number of rows",
           " (must be equal to number of treatments).",
           call. = FALSE)
    ##
    ## Row names must be a permutation of treatments
    ##
    wrong.labels <- FALSE
    ##
    if (is.null(rownames(C.matrix)))
      wrong.labels <- TRUE
    else {
      if (length(unique(x$trts)) ==
          length(unique(tolower(x$trts))) &&
          length(unique(rownames(C.matrix))) ==
          length(unique(tolower(rownames(C.matrix))))
          )
        idx <- charmatch(tolower(rownames(C.matrix)),
                         tolower(x$trts), nomatch = NA)
      else
        idx <- charmatch(rownames(C.matrix), x$trts, nomatch = NA)
      ##
      if (any(is.na(idx)) || any(idx == 0))
        wrong.labels <- TRUE
    }
    ##
    if (wrong.labels)
      stop("Row names of argument 'C.matrix' must be a ",
           "permutation of treatment names:\n  ",
           paste(paste0("'", x$trts, "'"),
                 collapse = " - "),
           call. = FALSE)
    ##
    C.matrix <- C.matrix[x$trts, , drop = FALSE]
    ##
    if (is.data.frame(C.matrix))
      C.matrix <- as.matrix(C.matrix)
  }
  ##
  c <- ncol(C.matrix) # number of components
  
  
  ##
  ##
  ## (3) Create design matrix based on treatment components
  ##
  ##
  X.matrix <- x$B.matrix %*% C.matrix
  colnames(X.matrix) <- colnames(C.matrix)
  rownames(X.matrix) <- x$studlab
  ##
  ## Check for and warn about not uniquely identifiable components
  ##
  comps.unident <- character(0)
  ##
  if (qr(X.matrix)$rank < c) {
    Xplus <- ginv(X.matrix)
    colnames(Xplus) <- rownames(X.matrix)
    rownames(Xplus) <- colnames(X.matrix)
    e <- eigen(Xplus %*% X.matrix)$values
    E <- eigen(Xplus %*% X.matrix)$vectors
    rownames(E) <- rownames(Xplus %*% X.matrix)
    M <- as.matrix(E[, is.zero(e, n = 100)])
    ##
    if (dim(M)[2] > 0) {
      for (i in 1:dim(M)[2])
        comps.unident <-
          c(comps.unident, names(M[, i])[!is.zero(M[, i], n = 100)])
      ##
      comps.unident <- unique(sort(comps.unident))
      warning(paste0("The following component",
                     if (length(comps.unident) > 1)
                       "s are " else " is ",
                     "not uniquely identifiable: ",
                     paste(paste0("'", comps.unident, "'"),
                           collapse = ", ")),
              call. = FALSE)
      ##
      if (details.chkident) {
        M[is.zero(M, n = 100)] <- 0
        prmatrix(M, quote = FALSE, right = TRUE)
      }
    }
  }
  
  
  ##
  ##
  ## (4) Conduct CNMA
  ##
  ##
  df.Q.additive <- x$df.Q + x$n - 1 - qr(X.matrix)$rank
  df.Q.diff <- x$n - 1 - qr(X.matrix)$rank
  ##
  p0 <- prepare(x$TE, x$seTE, x$treat1, x$treat2, x$studlab,
                func.inverse = func.inverse)
  ##
  o <- order(p0$order)
  ##
  ## Common effects models
  ##
  res.c <- nma.additive(p0$TE[o], p0$weights[o], p0$studlab[o],
                        p0$treat1[o], p0$treat2[o], x$level.ma,
                        X.matrix, C.matrix, x$B.matrix,
                        x$Q, df.Q.additive, df.Q.diff,
                        x$n, x$sep.trts)
  ##
  ## Calculate heterogeneity statistics (additive model)
  ##
  Q.additive <- res.c$Q.additive
  ##
  if (!is.null(tau.preset))
    tau <- tau.preset
  else
    tau <- res.c$tau
  ##
  I2 <- res.c$I2
  lower.I2 <- res.c$lower.I2
  upper.I2 <- res.c$upper.I2
  ##
  ## Random effects models
  ##
  p1 <- prepare(x$TE, x$seTE, x$treat1, x$treat2, x$studlab,
                if (is.na(tau)) 0 else tau, invmat)
  ##
  res.r <- nma.additive(p1$TE[o], p1$weights[o], p1$studlab[o],
                        p1$treat1[o], p1$treat2[o], x$level.ma,
                        X.matrix, C.matrix, x$B.matrix,
                        x$Q, df.Q.additive, df.Q.diff,
                        x$n, x$sep.trts)
  ##
  ## Set unidentifiable components and combinations to NA
  ##
  if (na.unident & length(comps.unident) > 0) {
    setNA <- function(x, select) {
      x[names(x) %in% select] <- NA
      x
    }
    ##
    trts <- names(res.c$combinations$TE)
    unident.pattern <- paste0("^", comps.unident, "$", collapse = "|")
    unident.combs <-
      trts[unlist(lapply(lapply(compsplit(trts, sep.comps),
                                grepl,
                                pattern = unident.pattern), any))]
    ##
    res.c$components <- lapply(res.c$components, setNA, comps.unident)
    res.c$combinations <- lapply(res.c$combinations, setNA, unident.combs)
    ##
    res.r$components <- lapply(res.r$components, setNA, comps.unident)
    res.r$combinations <- lapply(res.r$combinations, setNA, unident.combs)
  }
  ##
  if (length(comps.unident) == 0)
    comps.unident <- NULL
  
  
  ##
  ##
  ## (5) Generate CNMA object
  ##
  ##
  res <- list(studlab = x$studlab,
              treat1 = x$treat1,
              treat2 = x$treat2,
              ##
              TE = p0$TE[o],
              seTE = p0$seTE[o],
              seTE.adj = sqrt(1 / p0$weights[o]),
              seTE.adj.common = sqrt(1 / p0$weights[o]),
              seTE.adj.random = sqrt(1 / p1$weights[o]),
              ##
              design = x$design,
              ##
              event1 = x$event1,
              event2 = x$event2,
              n1 = x$n1,
              n2 = x$n2,
              incr = x$incr,
              ##
              k = x$k,
              m = x$m,
              n = x$n,
              d = x$d,
              c = c,
              ##
              trts = x$trts,
              k.trts = x$k.trts,
              n.trts = x$n.trts,
              events.trts = x$events.trts,
              ##
              n.arms = x$n.arms,
              multiarm = x$multiarm,
              ##
              studies = x$studies,
              narms = x$narms,
              ##
              designs = x$designs,
              ##
              comps = names(res.c$components$TE),
              k.comps = NA,
              n.comps = NA,
              events.comps = NA,
              na.unident = na.unident,
              comps.unident = comps.unident,
              ##
              TE.nma.common = x$TE.nma.common,
              seTE.nma.common = x$seTE.nma.common,
              lower.nma.common = x$lower.nma.common,
              upper.nma.common = x$upper.nma.common,
              statistic.nma.common = x$statistic.nma.common,
              pval.nma.common = x$pval.nma.common,
              ##
              TE.cnma.common = res.c$comparisons$TE,
              seTE.cnma.common = res.c$comparisons$seTE,
              lower.cnma.common = res.c$comparisons$lower,
              upper.cnma.common = res.c$comparisons$upper,
              statistic.cnma.common = res.c$comparisons$statistic,
              pval.cnma.common = res.c$comparisons$p,
              ##
              TE.common = res.c$all.comparisons$TE,
              seTE.common = res.c$all.comparisons$seTE,
              lower.common = res.c$all.comparisons$lower,
              upper.common = res.c$all.comparisons$upper,
              statistic.common = res.c$all.comparisons$statistic,
              pval.common = res.c$all.comparisons$p,
              ##
              TE.nma.random = x$TE.nma.random,
              seTE.nma.random = x$seTE.nma.random,
              lower.nma.random = x$lower.nma.random,
              upper.nma.random = x$upper.nma.random,
              statistic.nma.random = x$statistic.nma.random,
              pval.nma.random = x$pval.nma.random,
              ##
              TE.cnma.random = res.r$comparisons$TE,
              seTE.cnma.random = res.r$comparisons$seTE,
              lower.cnma.random = res.r$comparisons$lower,
              upper.cnma.random = res.r$comparisons$upper,
              statistic.cnma.random = res.r$comparisons$statistic,
              pval.cnma.random = res.r$comparisons$p,
              ##
              TE.random = res.r$all.comparisons$TE,
              seTE.random = res.r$all.comparisons$seTE,
              lower.random = res.r$all.comparisons$lower,
              upper.random = res.r$all.comparisons$upper,
              statistic.random = res.r$all.comparisons$statistic,
              pval.random = res.r$all.comparisons$p,
              ##
              Comp.common = unname(res.c$components$TE),
              seComp.common = unname(res.c$components$seTE),
              lower.Comp.common = unname(res.c$components$lower),
              upper.Comp.common = unname(res.c$components$upper),
              statistic.Comp.common = unname(res.c$components$statistic),
              pval.Comp.common = unname(res.c$components$p),
              ##
              Comp.random = unname(res.r$components$TE),
              seComp.random = unname(res.r$components$seTE),
              lower.Comp.random = unname(res.r$components$lower),
              upper.Comp.random = unname(res.r$components$upper),
              statistic.Comp.random = unname(res.r$components$statistic),
              pval.Comp.random = unname(res.r$components$p),
              ##
              Comb.common = unname(res.c$combinations$TE),
              seComb.common = unname(res.c$combinations$seTE),
              lower.Comb.common = unname(res.c$combinations$lower),
              upper.Comb.common = unname(res.c$combinations$upper),
              statistic.Comb.common = unname(res.c$combinations$statistic),
              pval.Comb.common = unname(res.c$combinations$p),
              ##
              Comb.random = unname(res.r$combinations$TE),
              seComb.random = unname(res.r$combinations$seTE),
              lower.Comb.random = unname(res.r$combinations$lower),
              upper.Comb.random = unname(res.r$combinations$upper),
              statistic.Comb.random = unname(res.r$combinations$statistic),
              pval.Comb.random = unname(res.r$combinations$p),
              ##
              Q.additive = Q.additive,
              df.Q.additive = df.Q.additive,
              pval.Q.additive = res.c$pval.Q.additive,
              tau = tau,
              I2 = I2, lower.I2 = lower.I2, upper.I2 = upper.I2,
              ##
              Q.standard = x$Q,
              df.Q.standard = x$df.Q,
              pval.Q.standard = x$pval.Q,
              ##
              Q.diff = res.c$Q.diff,
              df.Q.diff = df.Q.diff,
              pval.Q.diff = res.c$pval.Q.diff, 
              ##
              A.matrix = x$A.matrix,
              X.matrix = X.matrix,
              B.matrix = x$B.matrix,
              C.matrix = C.matrix,
              ##
              L.matrix.common = res.c$L.matrix,
              Lplus.matrix.common = res.c$Lplus.matrix,
              L.matrix.random = res.r$L.matrix,
              Lplus.matrix.random = res.r$Lplus.matrix,
              ##
              H.matrix.common = res.c$H.matrix[o, o],
              H.matrix.random = res.r$H.matrix[o, o],
              ##
              n.matrix = x$n.matrix,
              events.matrix = x$events.matrix,
              ##
              Cov.common = res.c$Cov,
              Cov.random = res.r$Cov,
              ##
              sm = x$sm,
              method = "Inverse",
              level = x$level,
              level.ma = x$level.ma,
              common = x$common,
              random = x$random,
              ##
              reference.group = x$reference.group,
              baseline.reference = x$baseline.reference,
              all.treatments = NULL,
              seq = x$seq,
              ##
              tau.preset = tau.preset,
              ##
              tol.multiarm = x$tol.multiarm,
              details.chkmultiarm = x$details.chkmultiarm,
              ##
              sep.trts = x$sep.trts,
              sep.comps = sep.comps,
              nchar.comps = nchar.comps,
              ##
              inactive = inactive,
              ##
              func.inverse = deparse(substitute(func.inverse)),
              #
              overall.hetstat = overall.hetstat,
              backtransf = backtransf,
              #
              title = x$title,
              ##
              x = x,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  #
  if (!is.null(x$data))
    res$data <- x$data
  ##
  ## Backward compatibility
  ##
  res$fixed <- res$common
  res$comb.fixed <- res$common
  res$comb.random <- res$random
  ##
  res$seTE.adj.fixed <- res$seTE.adj.common
  res$TE.nma.fixed <- res$TE.nma.common
  res$seTE.nma.fixed <- res$seTE.nma.common
  res$lower.nma.fixed <- res$lower.nma.common
  res$upper.nma.fixed <- res$upper.nma.common
  res$statistic.nma.fixed <- res$statistic.nma.common
  res$pval.nma.fixed <- res$pval.nma.common
  res$TE.cnma.fixed <- res$TE.cnma.common
  res$seTE.cnma.fixed <- res$seTE.cnma.common
  res$lower.cnma.fixed <- res$lower.cnma.common
  res$upper.cnma.fixed <- res$upper.cnma.common
  res$statistic.cnma.fixed <- res$statistic.cnma.common
  res$pval.cnma.fixed <- res$pval.cnma.common
  ##
  res$TE.fixed <- res$TE.common
  res$seTE.fixed <- res$seTE.common
  res$lower.fixed <- res$lower.common
  res$upper.fixed <- res$upper.common
  res$statistic.fixed <- res$statistic.common
  res$pval.fixed <- res$pval.common
  ##
  res$Comp.fixed <- res$Comp.common
  res$seComp.fixed <- res$seComp.common
  res$lower.Comp.fixed <- res$lower.Comp.common
  res$upper.Comp.fixed <- res$upper.Comp.common
  res$statistic.Comp.fixed <- res$statistic.Comp.common
  res$pval.Comp.fixed <- res$pval.Comp.common
  ##
  res$Comb.fixed <- res$Comb.common
  res$seComb.fixed <- res$seComb.common
  res$lower.Comb.fixed <- res$lower.Comb.common
  res$upper.Comb.fixed <- res$upper.Comb.common
  res$statistic.Comb.fixed <- res$statistic.Comb.common
  res$pval.Comb.fixed <- res$pval.Comb.common
  ##
  res$L.matrix.fixed <- res$L.matrix.common
  res$Lplus.matrix.fixed <- res$Lplus.matrix.common
  res$H.matrix.fixed <- res$H.matrix.common        
  ##
  class(res) <- "netcomb"
  
  
  res
}
