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
#'   (see Details).
#' @param sep.comps A single character to define separator between
#'   treatment components.
#' @param C.matrix C matrix (see Details).
#' @param comb.fixed A logical indicating whether a fixed effects
#'   (common effects) network meta-analysis should be conducted.
#' @param comb.random A logical indicating whether a random effects
#'   network meta-analysis should be conducted.
#' @param tau.preset An optional value for the square-root of the
#'   between-study variance \eqn{\tau^2}.
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
#' the additive model in a frequentist way (Rücker et al., 2019).
#' 
#' The underlying multivariate model is given by
#' 
#' \deqn{\bold{\delta} = \bold{B} \bold{\theta}, \bold{\theta} =
#' \bold{C} \bold{\beta}}
#' 
#' with \itemize{
#' \item[\eqn{\bold{\delta}}] vector of true treatment effects
#'   (differences) from individual studies,
#' \item[\eqn{\bold{B}}] is a design matrix describing the structure
#'   of the network,
#' \item[\eqn{\bold{\theta}}] parameter vector that represents the
#'   existing combined treatments,
#' \item[\eqn{\bold{C}}] matrix describing how the treatments are
#'   composed,
#' \item[\eqn{\bold{\beta}}] is a parameter vector representing the
#' treatment components.
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
#' An object of class \code{netcomb} with corresponding \code{print},
#' \code{summary}, and \code{forest} functions. The object is a list
#' containing the following components:
#' \item{studlab}{Study labels.}
#' \item{treat1}{Label/Number for first treatment.}
#' \item{treat2}{Label/Number for second treatment.}
#' \item{TE}{Estimate of treatment effect, i.e. difference between
#'   first and second treatment.}
#' \item{seTE}{Standard error of treatment estimate.}
#' \item{seTE.adj}{Standard error of treatment estimate, adjusted for
#'   multi-arm studies.}
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
#' \item{TE.nma.fixed, TE.nma.random}{A vector of length \emph{m} of
#'   consistent treatment effects estimated by network meta-analysis
#'   (nma) (fixed and random effects model).}
#' \item{seTE.nma.fixed, seTE.nma.random}{A vector of length \emph{m}
#'   of effective standard errors estimated by network meta-analysis
#'   (fixed and random effects model).}
#' \item{lower.nma.fixed, lower.nma.random}{A vector of length
#'   \emph{m} of lower confidence interval limits for consistent
#'   treatment effects estimated by network meta-analysis (fixed
#'   and random effects model).}
#' \item{upper.nma.fixed, upper.nma.random}{A vector of length
#'   \emph{m} of upper confidence interval limits for the consistent
#'   treatment effects estimated by network meta-analysis (fixed
#'   and random effects model).}
#' \item{zval.nma.fixed, zval.nma.random}{A vector of length \emph{m}
#'   of z-values for test of treatment effect for individual
#'   comparisons (fixed and random effects model).}
#' \item{pval.nma.fixed, pval.nma.random}{A vector of length \emph{m}
#'   of p-values for test of treatment effect for individual
#'   comparisons (fixed and random effects model).}
#' \item{TE.cnma.fixed, TE.cnma.random}{A vector of length \emph{m} of
#'   consistent treatment effects estimated by the additive (fixed and
#'   random effects) model.}
#' \item{seTE.cnma.fixed, seTE.cnma.random}{A vector of length
#'   \emph{m} with standard errors estimated by the additive (fixed
#'   and random effects) model.}
#' \item{lower.cnma.fixed, lower.cnma.random}{A vector of length
#'   \emph{m} of lower confidence interval limits for consistent
#'   treatment effects estimated by the additive (fixed and random
#'   effects) model.}
#' \item{upper.cnma.fixed, upper.cnma.random}{A vector of length
#'   \emph{m} of upper confidence interval limits for consistent
#'   treatment effects estimated by the additive (fixed and random
#'   effects) model.}
#' \item{zval.cnma.fixed, zval.cnma.random}{A vector of length
#'   \emph{m} of z-values for the test of an overall effect estimated
#'   by the additive (fixed and random effects) model.}
#' \item{pval.cnma.fixed, zval.cnma.random}{A vector of length
#'   \emph{m} of p-values for the test of an overall effect estimated
#'   by the additive (fixed and random effects) model.}
#' \item{TE.fixed, TE.random}{\emph{n}x\emph{n} matrix with overall
#'   treatment effects estimated by the additive (fixed and random
#'   effects) model.}
#' \item{seTE.fixed, seTE.random}{\emph{n}x\emph{n} matrix with
#'   standard errors estimated by the additive (fixed and random
#'   effects) model.}
#' \item{lower.fixed, upper.fixed, lower.random,
#'   upper.random}{\emph{n}x\emph{n} matrices with lower and upper
#'   confidence interval limits estimated by the additive (fixed and
#'   random effects) model.}
#' \item{zval.fixed, pval.fixed, zval.random,
#'   pval.random}{\emph{n}x\emph{n} matrices with z-values and
#'   p-values for test of overall effect estimated by the additive
#'   (fixed and random effects) model.}
#' \item{Comp.fixed, Comp.random}{A vector of component effects (fixed
#'   and random effects model).}
#' \item{seComp.fixed, seComp.random}{A vector with corresponding
#'   standard errors (fixed and random effects model).}
#' \item{lower.Comp.fixed, lower.Comp.random}{A vector with lower
#'   confidence limits for components (fixed and random effects
#'   model).}
#' \item{upper.Comp.fixed, upper.Comp.random}{A vector with upper
#'   confidence limits for components (fixed and random effects
#'   model).}
#' \item{zval.Comp.fixed, zval.Comp.random}{A vector with z-values for
#'   the overall effect of components (fixed and random effects
#'   model).}
#' \item{pval.Comp.fixed, pval.Comp.random}{A vector with p-values for
#'   the overall effect of components (fixed and random effects
#'   model).}
#' \item{Comb.fixed, Comb.random}{A vector of combination effects (fixed
#'   and random effects model).}
#' \item{seComb.fixed, seComb.random}{A vector with corresponding
#'   standard errors (fixed and random effects model).}
#' \item{lower.Comb.fixed, lower.Comb.random}{A vector with lower
#'   confidence limits for combinations (fixed and random effects
#'   model).}
#' \item{upper.Comb.fixed, upper.Comb.random}{A vector with upper
#'   confidence limits for combinations (fixed and random effects
#'   model).}
#' \item{zval.Comb.fixed, zval.Comb.random}{A vector with z-values for
#'   the overall effect of combinations (fixed and random effects
#'   model).}
#' \item{pval.Comb.fixed, pval.Comb.random}{A vector with p-values for
#'   the overall effect of combinations (fixed and random effects
#'   model).}
#' \item{Q.additive}{Overall heterogeneity / inconsistency statistic
#'   (additive model).}
#' \item{df.Q.additive}{Degrees of freedom for test of heterogeneity /
#'   inconsistency (additive model).}
#' \item{pval.Q.additive}{P-value for test of heterogeneity /
#'   inconsistency (additive model).}
#' \item{tau}{Square-root of between-study variance (additive model).}
#' \item{I2}{I-squared (additive model).}
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
#' \item{level.comb}{Level for confidence intervals.}
#' \item{comb.fixed, comb.random, tau.preset}{As defined above.}
#' \item{sep.trts}{A character used in comparison names as separator
#'   between treatment labels.}
#' \item{nchar.trts}{A numeric defining the minimum number of
#'   characters used to create unique treatment and component names.}
#' \item{inactive, sep.comps}{As defined above.}
#' \item{backtransf}{A logical indicating whether results should be
#'   back transformed in printouts and forest plots.}
#' \item{title}{Title of meta-analysis / systematic review.}
#' \item{x}{As defined above.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create
#'   object.}
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \link{discomb}, \link{netmeta}, \link{forest.netcomb},
#'   \link{print.netcomb}
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
#' Rücker G, Petropoulou M, Schwarzer G (2019):
#' Network meta-analysis of multicomponent interventions.
#' \emph{Biometrical Journal},
#' 1--14, https://doi.org/10.1002/bimj.201800167
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
#'                 data = face, ref = "placebo",
#'                 sm = "OR", comb.fixed = FALSE)
#' summary(net1)
#' forest(net1, xlim = c(0.2, 50))
#' 
#' # Additive model for treatment components (with placebo as inactive
#' # treatment)
#' #
#' nc1 <- netcomb(net1, inactive = "placebo")
#' summary(nc1)
#' forest(nc1, xlim = c(0.2, 50))
#' 
#' \dontrun{
#' # Specify, order of treatments
#' #
#' trts <- c("TCA", "SSRI", "SNRI", "NRI", "Low-dose SARI", "NaSSa",
#'           "rMAO-A", "Ind drug", "Hypericum", "Face-to-face CBT",
#'           "Face-to-face PST", "Face-to-face interpsy", "Face-to-face psychodyn",
#'           "Other face-to-face", "Remote CBT", "Self-help CBT", "No contact CBT",
#'           "Face-to-face CBT + SSRI", "Face-to-face interpsy + SSRI",
#'           "Face-to-face PST + SSRI", "UC", "Placebo")
#' #
#' # Note, three treatments are actually combinations of 'SSRI' with
#' # other components:
#' # "Face-to-face CBT + SSRI",
#' # "Face-to-face interpsy + SSRI",
#' # "Face-to-face PST + SSRI"
#' 
#' # Conduct random effects network meta-analysis
#' #
#' net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
#'                 data = Linde2016, ref = "placebo",
#'                 seq = trts,
#'                 sm = "OR", comb.fixed = FALSE)
#' summary(net1)
#' forest(net1, xlim = c(0.2, 50))
#' 
#' # Additive model for treatment components (with placebo as inactive
#' # treatment)
#' #
#' nc1 <- netcomb(net1, inactive = "placebo")
#' summary(nc1)
#' forest(nc1, xlim = c(0.2, 50))
#' }
#' 
#' @export netcomb


netcomb <- function(x,
                    inactive = NULL,
                    sep.comps = "+",
                    C.matrix,
                    comb.fixed = x$comb.fixed,
                    comb.random = x$comb.random | !is.null(tau.preset),
                    tau.preset = NULL) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  meta:::chkclass(x, "netmeta")
  ##
  meta:::chkchar(sep.comps, nchar = 1)
  ##
  meta:::chklogical(comb.fixed)
  meta:::chklogical(comb.random)
  ##
  if (!is.null(tau.preset))
    meta:::chknumeric(tau.preset, min = 0, single = TRUE)
  
  
  ##
  ##
  ## (2) Get data
  ##
  ##
  n <- x$n # number of treatments / combinations
  m <- x$m # number of comparisons
  ##
  trts <- x$trts
  ##
  TE <- x$TE
  seTE <- x$seTE
  treat1 <- x$treat1
  treat2 <- x$treat2
  studlab <- x$studlab
  ##
  B.matrix <- x$B.matrix
  
  
  ##
  ##
  ## (3) Create C.matrix
  ##
  ##
  if (missing(C.matrix)) {
    ##
    ## Create C-matrix from netmeta object
    ##
    C.matrix <- as.matrix(createC(x, sep.comps, inactive))
    C.matrix <- C.matrix[trts, , drop = FALSE]
  }
  else {
    ##
    ## Check if appropriate C.matrix is provided
    ##
    if (!(is.matrix(C.matrix) | is.data.frame(C.matrix)))
      stop("Argument 'C.matrix' must be a matrix or data frame.",
           call. = FALSE)
    ##
    if (nrow(C.matrix) != n)
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
      if (length(unique(trts)) == length(unique(tolower(trts))))
        idx <- charmatch(tolower(rownames(C.matrix)),
                         tolower(trts), nomatch = NA)
      else
        idx <- charmatch(rownames(C.matrix), trts, nomatch = NA)
      ##
      if (any(is.na(idx)) || any(idx == 0))
        wrong.labels <- TRUE
    }
    ##
    if (wrong.labels)
      stop(paste("Row names of argument 'C.matrix' must be a ",
                 "permutation of treatment names:\n  ",
                 paste(paste("'", trts, "'", sep = ""),
                       collapse = " - "), sep = ""),
           call. = FALSE)
    ##
    C.matrix <- C.matrix[trts, , drop = FALSE]
    ##
    if (is.data.frame(C.matrix))
      C.matrix <- as.matrix(C.matrix)
  }
  ##
  c <- ncol(C.matrix) # number of components
  
  
  p0 <- prepare(TE, seTE, treat1, treat2, studlab)
  ##
  o <- order(p0$order)
  ##
  B.matrix <- createB(p0$treat1.pos[o], p0$treat2.pos[o])
  ##
  colnames(B.matrix) <- trts
  rownames(B.matrix) <- studlab
  
  
  ##
  ## Design matrix based on treatment components
  ##
  X <- B.matrix %*% C.matrix
  ##
  colnames(X) <- colnames(C.matrix)
  rownames(X) <- studlab
  
  
  ##
  ## Fixed effects models
  ##
  df.Q.additive <- x$df.Q + x$n - 1 - qr(X)$rank
  ##
  net <- netmeta(TE, seTE, treat1, treat2, studlab)
  ##
  Q <- net$Q
  df.Q <- net$df.Q
  pval.Q <- net$pval.Q
  ##
  df.Q.diff <- x$n - 1 - qr(X)$rank
  
  
  res.f <- nma.additive(p0$TE[o], p0$weights[o], p0$studlab[o],
                        p0$treat1[o], p0$treat2[o], x$level.comb,
                        X, C.matrix, B.matrix,
                        Q, df.Q.additive, df.Q.diff,
                        x$n, x$sep.trts)
  
  
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
  tau2.calc <- if (is.na(tau)) 0 else tau^2
  
  
  ##
  ## Random effects models
  ##
  p1 <- prepare(TE, seTE, treat1, treat2, studlab, tau)
  res.r <- nma.additive(p1$TE[o], p1$weights[o], p1$studlab[o],
                        p1$treat1[o], p1$treat2[o], x$level.comb,
                        X, C.matrix, B.matrix,
                        Q, df.Q.additive, df.Q.diff,
                        x$n, x$sep.trts)
  
  
  res <- list(studlab = x$studlab,
              treat1 = x$treat1,
              treat2 = x$treat2,
              ##
              TE = x$TE,
              seTE = x$seTE,
              seTE.adj = x$seTE.adj,
              ##
              event1 = x$event1,
              event2 = x$event2,
              n1 = x$n1,
              n2 = x$n2,
              ##
              k = x$k,
              m = m,
              n = n,
              d = x$d,
              c = c,
              ##
              trts = trts,
              k.trts = x$k.trts,
              n.trts = x$n.trts,
              events.trts = x$events.trts,
              ##
              studies = x$studies,
              narms = x$narms,
              ##
              designs = x$designs,
              ##
              comps = names(res.f$components$TE),
              k.comps = NA,
              n.comps = NA,
              events.comps = NA,
              ##
              TE.nma.fixed = x$TE.nma.fixed,
              seTE.nma.fixed = x$seTE.nma.fixed,
              lower.nma.fixed = x$lower.nma.fixed,
              upper.nma.fixed = x$upper.nma.fixed,
              zval.nma.fixed = x$zval.nma.fixed,
              pval.nma.fixed = x$pval.nma.fixed,
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
              TE.nma.random = x$TE.nma.random,
              seTE.nma.random = x$seTE.nma.random,
              lower.nma.random = x$lower.nma.random,
              upper.nma.random = x$upper.nma.random,
              zval.nma.random = x$zval.nma.random,
              pval.nma.random = x$pval.nma.random,
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
              Q.standard = x$Q,
              df.Q.standard = x$df.Q,
              pval.Q.standard = x$pval.Q,
              ##
              Q.diff = res.f$Q.diff,
              df.Q.diff = df.Q.diff,
              pval.Q.diff = res.f$pval.Q.diff, 
              ##
              A.matrix = x$A.matrix,
              B.matrix = B.matrix,
              C.matrix = C.matrix,
              ##
              n.matrix = x$n.matrix,
              events.matrix = x$events.matrix,
              ##
              sm = x$sm,
              method = "Inverse",
              level = x$level,
              level.comb = x$level.comb,
              comb.fixed = x$comb.fixed,
              comb.random = x$comb.random,
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
              nchar.trts = x$nchar.trts,
              ##
              inactive = inactive,
              sep.comps = sep.comps,
              ##
              backtransf = x$backtransf,
              ##
              title = x$title,
              ##
              x = x,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- "netcomb"
  ##
  res
}
