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
#'   or log hazard ratio). Or an R object created with
#'   \code{\link[meta]{pairwise}}.
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
#'   component (see Details).
#' @param sep.comps A single character to define separator between
#'   treatment components.
#' @param C.matrix C matrix (see Details).
#' @param sm A character string indicating underlying summary measure,
#'   e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
#'   \code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param level The level used to calculate confidence intervals for
#'   individual comparisons.
#' @param level.ma The level used to calculate confidence intervals
#'   for network estimates.
#' @param common A logical indicating whether a common effects /
#'   common effects network meta-analysis should be conducted.
#' @param random A logical indicating whether a random effects network
#'   meta-analysis should be conducted.
#' @param reference.group Reference treatment (first treatment is used
#'   if argument is missing).
#' @param baseline.reference A logical indicating whether results
#'   should be expressed as comparisons of other treatments versus the
#'   reference treatment (default) or vice versa. This argument is
#'   only considered if \code{reference.group} has been specified.
#' @param seq A character or numerical vector specifying the sequence
#'   of treatments in printouts.
#' @param tau.preset An optional value for the square-root of the
#'   between-study variance \eqn{\tau^2}.
#' @param tol.multiarm A numeric for the tolerance for consistency of
#'   treatment estimates in multi-arm studies which are consistent by
#'   design.
#' @param tol.multiarm.se A numeric for the tolerance for consistency
#'   of standard errors in multi-arm studies which are consistent by
#'   design. This check is not conducted if the argument is
#'   \code{NULL}.
#' @param details.chkmultiarm A logical indicating whether treatment
#'   estimates and / or variances of multi-arm studies with
#'   inconsistent results or negative multi-arm variances should be
#'   printed.
#' @param details.chkident A logical indicating whether details on
#'   unidentifiable components should be printed.
#' @param sep.trts A character used in comparison names as separator
#'   between treatment labels.
#' @param overall.hetstat A logical indicating whether to print heterogeneity
#'   measures.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and forest plots. If
#'   \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#'   presented as odds ratios rather than log odds ratios, for
#'   example.
#' @param nchar.comps A numeric defining the minimum number of
#'   characters used to create unique names for components (see
#'   Details).
#' @param sep.ia A single character to define separator for interactions.
#' @param func.inverse R function used to calculate the pseudoinverse
#'   of the Laplacian matrix L (see \code{\link{netmeta}}).
#' @param n1 Number of observations in first treatment group.
#' @param n2 Number of observations in second treatment group.
#' @param event1 Number of events in first treatment group.
#' @param event2 Number of events in second treatment group.
#' @param incr Numerical value added to cell frequencies.
#' @param na.unident A logical indicating whether unidentifiable
#'   components and combinations should be set to missing values.
#' @param keepdata A logical indicating whether original data(set)
#'   should be kept in netmeta object.
#' @param title Title of meta-analysis / systematic review.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if studies are excluded from meta-analysis due to zero
#'   standard errors).
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param nchar.trts Deprecated argument (replaced by
#'   \code{nchar.comps}).
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
#' This R function can be used for disconnected networks. Use
#' \code{\link{netmeta}} and \code{\link{netcomb}} for connected
#' networks.
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
#' with
#' \describe{
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
#' An object of classes \code{discomb} and \code{netcomb} with
#' corresponding \code{print}, \code{summary}, and \code{forest}
#' functions. The object is a list containing the following
#' components:
#' \item{studlab}{Study labels.}
#' \item{treat1}{Label/Number for first treatment.}
#' \item{treat2}{Label/Number for second treatment.}
#' \item{TE}{Estimate of treatment effect, i.e. difference between
#'   first and second treatment.}
#' \item{seTE}{Standard error of treatment estimate.}
#' \item{seTE.adj.common, seTE.adj.random}{Standard error of treatment
#'   estimate, adjusted for multi-arm studies.}
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
#' \item{comps}{Unique list of components present in the network.}
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
#' \item{statistic.cnma.common, statistic.cnma.random}{A vector of
#'   length \emph{m} of z-values for the test of an overall effect
#'   estimated by the additive (common and random effects) model.}
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
#' \item{statistic.Comp.common, statistic.Comp.random}{A vector with
#'   z-values for the overall effect of components (common and random
#'   effects model).}
#' \item{pval.Comp.common, pval.Comp.random}{A vector with p-values for
#'   the overall effect of components (common and random effects
#'   model).}
#' \item{Comb.common, Comb.random}{A vector of combination effects (common
#'   and random effects model).}
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
#' \item{X.matrix}{Design matrix (\emph{m}x\emph{n}).}
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
#' @seealso \code{\link{netcomb}}, \code{\link{forest.netcomb}},
#'   \code{\link{summary.netcomb}}, \code{\link{netmeta}},
#'   \code{\link{netconnection}}, \code{\link{netcomplex}},
#'   \code{\link{netcomparison}}
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
#'   stringsAsFactors = FALSE)
#' #
#' trts <- c("A", "A + B", "A + C", "A + D",
#'   "B", "B + C", "C", "D", "E", "F")
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
#'   c(1, 1, 0, 0, 0, 0),  # A + B
#'   c(1, 0, 1, 0, 0, 0),  # A + C
#'   c(1, 0, 0, 1, 0, 0),  # A + D
#'   c(0, 1, 0, 0, 0, 0),  # B
#'   c(0, 1, 1, 0, 0, 0),  # B + C
#'   c(0, 0, 1, 0, 0, 0),  # C
#'   c(0, 0, 0, 1, 0, 0),  # D
#'   c(0, 0, 0, 0, 1, 0),  # E
#'   c(0, 0, 0, 0, 0, 1))  # F
#' #                  
#' colnames(C) <- comps
#' rownames(C) <- trts
#' #
#' dc2 <- discomb(mean, se.mean, t1, t2, study, seq = trts,
#'   C.matrix = C)
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
                    sep.comps = gs("sep.comps"),
                    C.matrix,
                    ##
                    sm,
                    level = gs("level"),
                    level.ma = gs("level.ma"),
                    common = gs("common"),
                    random = gs("random") | !is.null(tau.preset),
                    ##
                    reference.group,
                    baseline.reference = gs("baseline.reference"),
                    seq = gs("sep"),
                    ##
                    tau.preset = NULL,
                    ##
                    tol.multiarm = gs("tol.multiarm"),
                    tol.multiarm.se = gs("tol.multiarm.se"),
                    details.chkmultiarm = gs("details.chkmultiarm"),
                    ##
                    details.chkident = FALSE,
                    ##
                    sep.trts = gs("sep.trts"),
                    nchar.comps = gs("nchar.comps"),
                    #
                    sep.ia = gs("sep.ia"),
                    #
                    func.inverse = invmat,
                    #
                    n1 = NULL,
                    n2 = NULL,
                    event1 = NULL,
                    event2 = NULL,
                    incr = NULL,
                    #
                    overall.hetstat = gs("overall.hetstat"),
                    backtransf = gs("backtransf"),
                    #
                    na.unident = gs("na.unident"),
                    ##
                    title = gs("title"),
                    keepdata = gs("keepdata"),
                    #
                    warn = gs("warn"), warn.deprecated = gs("warn.deprecated"),
                    nchar.trts = nchar.comps,
                    ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkchar(sep.comps, nchar = 1, length = 1)
  ##
  chklevel(level)
  ##
  missing.reference.group <- missing(reference.group) 
  chklogical(baseline.reference)
  ##
  if (!is.null(tau.preset))
    chknumeric(tau.preset, min = 0, length = 1)
  ##
  chknumeric(tol.multiarm, min = 0, length = 1)
  if (!is.null(tol.multiarm.se))
    chknumeric(tol.multiarm.se, min = 0, length = 1)
  chklogical(details.chkmultiarm)
  ##
  chklogical(details.chkident)
  ##
  missing.sep.trts <- missing(sep.trts)
  chkchar(sep.trts)
  #
  overall.hetstat <- replaceNULL(overall.hetstat, TRUE)
  chklogical(overall.hetstat)
  chklogical(backtransf)
  chklogical(na.unident)
  ##
  chkchar(title)
  chklogical(keepdata)
  chklogical(warn)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  level.ma <- deprecated(level.ma, missing(level.ma), args, "level.comb",
                         warn.deprecated)
  chklevel(level.ma)
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
  nchar.comps <-
    deprecated2(nchar.comps, missing(nchar.comps),
                nchar.trts, missing(nchar.trts),
                warn.deprecated)
  nchar.comps <- replaceNULL(nchar.comps, 666)
  chknumeric(nchar.comps, min = 1, length = 1)
  #
  missing.sep.ia <- missing(sep.ia)
  chkchar(sep.ia, nchar = 0:1, length = 1)
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  if (nulldata)
    data <- sfsp
  ##
  ## Catch TE, treat1, treat2, seTE, studlab from data:
  ##
  TE <- catch("TE", mc, data, sfsp)
  ##
  avail.reference.group.pairwise <- FALSE
  ##
  if (is.data.frame(TE) & !is.null(attr(TE, "pairwise"))) {
    is.pairwise <- TRUE
    ##
    sm <- attr(TE, "sm")
    #
    if (missing.reference.group) {
      reference.group <- attr(TE, "reference.group")
      #
      if (is.null(reference.group))
        reference.group <- ""
      else
        avail.reference.group.pairwise <- TRUE
    }
    #
    keep.all.comparisons <- attr(TE, "keep.all.comparisons")
    if (!is.null(keep.all.comparisons) && !keep.all.comparisons)
      stop("First argument is a pairwise object created with ",
           "'keep.all.comparisons = FALSE'.",
           call. = TRUE)
    #
    if (is.null(attr(TE, "varnames")))
      seTE <- TE$seTE
    else
      seTE <- TE[[attr(TE, "varnames")[2]]]
    #
    treat1 <- TE$treat1
    treat2 <- TE$treat2
    studlab <- TE$studlab
    #
    if (!is.null(TE$n1))
      n1 <- TE$n1
    if (!is.null(TE$n2))
      n2 <- TE$n2
    if (!is.null(TE$event1))
      event1 <- TE$event1
    if (!is.null(TE$event2))
      event2 <- TE$event2
    if (!is.null(TE$incr))
      incr <- TE$incr
    #
    pairdata <- TE
    data <- TE
    #
    if (is.null(attr(TE, "varnames")))
      TE <- TE$TE
    else
      TE <- TE[[attr(TE, "varnames")[1]]]
  }
  else {
    is.pairwise <- FALSE
    if (missing(sm))
      if (!is.null(data) && !is.null(attr(data, "sm")))
        sm <- attr(data, "sm")
      else
        sm <- ""
    ##
    seTE <- catch("seTE", mc, data, sfsp)
    #
    treat1 <- catch("treat1", mc, data, sfsp)
    treat2 <- catch("treat2", mc, data, sfsp)
    #
    studlab <- catch("studlab", mc, data, sfsp)
    #
    if (length(studlab) == 0) {
      if (warn)
        warning("No information given for argument 'studlab'. ",
                "Assuming that comparisons are from independent studies.",
                call. = FALSE)
      studlab <- seq(along = TE)
    }
    studlab <- as.character(studlab)
    #
    n1 <- catch("n1", mc, data, sfsp)
    n2 <- catch("n2", mc, data, sfsp)
    #
    event1 <- catch("event1", mc, data, sfsp)
    event2 <- catch("event2", mc, data, sfsp)
    #
    incr <- catch("incr", mc, data, sfsp)
  }
  #
  subset <- catch("subset", mc, data, sfsp)
  missing.subset <- is.null(subset)
  #
  studlab <- as.character(studlab)
  #
  chknumeric(TE)
  chknumeric(seTE)
  ##
  if (!any(!is.na(TE) & !is.na(seTE)))
    stop("Missing data for estimates (argument 'TE') and ",
         "standard errors (argument 'seTE') in all studies.\n  ",
         "No network meta-analysis possible.",
         call. = FALSE)
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
  ## Remove leading and trailing whitespace
  ##
  treat1 <- rmSpace(rmSpace(treat1, end = TRUE))
  treat2 <- rmSpace(rmSpace(treat2, end = TRUE))
  #
  #
  if (!is.null(event1) & !is.null(event2))
    available.events <- TRUE
  else
    available.events <- FALSE
  #
  if (!is.null(n1) & !is.null(n2))
    available.n <- TRUE
  else
    available.n <- FALSE
  #
  if (available.events & is.null(incr))
    incr <- rep(0, length(event2))
  #
  # Store complete dataset in list object data
  # (if argument keepdata is TRUE)
  #
  if (keepdata) {
    if (nulldata & !is.pairwise)
      data <- data.frame(.studlab = studlab, stringsAsFactors = FALSE)
    else if (nulldata & is.pairwise) {
      data <- pairdata
      data$.studlab <- studlab
    }
    else
      data$.studlab <- studlab
    #
    data$.order <- seq_along(studlab)
    #
    data$.treat1 <- treat1
    data$.treat2 <- treat2
    #
    data$.TE <- TE
    data$.seTE <- seTE
    #
    data$.event1 <- event1
    data$.n1 <- n1
    data$.event2 <- event2
    data$.n2 <- n2
    data$.incr <- incr
    #
    # Check for correct treatment order within comparison
    #
    wo <- data$.treat1 > data$.treat2
    #
    if (any(wo)) {
      data$.TE[wo] <- -data$.TE[wo]
      ttreat1 <- data$.treat1
      data$.treat1[wo] <- data$.treat2[wo]
      data$.treat2[wo] <- ttreat1[wo]
      #
      if (isCol(data, ".n1") & isCol(data, ".n2")) {
        tn1 <- data$.n1
        data$.n1[wo] <- data$.n2[wo]
        data$.n2[wo] <- tn1[wo]
      }
      #
      if (isCol(data, ".event1") & isCol(data, ".event2")) {
        tevent1 <- data$.event1
        data$.event1[wo] <- data$.event2[wo]
        data$.event2[wo] <- tevent1[wo]
      }
      #
      #if (isCol(data, ".mean1") & isCol(data, ".mean2")) {
      #  tmean1 <- data$.mean1
      #  data$.mean1[wo] <- data$.mean2[wo]
      #  data$.mean2[wo] <- tmean1[wo]
      #}
      #
      if (isCol(data, ".sd1") & isCol(data, ".sd2")) {
        tsd1 <- data$.sd1
        data$.sd1[wo] <- data$.sd2[wo]
        data$.sd2[wo] <- tsd1[wo]
      }
      #
      if (isCol(data, ".time1") & isCol(data, ".time2")) {
        ttime1 <- data$.time1
        data$.time1[wo] <- data$.time2[wo]
        data$.time2[wo] <- ttime1[wo]
      }
    }
    #
    if (!missing.subset) {
      if (length(subset) == dim(data)[1])
        data$.subset <- subset
      else {
        data$.subset <- FALSE
        data$.subset[subset] <- TRUE
      }
    }
  }
  
  
  ##
  ##
  ## (3) Use subset for analysis
  ##
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
    #
    if (!is.null(n1))
      n1 <- n1[subset]
    if (!is.null(n2))
      n2 <- n2[subset]
    if (!is.null(event1))
      event1 <- event1[subset]
    if (!is.null(event2))
      event2 <- event2[subset]
    if (!is.null(incr))
      incr <- incr[subset]
  }
  ##
  trts <- sort(unique(c(treat1, treat2)))
  ##
  if (compmatch(trts, sep.trts)) {
    if (!missing.sep.trts)
      warning("Separator '", sep.trts, "' used in at least ",
              "one treatment label. Try to use predefined separators: ",
              "':', '-', '_', '/', '+', '.', '|', '*'.",
              call. = FALSE)
    ##
    if (!compmatch(trts, ":"))
      sep.trts <- ":"
    else if (!compmatch(trts, "-"))
      sep.trts <- "-"
    else if (!compmatch(trts, "_"))
      sep.trts <- "_"
    else if (!compmatch(trts, "/"))
      sep.trts <- "/"
    else if (!compmatch(trts, "+"))
      sep.trts <- "+"
    else if (!compmatch(trts, "."))
      sep.trts <- "-"
    else if (!compmatch(trts, "|"))
      sep.trts <- "|"
    else if (!compmatch(trts, "*"))
      sep.trts <- "*"
    else
      stop("All predefined separators ",
           "(':', '-', '_', '/', '+', '.', '|', '*') ",
           "are used in at least one treatment label.",
           "\n   Please specify a different character that ",
           "should be used as separator (argument 'sep.trts').",
           call. = FALSE)
  }
  ##
  if (!is.null(seq))
    seq <- setseq(seq, trts)
  else {
    seq <- trts
    if (is.numeric(seq))
      seq <- as.character(seq)
  }
  
  
  ##
  ##
  ## (4) Additional checks
  ##
  ##
  if (!(any(grepl(sep.comps, treat1, fixed = TRUE)) |
        any(grepl(sep.comps, treat2, fixed = TRUE))))
    warning("No treatment contains the component separator '", sep.comps, "'.",
            call. = FALSE)
  ##
  if (any(treat1 == treat2))
    stop("Treatments must be different (arguments 'treat1' and 'treat2').",
         call. = FALSE)
  ##
  if (length(studlab) != 0)
    studlab <- as.character(studlab)
  else {
    if (warn)
      warning("No information given for argument 'studlab'. ",
              "Assuming that comparisons are from independent studies.",
              call. = FALSE)
    studlab <- seq(along = TE)
  }
  ##
  ## Check for correct number of comparisons
  ##
  tabnarms <- table(studlab)
  sel.narms <- !is_wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  ##
  if (sum(sel.narms) == 1)
    stop(paste0("Study '", names(tabnarms)[sel.narms],
                "' has a wrong number of comparisons.",
                "\n  Please provide data for all treatment comparisons",
                " (two-arm: 1; three-arm: 3; four-arm: 6, ...)."))
  if (sum(sel.narms) > 1)
    stop(paste0("The following studies have a wrong number of comparisons: ",
                paste(paste0("'", names(tabnarms)[sel.narms], "'"),
                      collapse = ", "),
                "\n  Please provide data for all treatment comparisons",
                " (two-arm: 1; three-arm: 3; four-arm: 6, ...)."))
  ##
  ## Check NAs and zero standard errors
  ##
  excl <- is.na(TE) | is.na(seTE) | seTE <= 0
  ##
  if (any(excl)) {
    if (keepdata)
      data$.excl <- excl
    #
    dat.NAs <- data.frame(studlab = studlab[excl],
                          treat1 = treat1[excl],
                          treat2 = treat2[excl],
                          TE = format(round(TE[excl], 4)),
                          seTE = format(round(seTE[excl], 4))
                          )
    warning("Comparison",
            if (sum(excl) > 1) "s",
            " with missing TE / seTE or zero seTE not considered in",
            " network meta-analysis.",
            call. = FALSE)
    cat("Comparison", if (sum(excl) > 1) "s",
        " not considered in network meta-analysis:\n",
        sep = "")
    prmatrix(dat.NAs, quote = FALSE, right = TRUE,
             rowlab = rep("", sum(excl)))
    cat("\n")
    ##
    studlab <- studlab[!excl]
    treat1  <- treat1[!excl]
    treat2  <- treat2[!excl]
    TE      <- TE[!excl]
    seTE    <- seTE[!excl]
    #
    if (!is.null(n1))
      n1 <- n1[!excl]
    if (!is.null(n2))
      n2 <- n2[!excl]
    if (!is.null(event1))
      event1 <- event1[!excl]
    if (!is.null(event2))
      event2 <- event2[!excl]
    if (!is.null(incr))
      incr <- incr[!excl]
    #
    seq <- seq[seq %in% unique(c(treat1, treat2))]
    trts <- trts[trts %in% unique(c(treat1, treat2))]
  }
  ##
  ## Check for correct number of comparisons (after removing
  ## comparisons with missing data)
  ##
  tabnarms <- table(studlab)
  sel.narms <- !is_wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  ##
  if (sum(sel.narms) == 1)
    stop(paste0("After removing comparisons with missing treatment effects",
                " or standard errors,\n  study '",
                names(tabnarms)[sel.narms],
                "' has a wrong number of comparisons.",
                " Please check data and\n  consider to remove study",
                " from network meta-analysis."))
  if (sum(sel.narms) > 1)
    stop(paste0("After removing comparisons with missing treatment effects",
                " or standard errors,\n  the following studies have",
                " a wrong number of comparisons: ",
                paste(paste0("'", names(tabnarms)[sel.narms], "'"),
                      collapse = ", "),
                "\n  Please check data and consider to remove studies",
                " from network meta-analysis."))
  ##
  ## Check for correct treatment order within comparison
  ##
  wo <- treat1 > treat2
  ##
  if (any(wo)) {
    TE[wo] <- -TE[wo]
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
    #
    if (available.n) {
      tn1 <- n1
      n1[wo] <- n2[wo]
      n2[wo] <- tn1[wo]
    }
    #
    if (available.events) {
      tevent1 <- event1
      event1[wo] <- event2[wo]
      event2[wo] <- tevent1[wo]
    }
  }
  #
  # Set reference group
  #
  if (missing.reference.group & !avail.reference.group.pairwise) {
    go.on <- TRUE
    i <- 0
    while (go.on) {
      i <- i + 1
      sel.i <-
        !is.na(TE) & !is.na(seTE) &
        (treat1 == trts[i] | treat2 == trts[i])
      if (sum(sel.i) > 0) {
        go.on <- FALSE
        reference.group <- trts[i]
      }
      else if (i == length(trts)) {
        go.on <- FALSE
        reference.group <- ""
      }
    }
  }
  #
  # Check reference group
  #
  if (reference.group != "")
    reference.group <- setref(reference.group, trts)
  
  
  ##
  ##
  ## (5) Create C.matrix, B.matrix and Design matrix
  ##
  ##
  netc <- netconnection(treat1, treat2, studlab)
  ##
  if (missing(C.matrix)) {
    C.matrix <- createC(netc, inactive = inactive, sep.comps = sep.comps)
    inactive <- attr(C.matrix, "inactive")
    C.matrix <- C.matrix[trts, , drop = FALSE]
  }
  else {
    ##
    if (!(is.matrix(C.matrix) | is.data.frame(C.matrix))) 
      stop("Argument 'C.matrix' must be a matrix or data frame.", 
           call. = FALSE)
    ##
    wrong.trts <- FALSE
    if (is.null(rownames(C.matrix)))
      wrong.trts <- TRUE
    else {
      if (length(unique(trts)) ==
          length(unique(tolower(trts))) &&
          length(unique(rownames(C.matrix))) ==
          length(unique(tolower(rownames(C.matrix))))
          )
        idx <- charmatch(tolower(rownames(C.matrix)), 
                         tolower(trts), nomatch = NA)
      else
        idx <- charmatch(rownames(C.matrix), trts, nomatch = NA)
      ##
      if (any(is.na(idx)) || any(idx == 0)) 
        wrong.trts <- TRUE
    }
    if (wrong.trts) 
      stop(paste0("Row names of argument 'C.matrix' must be a ", 
                  "permutation of treatment names:\n  ",
                  paste(paste0("'", trts, "'"), collapse = " - ")),
           call. = FALSE)
    ##
    C.matrix <- C.matrix[trts, , drop = FALSE]
  }
  if (is.data.frame(C.matrix))
    C.matrix <- as.matrix(C.matrix)
  ##
  c <- ncol(C.matrix) # number of components
  ##
  ## Create B.matrix
  ##
  p0 <- prepare(TE, seTE, treat1, treat2, studlab,
                func.inverse = func.inverse)
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
  X.matrix <- B.matrix %*% C.matrix
  ##
  colnames(X.matrix) <- colnames(C.matrix)
  rownames(X.matrix) <- studlab
  #
  A.matrix <- diag(diag(t(B.matrix) %*% B.matrix)) - t(B.matrix) %*% B.matrix
  #
  # Identify and warn about not uniquely identifiable components
  #
  if (qr(X.matrix)$rank < c)
    comps.unident <- get_unident(X.matrix, details.chkident)
  else
    comps.unident <- NULL
  
    
  
  ##
  ##
  ## (6) Conduct network meta-analyses
  ##
  ##
  tdata <- data.frame(studies = p0$studlab, narms = p0$narms,
                      order = p0$order,
                      stringsAsFactors = FALSE)
  #
  tdata <- tdata[!duplicated(tdata[, c("studies", "narms")]), , drop = FALSE]
  studies <- tdata$studies[order(tdata$order)]
  narms <- tdata$narms[order(tdata$order)]
  n.a <- sum(narms)  
  ##
  comps <- colnames(C.matrix) # treatment components
  ##
  n <- length(trts)
  m <- length(TE)
  k <- length(unique(studlab))
  ##
  ## Common effects models
  ##
  df.Q.additive <- n.a - k - qr(X.matrix)$rank
  ##
  if (netc$n.subnets == 1) {
    net <- netmeta(TE, seTE, treat1, treat2, studlab,
                   tol.multiarm = tol.multiarm,
                   tol.multiarm.se = tol.multiarm.se,
                   details.chkmultiarm = details.chkmultiarm)
    ##
    Q <- net$Q
    df.Q <- net$df.Q
    pval.Q <- net$pval.Q
    ##
    df.Q.diff <- n - 1 - qr(X.matrix)$rank
  }
  else {
    Q <- df.Q <- pval.Q <- df.Q.diff <- NA
  }
  ##
  res.c <- nma_additive(p0$TE[o], p0$weights[o], p0$studlab[o],
                        p0$treat1[o], p0$treat2[o], level.ma,
                        X.matrix, C.matrix, B.matrix,
                        df.Q.additive, n, sep.trts)
  ##
  ## Random effects models
  ##
  if (!is.null(tau.preset))
    tau <- tau.preset
  else
    tau <- res.c$tau
  ##
  p1 <- prepare(TE, seTE, treat1, treat2, studlab, tau, func.inverse)
  ##
  res.r <- nma_additive(p1$TE[o], p1$weights[o], p1$studlab[o],
                        p1$treat1[o], p1$treat2[o], level.ma,
                        X.matrix, C.matrix, B.matrix,
                        df.Q.additive,
                        n, sep.trts)
  #
  # Difference to standard network meta-analysis model
  #
  Q.diff <- res.c$Q.additive - Q
  if (!is.na(Q.diff) && abs(Q.diff) < .Machine$double.eps^0.75)
    Q.diff <- 0
  #
  if (is.na(df.Q.diff) | df.Q.diff == 0)
    pval.Q.diff <- NA
  else
    pval.Q.diff <- 1 - pchisq(Q.diff, df.Q.diff)
  #
  # Set unidentifiable components and combinations to NA
  #
  if (na.unident & length(comps.unident) > 0) {
    trts <- names(res.c$combinations$TE)
    unident.pattern <- paste0("^", comps.unident, "$", collapse = "|")
    unident.combs <-
      trts[unlist(lapply(lapply(compsplit(trts, sep.comps),
                                grepl,
                                pattern = unident.pattern), any))]
    ##
    res.c$components <- lapply(res.c$components, setNA_vars, comps.unident)
    res.c$combinations <- lapply(res.c$combinations, setNA_vars, unident.combs)
    ##
    res.r$components <- lapply(res.r$components, setNA_vars, comps.unident)
    res.r$combinations <- lapply(res.r$combinations, setNA_vars, unident.combs)
  }
  ##
  if (length(comps.unident) == 0)
    comps.unident <- NULL
  #
  comps <- names(res.c$components$TE)
  #
  if (sep.comps == sep.ia)
    stop("Input for arguments 'sep.comps' and 'sep.ia' must be different.",
         call. = FALSE)
  #
  sep.ia <- setsep(comps, sep.ia, missing = missing.sep.ia)
  
  
  ##
  ##
  ## (7) Generate CNMA object
  ##
  ##
  n.comps <- table(p0$studlab)
  designs <- designs(p0$treat1, p0$treat2, p0$studlab,
                     sep.trts = sep.trts)
  ##  
  NAs <- rep(NA, length(res.c$comparisons$TE))
  ##
  res <- list(studlab = p0$studlab[o],
              treat1 = p0$treat1[o],
              treat2 = p0$treat2[o],
              ##
              TE = p0$TE[o],
              seTE = p0$seTE[o],
              seTE.adj = sqrt(1 / p0$weights[o]),
              seTE.adj.common = sqrt(1 / p0$weights[o]),
              seTE.adj.random = sqrt(1 / p1$weights[o]),
              ##
              design = designs$design[o],
              ##
              event1 = event1,
              event2 = event2,
              n1 = n1,
              n2 = n2,
              incr = incr,
              ##
              k = k,
              m = m,
              n = n,
              d = length(unique(designs$design)),
              c = c,
              s = netc$n.subnets,
              ##
              trts = trts,
              k.trts = NA,
              n.trts = if (available.n) NA else NULL,
              events.trts = if (available.events) NA else NULL,
              ##
              n.arms = NA,
              multiarm = NA,
              ##
              studies = studies,
              narms = narms,
              ##
              designs = unique(sort(designs$design)),
              ##
              comps = comps,
              k.comps = NA,
              n.comps = NA,
              events.comps = NA,
              na.unident = na.unident,
              comps.unident = comps.unident,
              ##
              TE.nma.common = NAs,
              seTE.nma.common = NAs,
              lower.nma.common = NAs,
              upper.nma.common = NAs,
              statistic.nma.common = NAs,
              pval.nma.common = NAs,
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
              TE.nma.random = NAs,
              seTE.nma.random = NAs,
              lower.nma.random = NAs,
              upper.nma.random = NAs,
              statistic.nma.random = NAs,
              pval.nma.random = NAs,
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
              Q.additive = res.c$Q.additive, 
              df.Q.additive = df.Q.additive, 
              pval.Q.additive = res.c$pval.Q.additive,
              tau = tau,
              I2 = res.c$I2,
              lower.I2 = res.c$lower.I2, upper.I2 = res.c$upper.I2,
              ##
              Q.standard = Q,
              df.Q.standard = df.Q,
              pval.Q.standard = pval.Q, 
              ##
              Q.diff = Q.diff,
              df.Q.diff = df.Q.diff,
              pval.Q.diff = pval.Q.diff, 
              ##
              A.matrix = A.matrix,
              X.matrix = X.matrix,
              B.matrix = B.matrix,
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
              n.matrix = NULL,
              events.matrix = NULL,
              #
              Cov.common = res.c$Cov,
              Cov.random = res.r$Cov,
              #
              sm = sm,
              method = "Inverse",
              level = level,
              level.ma = level.ma,
              common = common,
              random = random, 
              ##
              reference.group = reference.group,
              baseline.reference = baseline.reference,
              all.treatments = NULL,
              seq = seq,
              ##
              method.tau = "DL",
              tau.preset = tau.preset,
              #
              tol.multiarm = tol.multiarm,
              details.chkmultiarm = details.chkmultiarm,
              #
              sep.trts = sep.trts,
              sep.comps = sep.comps,
              nchar.comps = nchar.comps,
              sep.ia = sep.ia,
              #
              inactive = inactive,
              #
              func.inverse = deparse(substitute(func.inverse)),
              #
              overall.hetstat = overall.hetstat,
              backtransf = backtransf, 
              ##
              title = title,
              ##
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )
  #
  class(res) <- c("discomb", "netcomb")
  #
  # Add information on multi-arm studies
  #
  if (any(res$narms > 2)) {
    tdata1 <- data.frame(studlab = res$studlab,
                         .order = seq(along = res$studlab))
    tdata2 <- data.frame(studlab = res$studies, narms = res$narms)
    ##
    tdata12 <- merge(tdata1, tdata2,
                     by = "studlab", all.x = TRUE, all.y = FALSE,
                     sort = FALSE)
    tdata12 <- tdata12[order(tdata12$.order), ]
    res$n.arms <- tdata12$narms
    res$multiarm <- tdata12$narms > 2
  }
  else {
    res$n.arms <- rep(2, length(res$studlab))
    res$multiarm <- rep(FALSE, length(res$studlab))
  }
  #
  # Additional assignments
  #
  l1 <- length(res$treat1)
  tab.trts <-
    table(longarm(res$treat1, res$treat2,
                  rep(1, l1), rep(2, l1), rep(1, l1), rep(2, l1),
                  studlab = res$studlab)$treat)
  res$k.trts <- as.numeric(tab.trts)
  names(res$k.trts) <- names(tab.trts)
  #
  if (keepdata) {
    data$.design <- designs(data$.treat1, data$.treat2, data$.studlab,
                            sep = sep.trts)$design
    #
    res$data <- merge(data,
                      data.frame(.studlab = res$studies, .narms = res$narms),
                      by = ".studlab",
                      stringsAsFactors = FALSE)
    #
    # Store adjusted standard errors in dataset
    #
    if (isCol(res$data, ".subset"))
      sel.s <- res$data$.subset
    else
      sel.s <- rep(TRUE, nrow(res$data))
    #
    res$data$.seTE.adj.common <- NA
    res$data$.seTE.adj.random <- NA
    #
    res$data$.seTE.adj.common[sel.s] <- res$seTE.adj.common
    res$data$.seTE.adj.random[sel.s] <- res$seTE.adj.random
    #
    res$data <- res$data[order(res$data$.order), ]
    res$data$.order <- NULL
  }
  #
  # Add information on events and sample sizes
  #
  if (available.n) {
    res$n.matrix <- netmatrix(res, n1 + n2, func = "sum")
    #
    dat.n <- data.frame(studlab = c(studlab, studlab),
                        treat = c(treat1, treat2),
                        n = c(n1, n2))
    dat.n <- dat.n[!duplicated(dat.n[, c("studlab", "treat")]), ]
    dat.n <- by(dat.n$n, dat.n$treat, sum, na.rm = TRUE)
    res$n.trts <- as.vector(dat.n[trts])
    names(res$n.trts) <- trts
  }
  #
  if (available.events) {
    res$events.matrix <- netmatrix(res, event1 + event2, func = "sum")
    #
    dat.e <- data.frame(studlab = c(studlab, studlab),
                        treat = c(treat1, treat2),
                        n = c(event1, event2))
    dat.e <- dat.e[!duplicated(dat.e[, c("studlab", "treat")]), ]
    dat.e <- by(dat.e$n, dat.e$treat, sum, na.rm = TRUE)
    res$events.trts <- as.vector(dat.e[trts])
    names(res$events.trts) <- trts
  }
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
  
  
  res
}
