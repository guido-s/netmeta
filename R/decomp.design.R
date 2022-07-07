#' Design-based decomposition of Cochran's Q in network meta-analysis
#' 
#' @description
#' This function performs a design-based decomposition of Cochran's Q
#' for assessing the homogeneity in the whole network, the homogeneity
#' within designs, and the homogeneity/consistency between designs. It
#' allows also an assessment of the consistency assumption after
#' detaching the effect of single designs.
#' 
#' @details
#' In the context of network meta-analysis and the assessment of the
#' homogeneity and consistency assumption, a generalized Cochran's Q
#' statistic for multivariate meta-analysis can be used as shown in
#' Krahn et al. (2013). This Q statistic can be decomposed in a sum
#' of within-design Q statistics and one between-designs Q statistic
#' that incorporates the concept of design inconsistency, see Higgins
#' et al. (2012).
#' 
#' For assessing the inconsistency in a random effects model, the
#' between-designs Q statistic can be calculated based on a full
#' design-by-treatment interaction random effects model (see Higgins
#' et al., 2012). This Q statistic will be automatically given in the
#' output (\eqn{\tau^2} estimated by the method of moments (see Jackson
#' et al., 2012). Alternatively, the square-root of the between-study
#' variance can be prespecified by argument \code{tau.preset} to
#' obtain a between-designs Q statistic (in \code{Q.inc.random}), its
#' design-specific contributions \code{Q.inc.design.random.preset}) as
#' well as residuals after detaching of single designs
#' (\code{residuals.inc.detach.random.preset}).
#' 
#' Since an inconsistent treatment effect of one design can
#' simultaneously inflate several residuals, Krahn et al. (2013)
#' suggest for locating the inconsistency in a network to fit a set of
#' extended models allowing for example for a deviating effect of each
#' study design in turn. The recalculated between-designs Q statistics
#' are given in list component \code{Q.inc.detach}. The change of the
#' inconsistency contribution of single designs can be investigated in
#' more detail by a net heat plot (see function
#' \link{netheat}). Designs where only one treatment is involved in
#' other designs of the network or where the removal of corresponding
#' studies would lead to a splitting of the network do not contribute
#' to the inconsistency assessment. These designs are not included in
#' \code{Q.inc.detach}.
#' 
#' @param x An object of class \code{netmeta}.
#' @param tau.preset An optional value for the square-root of the
#'   between-study variance \eqn{\tau^2} (see Details).
#' @param warn A logical indicating whether warnings should be
#'   printed.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' 
#' @return
#' Network meta-analysis with a single design: \code{NULL}. Otherwise,
#' a list containing the following components:
#' \item{Q.decomp}{Data frame with Q statistics (variable \code{Q})
#'   based on the common effects model to assess the
#'   homogeneity/consistency in the whole network, within designs, and
#'   between designs. Corresponding degrees of freedom (\code{df}) and
#'   p-values (\code{p.val}) are also given.}
#' \item{Q.het.design}{Data frame with design-specific decomposition
#'   of the within-designs Q statistic (\code{Q}) of the common effects
#'   model, corresponding degrees of freedom (\code{df}) and p-values
#'   (\code{p.val}) are given.}
#' \item{Q.inc.detach}{Data frame with between-designs Q statistics
#'   (\code{Q}) of the common effects model after detaching of single
#'   designs, corresponding degrees of freedom (\code{df}) and
#'   p-values (\code{p.val}) are given.}
#' \item{Q.inc.design}{A named vector with contributions of single
#'   designs to the between design Q statistic given in
#'   \code{Q.decomp}.}
#' \item{Q.inc.random}{Data frame with between-designs Q statistic
#'   (\code{Q}) based on a random effects model with square-root of
#'   between-study variance \code{tau.within} estimated embedded in a
#'   full design-by-treatment interaction model, corresponding degrees
#'   of freedom (\code{df}) and p-value (\code{p.val}).}
#' \item{Q.inc.random.preset}{Data frame with between-designs Q
#'   statistic (\code{Q}) based on a random effects model with
#'   prespecified square-root of between-study variance
#'   \code{tau.preset} in the case if argument \code{tau.preset} is
#'   not NULL, corresponding degrees of freedom (\code{df}) and
#'   p-value (\code{p.val}).}
#' \item{Q.inc.design.random.preset}{A named vector with contributions
#'   of single designs to the between design Q statistic based on a
#'   random effects model with prespecified square-root of
#'   between-study variance \code{tau.preset} in the case if argument
#'   \code{tau.preset} is given.}
#' \item{residuals.inc.detach}{Matrix with residuals,
#'   i.e. design-specific direct estimates minus the corresponding
#'   network estimates after detaching the design of the column.}
#' \item{residuals.inc.detach.random.preset}{Matrix with residuals
#'   analogous to \code{residuals.inc.detach} but based on a random
#'   effects model with prespecified square-root of between-study
#'   variance \code{tau.preset} in the case if argument
#'   \code{tau.preset} is not NULL.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package netmeta used to create object.}
#' 
#' @author Ulrike Krahn \email{ulrike.krahn@@bayer.com}, Jochem König
#'   \email{koenigjo@@uni-mainz.de}
#' 
#' @seealso \link{netmeta}, \link{netheat}
#' 
#' @references
#' Higgins JPT, Jackson D, Barrett JK, Lu G, Ades AE, White IR (2012):
#' Consistency and inconsistency in network meta-analysis: concepts
#' and models for multi-arm studies.
#' \emph{Research Synthesis Methods},
#' \bold{3}, 98--110
#' 
#' Krahn U, Binder H, König J (2013):
#' A graphical tool for locating inconsistency in network meta-analyses.
#' \emph{BMC Medical Research Methodology},
#' \bold{13}, 35
#' 
#' Jackson D, White IR and Riley RD (2012):
#' Quantifying the impact of between-study heterogeneity in
#' multivariate meta-analyses.
#' \emph{Statistics in Medicine},
#' \bold{31}, 3805--20
#' 
#' @examples
#' data(Senn2013)
#' 
#' # Only consider first five studies (to reduce runtime of example)
#' #
#' studies <- unique(Senn2013$studlab)
#' Senn2013.5 <- subset(Senn2013, studlab %in% studies[1:5])
#' 
#' # Conduct network meta-analysis with placebo as reference treatment
#' #
#' net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
#'   data = Senn2013.5, sm = "MD", reference = "plac")
#' 
#' # Decomposition of Cochran's Q
#' #
#' decomp.design(net1)
#' 
#' @export decomp.design


decomp.design <- function(x, tau.preset = x$tau.preset, warn = TRUE,
                          nchar.trts = x$nchar.trts) {
  
  
  chkclass(x, "netmeta")
  x <- updateversion(x)
  ##
  if (x$n == 2) {
    warning("No decomposition possible for network meta-analysis ",
            "with only two treatments.",
            call. = FALSE)
    return(NULL)
  }
  ##
  if (!is.null(tau.preset))
    chknumeric(tau.preset, min = 0, length = 1)
  chklogical(warn)
  chknumeric(nchar.trts, min = 1, length = 1)
  ##
  if (inherits(x, "netmetabin")) {
    warning("Decomposition of designs not implemented for ",
            "objects created with netmetabin().",
            call. = FALSE)
    ##
    Q.decomp <- data.frame(Q = x$Q.inconsistency,
                           df = x$df.Q.inconsistency,
                           pval = x$pval.Q.inconsistency)
    ##
    Q.decomp$pval[Q.decomp$df == 0] <- NA
    ##
    rownames(Q.decomp) <- "Between designs"
    res <- list(Q.decomp = Q.decomp,
                call = match.call(),
                version = packageDescription("netmeta")$Version
                )
    ##
    attr(res, "netmetabin") <- TRUE
    class(res) <- "decomp.design"
    ##
    return(res)
  }
  
  
  tau.within <- tau.within(x)
  ##
  decomp.random <- decomp.tau(x, tau.preset = tau.within, warn = warn)
  ##
  Q.inc.random <- decomp.random$Q.decomp["Between designs", ]
  Q.inc.design.random <- decomp.random$Q.inc.design
  residuals.detach.random <- decomp.random$residuals.detach
  
  
  if (length(tau.preset) == 1) {
    decomp.random.preset <-
      decomp.tau(x, tau.preset = tau.preset, warn = warn)
    Q.inc.random.preset <-
      decomp.random.preset$Q.decomp["Between designs",]
    Q.inc.design.random.preset <- decomp.random.preset$Q.inc.design
    residuals.inc.detach.random.preset <-
      decomp.random.preset$residuals.inc.detach
  }
  else {
    tau.preset <- NULL
    Q.inc.random.preset <- NULL
    Q.inc.design.random.preset <- NULL
    residuals.inc.detach.random.preset <- NULL
  }
  
  
  dct <- decomp.tau(x, warn = warn)
  
  
  res <- list(
    Q.decomp = dct$Q.decomp,
    Q.het.design = dct$Q.het.design,
    Q.inc.detach = dct$Q.inc.detach,
    Q.inc.design = dct$Q.inc.design,
    ##
    Q.inc.random = data.frame(Q.inc.random, tau.within),
    Q.inc.random.preset = data.frame(Q.inc.random.preset, tau.preset),
    Q.inc.design.random.preset = Q.inc.design.random.preset,
    ##
    residuals.inc.detach = dct$residuals.inc.detach,
    residuals.inc.detach.random.preset = residuals.inc.detach.random.preset,
    ##
    tau.preset = tau.preset,
    nchar.trts = nchar.trts,
    ##
    x = x,
    ##
    call = match.call(),
    version = packageDescription("netmeta")$Version
    )
  ##
  class(res) <- "decomp.design"
  ##
  res
}
