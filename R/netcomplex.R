#' Calculate effect of arbitrary complex interventions in component
#' network meta-analysis
#' 
#' @description
#' Calculate effect of arbitrary complex interventions (i.e.,
#' combinations of several components) in component network
#' meta-analysis.
#' 
#' @param x An object of class \code{netcomb} or \code{netcomplex}
#'   (print function).
#' @param complex A matrix, vector or single numeric defining the
#'   complex intervention(s) (see Details).
#' @param fixed A logical indicating whether results for fixed effects
#'   / common effects model should be conducted.
#' @param random A logical indicating whether results for random
#'   effects model should be conducted.
#' @param level The level used to calculate confidence intervals for
#'   combinations of components.
#' @param nchar.comps A numeric defining the minimum number of
#'   characters used to create unique names for components (see
#'   Details).
#' @param backtransf A logical indicating whether printed results
#'   should be back transformed. If \code{backtransf=TRUE}, results
#'   for \code{sm="OR"} are printed as odds ratios rather than log
#'   odds ratios.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z-value
#'   of test for overall effect, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for
#'   p-values, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   combination effect should be printed according to JAMA reporting
#'   standards.
#' @param big.mark A character used as thousands separator.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' R functions \code{\link{netcomb}} and \code{\link{discomb}} only
#' report results for complex interventions present in the
#' network. This function can be used to calculate the effect for
#' arbitrary complex interventions.
#'
#' Complex interventions can be specified using argument \code{complex}:
#' \itemize{
#' \item a character vector with definition of complex interventions,
#' \item a single numeric defining the number of components to combine
#'   in a complex intervention,
#' \item a dedicated C matrix.
#' }
#' In order to calculate effects of arbitrary complex interventions, a
#' C matrix is needed which describes how the complex interventions
#' are composed by the components (Rücker et al., 2020, Section
#' 3.2). The C matrix is constructed internally if not provided by
#' argument \code{complex}. All complex interventions occuring in the
#' network are considered if argument \code{complex} is missing.

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
#' R function \code{\link{netcomparison}} can be used to calculate the
#' effect for comparisons of two arbitrary complex intervention in a
#' component network meta-analysis.
#' 
#' @return
#' A list is returned by the function \code{netcomplex} with the
#' following elements:
#' \item{complex}{Complex intervention(s).}
#' \item{Comb.fixed, Comb.random}{A vector of combination effects
#'   (fixed and random effects model).}
#' \item{seComb.fixed, seComb.random}{A vector with corresponding
#'   standard errors (fixed and random effects model).}
#' \item{lower.Comb.fixed, lower.Comb.random}{A vector with lower
#'   confidence limits for combinations (fixed and random effects
#'   model).}
#' \item{upper.Comb.fixed, upper.Comb.random}{A vector with upper
#'   confidence limits for combinations (fixed and random effects
#'   model).}
#' \item{statistic.Comb.fixed, statistic.Comb.random}{A vector with
#'   z-values for the overall effect of combinations (fixed and random
#'   effects model).}
#' \item{pval.Comb.fixed, pval.Comb.random}{A vector with p-values for
#'   the overall effect of combinations (fixed and random effects
#'   model).}
#' \item{fixed, random}{A defined above.}
#' \item{level, nchar.comps, backtransf, x}{A defined above.}
#' \item{C.matrix}{C matrix.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{netcomb}}, \code{\link{discomb}},
#'   \code{\link{netcomparison}}
#' 
#' @references
#' Rücker G, Petropoulou M, Schwarzer G (2020):
#' Network meta-analysis of multicomponent interventions.
#' \emph{Biometrical Journal},
#' \bold{62}, 808--21
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
#'   data = face, ref = "placebo", sm = "OR", fixed = FALSE)
#' 
#' # Additive model for treatment components (with placebo as inactive
#' # treatment)
#' #
#' nc1 <- netcomb(net1, inactive = "placebo")
#' 
#' # Result for combination Face-to-face PST + SSRI
#' netcomplex(nc1, "Face-to-face PST + SSRI", nchar.comps = 4)
#' netcomplex(nc1, "F + S", nchar.comps = 4)
#'
#' # Result for combination Face-to-face PST + SSRI + Placebo
#' netcomplex(nc1, "Face-to-face PST + SSRI + Plac", nchar.comps = 4)
#'
#' \dontrun{
#' # Artificial example
#' t1 <- rep("A", 3)
#' t2 <- c("B+C", "A+C", "C+D")
#' TE <- c(0, 1, 0)
#' seTE <- rep(1, 3)
#' # Conduct (C)NMA
#' net2 <- netmeta(TE, seTE, t1, t2, random = FALSE)
#' nc2 <- netcomb(net2)
#'
#' # Result for combination A + B + C
#' netcomplex(nc2, "A + B + C")
#' # Same results
#' netcomplex(nc2, "A+B+C")
#' netcomplex(nc2, "B+C+A")
#' netcomplex(nc2, "C+B+A")
#' netcomplex(nc2, "c+b+a")
#' 
#' # Generated C matrix
#' netcomplex(nc2, c(LETTERS[1:4], "A+B+C"))$C.matrix
#' 
#' # Results for all possible combinations of two components
#' netcomplex(nc2, 2)
#' # Results for all possible combinations of three components
#' netcomplex(nc2, 3)
#' }
#' 
#' @rdname netcomplex
#' @export
#' @export netcomplex


netcomplex <- function(x, complex,
                       fixed = x$fixed,
                       random = x$random,
                       level = x$level.ma,
                       nchar.comps = x$nchar.trts) {
  
  
  chkclass(x, "netcomb")
  x <- updateversion(x)
  ##
  chklogical(fixed)
  chklogical(random)
  chklevel(level)
  ##
  nchar.comps <- replaceNULL(nchar.comps, 666)
  chknumeric(nchar.comps, min = 1, length = 1)
  

  if (missing(complex)) {
    complex <- x$trts
    complex.orig <- complex
  }
  else
    complex.orig <- complex    
  ##
  n.comps <- length(x$comps)
  n.complex <- length(complex)
  ##
  comps <- x$comps
  add <- rep("", n.complex)
  ##
  if (is.matrix(complex)) {
    C.matrix <- complex
    complex <- rownames(complex)   
  }
  else if (is.numeric(complex)) {
    chknumeric(complex, min = 1, max = n.comps, length = 1)
    ##
    C.matrix <- createC(ncol = n.comps, ncomb = complex)
    colnames(C.matrix) <- nam.comps <- x$comps
    rownames(C.matrix) <- paste0("...", seq_len(nrow(C.matrix)))
    ##
    for (i in seq_len(nrow(C.matrix))) {
      rownames(C.matrix)[i] <-
        paste(colnames(C.matrix[i, C.matrix[i, ] == 1, drop = FALSE]),
              collapse = x$sep.comps)
    }
    ##
    complex <- rownames(C.matrix)
  }
  else if (is.vector(complex)) {
    ##
    complex.list <- compsplit(complex, x$sep.comps)
    ##
    comps.c <- setref(unique(unlist(complex.list)), c(comps, x$inactive),
                      error.text = "component names", length = 0)
    ##
    complex.list <-
      lapply(complex.list, setref, c(comps, x$inactive), length = 0)
    ##
    complex.orig <- complex
    ##
    add <- rep("", n.complex)
    ##
    for (i in seq_len(n.complex)) {
      add[i] <-
        if (attr(compsplit(complex[i], x$sep.comps), "withspace")) " " else ""
      complex[i] <-
        paste(complex.list[[i]],
              collapse = paste0(add[i], x$sep.comps, add[i]))
    }
    ##
    ## Generate C matrix
    ##
    C.matrix <- matrix(0, ncol = n.comps, nrow = n.complex)
    colnames(C.matrix) <- x$comps
    rownames(C.matrix) <- complex
    ##
    for (i in seq_len(n.complex))
      C.matrix[i, ] <- 1L * colnames(C.matrix) %in% complex.list[[i]]
  }
  else
    stop("Argument 'complex' must be a single numeric, a matrix, or a vector.",
         call. = FALSE)
  
  
  ##
  ## Calculate estimates for combinations
  ##
  Comb.fixed <- as.vector(C.matrix %*% x$Comp.fixed)
  seComb.fixed <-
    sqrt(diag(C.matrix %*% x$Lplus.matrix.fixed %*% t(C.matrix)))
  ##
  Comb.random <- as.vector(C.matrix %*% x$Comp.random)
  seComb.random <-
    sqrt(diag(C.matrix %*% x$Lplus.matrix.random %*% t(C.matrix)))
  ##
  ci.f <- ci(Comb.fixed, seComb.fixed, level = level)
  ci.r <- ci(Comb.random, seComb.random, level = level)
  
  
  res <- list(complex = complex,
              ##
              Comb.fixed = ci.f$TE,
              seComb.fixed = ci.f$seTE,
              lower.Comb.fixed = ci.f$lower,
              upper.Comb.fixed = ci.f$upper,
              statistic.Comb.fixed = ci.f$statistic,
              pval.Comb.fixed = ci.f$p,
              ##
              Comb.random = ci.r$TE,
              seComb.random = ci.r$seTE,
              lower.Comb.random = ci.r$lower,
              upper.Comb.random = ci.r$upper,
              statistic.Comb.random = ci.r$statistic,
              pval.Comb.random = ci.r$p,
              ##
              fixed = fixed,
              random = random,
              level = level,
              ##
              C.matrix = C.matrix,
              ##
              comps = colnames(C.matrix)[apply(C.matrix, 2, sum) > 0],
              inactive = x$inactive,
              nchar.comps = nchar.comps,
              ##
              backtransf = x$backtransf,
              ##
              x = x,
              ##
              add = add,
              ##
              complex.orig = complex.orig,
              ##
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- c("netcomplex", class(res))
  
  res
}





#' @rdname netcomplex
#' @method print netcomplex
#' @export

print.netcomplex <- function(x,
                             fixed = x$fixed,
                             random = x$random,
                             backtransf = x$backtransf,
                             nchar.comps = x$nchar.comps,
                             ##
                             digits = gs("digits"),
                             digits.stat = gs("digits.stat"),
                             digits.pval = gs("digits.pval"),
                             ##
                             scientific.pval = gs("scientific.pval"),
                                 zero.pval = gs("zero.pval"),
                             JAMA.pval = gs("JAMA.pval"),
                             ##
                             big.mark = gs("big.mark"),
                             ##
                             legend = TRUE,
                             ##
                             ...) {
  
  chkclass(x, "netcomplex")
  ##
  chklogical(fixed)
  chklogical(random)
  chklogical(backtransf)
  ##
  nchar.comps <- replaceNULL(nchar.comps, 666)
  chknumeric(nchar.comps, min = 1, length = 1)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 0, length = 1)
  ##
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  ##
  chklogical(legend)
  
  
  ##
  ## Abbreviated component labels
  ##
  n.complex <- length(x$complex)
  complex <- rep("", n.complex)
  ##
  if (fixed | random) {
    comps <- c(x$comps, x$inactive)
    comps.abbr <- treats(comps, nchar.comps)
    ##
    for (i in seq_len(n.complex))
      complex[i] <- compos(x$complex[i], comps, comps.abbr,
                           x$x$sep.comps, x$add[i] == " ")
  }
  
  
  relative <- is.relative.effect(x$x$sm)
  ##
  sm.lab <- x$x$sm
  ##
  if (!backtransf & relative)
    sm.lab <- paste0("log", sm.lab)
  ##  
  ci.lab <- paste0(round(100 * x$level, 1), "%-CI")
  
  
  if (fixed) {
    Comb.fixed <- x$Comb.fixed
    lower.Comb.fixed <- x$lower.Comb.fixed
    upper.Comb.fixed <- x$upper.Comb.fixed
    ##
    if (backtransf & relative) {
      Comb.fixed <- exp(Comb.fixed)
      lower.Comb.fixed <- exp(lower.Comb.fixed)
      upper.Comb.fixed <- exp(upper.Comb.fixed)
    }
  }
  ##
  if (random) {
    Comb.random <- x$Comb.random
    lower.Comb.random <- x$lower.Comb.random
    upper.Comb.random <- x$upper.Comb.random
    ##
    if (backtransf & relative) {
      Comb.random <- exp(Comb.random)
      lower.Comb.random <- exp(lower.Comb.random)
      upper.Comb.random <- exp(upper.Comb.random)
    }
  }
  
  
  if (fixed | random) {
    ##
    if (is.numeric(x$complex.orig) & !is.matrix(x$complex.orig))
      if (all(x$complex.orig == 1))
        msg <- "Results for single components "
      else {
        if (x$complex.orig <= 10)
          complex.orig <- c("two", "three", "four",
                            "five", "six", "seven", "eight", "nine",
                            "ten")[match(x$complex.orig, 2:10)]
        else
          complex.orig <- x$complex.orig
        ##
        msg <- paste("Results for all possible combinations of",
                     complex.orig, "components\n")
      }
    else
      msg <- "Results for combinations "
  }
  
  
  if (fixed) {
    pval.f <- formatPT(x$pval.Comb.fixed, digits = digits.pval,
                       scientific = scientific.pval,
                       zero = zero.pval, JAMA = JAMA.pval,
                       lab.NA = "")
    ##
    res.f <-
      cbind(Comb = formatN(Comb.fixed, digits = digits,
                           "NA", big.mark = big.mark),
            CI = formatCI(formatN(round(lower.Comb.fixed, digits),
                                  digits, "NA",
                                  big.mark = big.mark),
                          formatN(round(upper.Comb.fixed, digits),
                                  digits, "NA",
                                  big.mark = big.mark)),
            zval = formatN(x$statistic.Comb.fixed, digits = digits.stat,
                           "NA", big.mark = big.mark),
            pval = pval.f)
    ##
    cat(msg)
    ##
    dimnames(res.f)[[2]] <- c(sm.lab, ci.lab, "z", "p-value")
    ##
    cat("(additive CNMA model, fixed effects model):\n")
    prmatrix(res.f, quote = FALSE, right = TRUE, na.print = "--")
    ##
    if (random)
      cat("\n")
  }
  
  
  if (random) {
    ##
    pval.r <- formatPT(x$pval.Comb.random, digits = digits.pval,
                       scientific = scientific.pval,
                       zero = zero.pval, JAMA = JAMA.pval,
                       lab.NA = "")
    ##
    res.r <-
      cbind(Comb = formatN(Comb.random, digits = digits,
                           "NA", big.mark = big.mark),
            CI = formatCI(formatN(round(lower.Comb.random, digits),
                                  digits, "NA",
                                  big.mark = big.mark),
                          formatN(round(upper.Comb.random, digits),
                                  digits, "NA",
                                  big.mark = big.mark)),
            zval = formatN(x$statistic.Comb.random, digits = digits.stat,
                           "NA", big.mark = big.mark),
            pval = pval.r)
    ##
    dimnames(res.r)[[2]] <- c(sm.lab, ci.lab, "z", "p-value")
    ##
    cat(msg)
    ##
    cat("(additive CNMA model, random effects model):\n")
    prmatrix(res.r, quote = FALSE, right = TRUE, na.print = "--")
  }
  
  
  if (legend && (fixed | random)) {
    diff.comps <- comps != comps.abbr
    any.comps <- any()
    ##
    if (any(diff.comps)) {
      cat("\nLegend:\n")
      ##
      tmat <- data.frame(comps.abbr, comps)
      tmat <- tmat[diff.comps, ]
      names(tmat) <- c("Abbreviation", " Component name")
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(comps.abbr))) 
    }
  }
  
  
  invisible(NULL)
}
