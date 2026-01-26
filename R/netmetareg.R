#' Network meta-regression with a single continuous or binary covariate
#'
#' @description
#' Network meta-regression with a single continuous or binary covariate for
#' objects of class \code{netmeta}. This is a wrapper function for the R
#' function \code{\link[metafor]{rma.mv}} in the R package \bold{metafor}
#' (Viechtbauer 2010).
#'
#' @details
#' This R function is a wrapper function for R function
#' \code{\link[metafor]{rma.mv}} in the R package \bold{metafor}
#' (Viechtbauer 2010).
#'
#' Note, results are not back-transformed in printouts of
#' network meta-analyses using summary measures with transformations, e.g.,
#' log risk ratios are printed instead of the risk ratio if argument
#' \code{sm = "RR"}.
#'
#' Argument '\dots{}' can be used to pass additional arguments to R
#' function \code{\link[metafor]{rma.mv}}. For example, argument
#' \code{control} to provide a list of control values for the
#' iterative estimation algorithm. See help page of R function
#' \code{\link[metafor]{rma.mv}} for more details.
#'
#' @param x An object of class \code{netmeta}.
#' @param covar Continuous or binary covariate.
#' @param consistency A logical indicating whether a consistency or
#'   inconsistency model should be assumed.
#' @param assumption A character string indicating which assumption is done
#'   for the covariate; either "independent" or "common" (can be abbreviated).
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance tau-squared. Either
#'   \code{"FE"}, \code{"REML"}, or \code{"ML"}.
#' @param level The level used to calculate confidence intervals for regression
#'   coefficients.
#' @param reference.group Reference treatment.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard
#'   errors.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-value, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param print.se A logical specifying whether standard errors should be
#'   printed.
#' @param details.methods A logical specifying whether details on statistical
#'   methods should be printed.
#' @param \dots Additional arguments passed to R function
#'   \code{\link[metafor]{rma.uni}}.
#'
#' @return
#' An object of class \code{c("netmetareg", "rma.mv", "rma")}. Please
#' look at the help page of R function \code{\link[metafor]{rma.mv}}
#' for more details on the output from this function.
#'
#' In addition, a list \code{.netmeta} is added to the output containing
#' the following components:
#' \item{x, covar, method.tau}{As defined above.}
#' \item{dots}{Information provided in argument '\dots{}'.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{netmeta} used to create
#'   object.}
#' \item{version.metafor}{Version of R package \bold{metafor} used to
#'   create object.}
#'
#' @author Nana-Adjoa Kwarteng
#'   \email{nana-adjoa.kwarteng@uniklinik-freiburg.de},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{netmeta}}
#'
#' @references
#' Viechtbauer W (2010):
#' Conducting Meta-Analyses in R with the Metafor Package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#'
#' @keywords models regression
#' 
#' @examples
#' \donttest{
#' data(smokingcessation)
#' # Add variable with (fictitious) risk of bias values
#' # with 1 = "low risk" and 2 = "high risk"
#' #
#' smokingcessation$rob <- rep(1:2, 12)
#' 
#' pw1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' 
#' net1 <- netmeta(pw1, common = FALSE, ref = "A")
#' 
#' # Network meta-regression with continuous covariate and assumption of
#' # independent slopes
#' nr1 <- netmetareg(net1, rob)
#' nr1
#' }
#' 
#' @rdname netmetareg
#' @method netmetareg netmeta
#' @export 

netmetareg.netmeta <- function(x, covar = NULL,
                               consistency = TRUE,
                               assumption = "independent",
                               method.tau = if (!x$random) "FE" else "REML",
                               level = x$level.ma,
                               reference.group = x$reference.group,
                               nchar.trts = x$nchar.trts, ...) {
  
  #
  #
  # (1) Checks and assignments
  #
  #
  
  chkclass(x, "netmeta")
  #
  x <- updateversion(x)
  #
  chklogical(consistency)
  assumption <- setchar(assumption, c("independent", "common"))
  #
  reference.group <- setref(reference.group, x$trts)
  method.tau <- setchar(method.tau, c("REML", "ML", "FE"))
  chklevel(level)
  #
  sm <- x$sm
  
  
  #
  #
  # (2) Return network meta-analysis object if covariate is missing
  #
  #

  if (missing(covar)) {
    warning("No network meta-regresssion conducted as argument 'covar'",
            "is missing.")
    return(x)
  }
  
  
  #
  #
  # (3) Create analysis data set
  #
  #
  
  data <- x$data
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  if (is.null(data))
    data <- sfsp
  #
  # Catch covariate
  #
  covar.name <- gsub("\"", "", deparse(substitute(covar)), fixed = TRUE)
  covar.name <- gsub("$", ".", covar.name, fixed = TRUE)
  covar.name <- gsub("+", "_", covar.name, fixed = TRUE)
  covar.name <- gsub("-", "_", covar.name, fixed = TRUE)
  covar.name <- gsub(":", "_", covar.name, fixed = TRUE)
  #
  covar <- catch("covar", mc, data, sfsp)
  #
  trts <- x$trts
  trts.abbr <- treats(x$trts, nchar.trts)
  #
  reference.group <- trts.abbr[trts == reference.group]
  
  dat <- data.frame(studlab = x$data$.studlab,
                    treat1 = x$data$.treat1,
                    treat2 = x$data$.treat2,
                    TE = x$data$.TE,
                    seTE = x$data$.seTE)
  #
  dat$treat1 <-
    as.character(factor(dat$treat1, levels = trts, labels = trts.abbr))
  dat$treat2 <-
    as.character(factor(dat$treat2, levels = trts, labels = trts.abbr))
  #
  dat$treat1 <- gsub(" ", "_", dat$treat1, fixed = TRUE)
  dat$treat2 <- gsub(" ", "_", dat$treat2, fixed = TRUE)
  #
  if (length(unique(c(dat$treat1, dat$treat2))) != length(x$trts))
    stop("Number of treatments changes after replacing whitespaces with '_'.",
         call. = FALSE)
  #
  if (is.factor(covar) && length(levels(covar)) > 2)
    stop("Network meta-regression not available for a covariate of ",
         "class 'factor' with more than two categories.",
         call. = FALSE)
  #
  # if (is.logical(covar)) # deprecated because model matrix is now constructed manually
  #  covar <- as.numeric(covar)
  #
  dat[[covar.name]] <- covar
  #
  available.n <- available.events <- available.times <- available.sds <- FALSE
  if (!is.null(x$data$.n1) & !is.null(x$data$.n2)) {
    dat$n1 <- x$data$.n1
    dat$n2 <- x$data$.n2
    available.n <- TRUE
  }
  #
  if (!is.null(x$data$.mean1) & !is.null(x$data$.mean2)) {
    dat$mean1 <- x$data$.mean1
    dat$mean2 <- x$data$.mean2
  }
  #
  if (!is.null(x$data$.sd1) & !is.null(x$data$.sd2)) {
    dat$sd1 <- x$data$.sd1
    dat$sd2 <- x$data$.sd2
    available.sds <- TRUE
  }
  #
  if (!is.null(x$data$.event1) & !is.null(x$data$.event2)) {
    dat$event1 <- x$data$.event1
    dat$event2 <- x$data$.event2
    available.events <- TRUE
  }
  #
  if (!is.null(x$data$.time1) & !is.null(x$data$.time2)) {
    dat$time1 <- x$data$.time1
    dat$time2 <- x$data$.time2
    available.times <- TRUE
  }
  #
  if (!is.null(x$data$.incr1) & !is.null(x$data$.incr2)) {
    dat$incr1 <- x$data$.incr1
    dat$incr2 <- x$data$.incr2
  }
  #
  dat <- dat[order(dat$studlab, dat$treat1, dat$treat2), , drop = FALSE]
  #
  keep <- logical(0)
  wo <- logical(0)
  
  # Check for constant covariate in reference treatment
  #
  if (consistency) {
    if (length(unique(dat[dat$treat1 == reference.group |
                          dat$treat2 == reference.group, covar.name])) == 1)
      stop("Invalid reference treatment for interaction. Insufficient ",
           "variation in observed covariate values for reference treatment.")
    
    for (i in unique(dat$studlab)) {
      d.i <- dat[dat$studlab == i, , drop = FALSE]
      trts.i <- unique(sort(c(d.i$treat1, d.i$treat2)))
      if (reference.group %in% trts.i)
        ref.i <- reference.group
      else
        ref.i <- rev(trts.i)[1]
      #
      keep.i <- !(d.i$treat1 != ref.i & d.i$treat2 != ref.i)
      wo.i <- d.i$treat1 == ref.i
      #
      keep <- c(keep, keep.i)
      wo <- c(wo, wo.i)
    }
  }
  dat <- dat[keep, , drop = FALSE]
  #
  wo <- wo[keep]
  if (sum(wo) > 0) {
    dat$TE[wo] <- -dat$TE[wo]
    #
    t2.i <- dat$treat2
    e2.i <- dat$event2
    n2.i <- dat$n2
    mean2.i <- dat$mean2
    sd2.i <- dat$sd2
    time2.i <- dat$time2
    incr2.i <- dat$incr2
    #
    dat$treat2[wo] <- dat$treat1[wo]
    dat$event2[wo] <- dat$event1[wo]
    dat$n2[wo] <- dat$n1[wo]
    dat$mean2[wo] <- dat$mean1[wo]
    dat$sd2[wo] <- dat$sd1[wo]
    dat$time2[wo] <- dat$time1[wo]
    dat$incr2[wo] <- dat$incr1[wo]
    #
    dat$treat1[wo] <- t2.i[wo]
    dat$event1[wo] <- e2.i[wo]
    dat$n1[wo] <- n2.i[wo]
    dat$mean1[wo] <- mean2.i[wo]
    dat$sd1[wo] <- sd2.i[wo]
    dat$time1[wo] <- time2.i[wo]
    dat$incr1[wo] <- incr2.i[wo]
  }
  #
  ncols1 <- ncol(dat)
  dat <- contrmat(dat, grp1 = "treat1", grp2 = "treat2")
  dat <- dat[order(dat$studlab), ]
  trts.all <- names(dat)[(ncols1 + 1):ncol(dat)]
  #
  # Ensure that the reference treatment is correctly assigned
  trts <- trts.all[trts.all != reference.group]
  
  
  #
  #
  # (4) Create model matrix
  #
  #
  
  if (consistency) {
    if (assumption == "independent") {
      formula.nmr_default <- # renamed for construction of manual matrix below
        as.formula(paste("~ 0 + ",
                         paste(trts, collapse = " + "),
                         if (!is.null(covar))
                           paste0( " + ",
                             paste(paste0(trts, ":", covar.name),
                                 collapse = " + "))))
    }
    else {
      #
      # The interaction is all non-reference treatments vs reference,
      # so we can define a variable indicating whether it is reference
      # or non-reference.
      # Then get the interaction using the colon(:). The output has a
      # colon which is in the same format as the independent. Helpful
      # for later extraction
      #
      dat$nonref <- as.numeric(dat[[make.names(reference.group)]] != 0)
      #
      formula.nmr_default <- # renamed for construction of manual matrix below
        as.formula(paste("~ 0 + ",
                         paste(trts, collapse = " + "),
                         if (!is.null(covar))
                           paste0( " + ",
                                   paste(paste0("nonref:", covar.name),
                                         collapse = " + "))))
    }
  }
  else {
    warning("Inconsistency models not yet implemented.")
    return(NULL)
  }
  #
  # Checks for non-numeric covariate
  #
  mm <- model.matrix(formula.nmr_default, data = dat) # default model matrix
  #
  # Drop interactions which contain the reference covariate level if covariate
  # is of mode factor, character, or logical
  #
  if (is.factor(covar))
    mm <- mm[, !grepl(paste0(levels(covar)[1], "$"), colnames(mm)),
             drop = FALSE]
  else if (is.character(covar))
    mm <- mm[, !grepl(paste0(min(covar, na.rm = TRUE), "$"), colnames(mm)),
             drop = FALSE]
  else if (is.logical(covar))
    mm <- mm[, !grepl(paste0(as.logical(min(covar)), "$"), colnames(mm)),
             drop = FALSE]
  #
  dat[colnames(mm)] <- mm # update model matrix columns
  #
  # The back ticks are necessary if the covariate levels / values contain
  # spaces or special characters
  #
  formula.nmr <-
    as.formula(paste("~ 0 ",
                     paste0(" + ",
                            paste0("`", colnames(mm), "`", collapse = " + "))))
  
  
  #
  #
  # (5) Run network meta-regression
  #
  #
  
  # Get rid of warning 'Undefined global functions or variables'
  #
  treat1 <- treat2 <-  comparison <- NULL
  #
  # Covariate 'x' makes problems without removing network meta-analysis object x
  #
  ..x <- x
  rm(x)
  #
  dat.TE <- dat$TE
  dat$comparison <- paste(dat$treat1, dat$treat2, sep = " vs ")
  #
  # Calculate Variance-Covariance matrix
  #
  if (available.n &
      (available.events | available.times | (available.sds))) {
    V <- bldiag(lapply(split(dat, dat$studlab), calcV, sm = sm))
  }
  else
    V <- dat$seTE^2
  #
  suppressWarnings(
    res <-
      runNN(rma.mv,
            list(yi = dat.TE, V = V,
                 data = dat,
                 mods = formula.nmr,
                 random = as.call(~ factor(comparison) | studlab),
                 rho = 0.5,
                 method = method.tau,
                 level = 100 * level,
                 ...)))
  #
  X <- res$X.f
  rownames(X) <- dat$studlab
  X <- X[rev(do.call(order, as.data.frame(X))), , drop = FALSE]
  #
  res$.netmeta <- list(x = ..x,
                       covar.name = covar.name,
                       covar = covar,
                       X = X,
                       consistency = consistency,
                       assumption = assumption,
                       method.tau = method.tau,
                       level = level,
                       reference.group = reference.group,
                       trts = ..x$trts,  # keep the long name of the treatment
                       trts.abbr = trts.abbr,
                       nchar.trts = nchar.trts,
                       dots = list(...),
                       call = match.call(),
                       version = packageDescription("netmeta")$Version,
                       version.metafor = packageDescription("metafor")$Version)
  #
  res$results <- nmr_results(res)
  res$full_results <- nmr_full_results(res)
  #
  class(res) <- c("netmetareg", class(res))
  res$call <- NULL
  #
  res
}


#' @rdname netmetareg
#' @export netmetareg

netmetareg <- function(x, ...)
  UseMethod("netmetareg")


#' @rdname netmetareg
#' @method netmetareg default
#' @export

netmetareg.default <- function(x, ...)
  stop("Network meta-regression not available for an object of class '",
       class(x)[1], "'.",
       call. = FALSE)


#' @rdname netmetareg
#' @method print netmetareg
#' @export

print.netmetareg <- function(x,
                             digits = gs("digits"),
                             digits.se = gs("digits.se"),
                             digits.stat = gs("digits.stat"),
                             digits.pval = gs("digits.pval"),
                             print.se = FALSE,
                             details.methods = TRUE, ...) {
  
  chkclass(x, "netmetareg")
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 0, length = 1)
  chklogical(print.se)
  chklogical(details.methods)
  #
  covar <- x$.netmeta$covar
  #
  if (is.null(covar))
    cat("Network Meta-Analysis\n\n")
  else
    cat("Network Meta-Regression\n\n")
  #
  x.orig <- x
  #
  consistency <- x$.netmeta$consistency
  assumption <- x$.netmeta$assumption
  #
  lower <- NULL
  #
  if (is.null(covar) || is.numeric(covar))
    dat <- x$results[ , c("coef", if (print.se) "se", "lower", "upper",
                          "z", "pval")]
  else
    dat <- x$results[ , c("cov_lvl", "cov_ref",
                          "coef", if (print.se) "se", "lower", "upper",
                          "z", "pval")]
  #
  dat$coef <- formatN(dat$coef, digits = digits)
  #
  if (print.se)
    dat$se <- formatN(dat$se, digits = digits.se, ...)
  #
  dat$lower <- formatCI(formatN(dat$lower, digits = digits),
                        formatN(dat$upper, digits = digits))
  dat$upper <- NULL
  #
  dat <- replaceNA(dat, ".")
  #
  names(dat)[names(dat) == "lower"] <-
    paste0(round(100 * x$.netmeta$level, 1), "%-CI")
  #
  dat$z <- formatN(dat$z, digits = digits.stat)
  dat$pval <- formatPT(dat$pval, digits = digits.pval)
  names(dat)[names(dat) == "pval"] <- "p-value"
  #
  prmatrix(dat, quote = FALSE, right = TRUE, ...)
  #
  # Get rid of warning 'Undefined global functions or variables'
  comparison <- NULL
  #
  if (!is.null(covar) & details.methods) {
    cat("\nDetails on network meta-regression methods:\n")
    #
    if (x$.netmeta$method.tau == "FE")
      cat("- Common effects model\n")
    else
      cat("- Random effects model\n")
    #
    if (consistency)
      cat(paste0("- ", if (consistency) "C" else "Inc", "onsistency model\n"))
    #
    if (assumption == "independent")
      cat("- Independent slopes\n")
    else
      cat("- Common slope\n")
    #
    cat(paste0("- Reference group: ", x$.netmeta$reference.group, "\n"))
  }
  else if (details.methods) {
    cat("\nDetails on network meta-analysis methods:\n")
    #
    if (x$.netmeta$method.tau == "FE")
      cat("- Common effects model\n")
    else
      cat("- Random effects model\n")
  }
  #
  invisible(NULL)
}
