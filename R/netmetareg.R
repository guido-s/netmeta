#' Network meta-regression
#'
#' @description
#' Network meta-regression for objects of class \code{netmeta}. This is a
#' wrapper function for the R function \code{\link[metafor]{rma.mv}} in the R
#' package \bold{metafor} (Viechtbauer 2010).
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
#' @param covar Covariate.
#' @param consistency A logical indicating whether a consistency or
#'   inconsistency model should be assumed.
#' @param consistency A logical indicating whether a consistency or
#'   inconsistency model should be assumed.
#' @param assumption A character string indicating which assumption is done
#'   for the covariate; either "independent" or "constant" (can be abbreviated).
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance tau-squared. Either
#'   \code{"FE"}, \code{"REML"}, or \code{"ML"},
#' @param reference.group Reference treatment.
#' @param nchar.trts A numeric defining the minimum number of
#'   characters used to create unique treatment names.
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
#' # Add examples
#' 
#' @rdname netmetareg
#' @method netmetareg netmeta
#' @export 

netmetareg.netmeta <- function(x, covar = NULL,
                       consistency = TRUE,
                       assumption = "independent",
                       method.tau = if (!x$random) "FE" else "REML",
                       reference.group = x$reference.group,
                       nchar.trts = x$nchar.trts, ...) {
  
  chkclass(x, "netmeta")
  #
  x <- updateversion(x)
  
  #
  # Checks and assignments
  #
  chklogical(consistency)
  assumption <- setchar(assumption, c("independent", "constant"))
  #
  reference.group <- setref(reference.group, x$trts)
  method.tau <- setchar(method.tau, c("REML", "ML", "FE"))
  #
  sm <- x$sm
  
  #
  # Return network meta-analysis object if covariate is missing
  #
  
  if (missing(covar)) {
    warning("No network meta-regresssion conducted as argument 'covar'",
            "is missing.")
    return(x)
  }
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
  dat[[covar.name]] <- covar
  #
  if (!is.null(x$data$.n1) & !is.null(x$data$.n2)) {
    dat$n1 <- x$data$.n1
    dat$n2 <- x$data$.n2
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
  }
  #
  if (!is.null(x$data$.time1) & !is.null(x$data$.time2)) {
    dat$time1 <- x$data$.time1
    dat$time2 <- x$data$.time2
  }
  #
  keep <- logical(0)
  wo <- logical(0)
  #
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
  #
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
    #
    dat$treat2[wo] <- dat$treat1[wo]
    dat$event2[wo] <- dat$event1[wo]
    dat$n2[wo] <- dat$n1[wo]
    dat$mean2[wo] <- dat$mean1[wo]
    dat$sd2[wo] <- dat$sd1[wo]
    dat$time2[wo] <- dat$time2[wo]
    #
    dat$treat1[wo] <- t2.i[wo]
    dat$event1[wo] <- e2.i[wo]
    dat$n1[wo] <- n2.i[wo]
    dat$mean1[wo] <- mean2.i[wo]
    dat$sd1[wo] <- sd2.i[wo]
    dat$time1[wo] <- time2.i[wo]
  }
  #
  ncols1 <- ncol(dat)
  dat <- contrmat(dat, grp1 = "treat1", grp2 = "treat2")
  ncols2 <- ncol(dat)
  varnames <- names(dat)[(ncols1 + 1):ncols2]
  #
  dat <- dat[order(dat$studlab), ]
  #
  trts.tau.all <- varnames
  trts.tau <- varnames[-length(varnames)]
  #
  formula.trts <-
    as.formula(paste("~ ", paste(trts.tau, collapse = " + "), " - 1"))
  #
  # Calculate Variance-Covariance matrix
  #
  V <- bldiag(lapply(split(dat, dat$studlab), calcV, sm = sm))
  
  
  if (consistency) {
    if (assumption == "independent")
      formula.nmr <-
        as.formula(paste("~ ",
                         paste(trts.tau, collapse = " + "), " + ",
                         paste(paste0(trts.tau, ":", covar.name),
                               collapse = " + "),
                         " - 1"))
    else {
      covar.pre <- dat[[covar.name]]
      #
      if (is.numeric(dat[[covar.name]]))
        dat[[covar.name]] <-
          ifelse(dat[[reference.group]] == 0, 0, dat[[covar.name]])
      #
      formula.nmr <-
        as.formula(paste("~ ",
                         paste(trts.tau, collapse = " + "), " + ",
                         covar.name,
                         " - 1"))
    }
  }
  else {
    warning("Inconsistency models not yet implemented.")
    return(NULL)
  }
  
  
  #
  # Covariate 'x' makes problems without removing network meta-analysis object x
  #
  ..x <- x
  rm(x)
  #
  dat.TE <- dat$TE
  dat$comparison <- paste(dat$treat1, dat$treat2, sep = " vs ")
  #
  res <-
    runNN(rma.mv,
          list(yi = dat.TE, V = V,
               data = dat,
               mods = formula.nmr,
               random = as.call(~ factor(comparison) | studlab),
               rho = 0.5,
               method = method.tau, ...))
  #
  res$.netmeta <- list(x = ..x,
                       covar = covar,
                       consistency = consistency,
                       assumption = assumption,
                       method.tau = method.tau,
                       reference.group = reference.group,
                       trts = trts,
                       trts.abbr = trts.abbr,
                       nchar.trts = nchar.trts,
                       dots = list(...),
                       call = match.call(),
                       version = packageDescription("netmeta")$Version,
                       version.metafor = packageDescription("metafor")$Version)
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


#' @method print netmetareg
#' @export

print.netmetareg <- function(x, details = TRUE, ...) {
  
  chkclass(x, "netmetareg")
  chklogical(details)
  #
  cat("Network Meta-Regression\n")
  #
  x.orig <- x
  #
  consistency <- x$.netmeta$consistency
  assumption <- x$.netmeta$assumption
  #
  class(x) <- class(x)[class(x) != "netmetareg"]
  #
  print(x, ...)
  #
  if (details) {
    cat("Details on network meta-regression methods:\n")
    if (consistency)
      cat(paste0("- ", if (consistency) "C" else "Inc", "onsistency model\n"))
    if (assumption == "independent")
      cat("- Independent slopes\n")
    else
      cat("- Constant slope\n")
  }
  #
  invisible(NULL)
}
