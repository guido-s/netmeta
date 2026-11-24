#' Network meta-regression with a single continuous or binary covariate
#'
#' @description
#' Network meta-regression with a single continuous or binary covariate for
#' objects of class \code{netmeta} (Kwarteng et al., 2026). This is a wrapper
#' function for \code{\link[metafor]{rma.mv}} in the R package \bold{metafor}
#' (Viechtbauer 2010).
#'
#' @details
#' This R function is a wrapper function for \code{\link[metafor]{rma.mv}} in
#' the R package \bold{metafor} (Viechtbauer 2010).
#'
#' Note, results are not back-transformed in printouts of network meta-analyses
#' using summary measures with transformations, e.g., log risk ratios are
#' printed instead of the risk ratio if argument \code{sm = "RR"}.
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
#'   unrelated mean interaction effect (UMIE) model should be assumed.
#' @param assumption A character string indicating which assumption is done
#'   for the covariate; either "independent", "common", or "exchangeable" (can be abbreviated).
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance tau-squared. Either
#'   \code{"FE"}, \code{"REML"}, or \code{"ML"}.
#' @param level The level used to calculate confidence intervals for regression
#'   coefficients.
#' @param reference.group Reference treatment.
#' @param direction1 Directionality parameter for the first treatment group
#'   when using UMIE models.
#' @param direction2 Directionality parameter for the second treatment group
#'   when using UMIE models.
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
#' @param \dots Additional arguments passed to \code{\link[metafor]{rma.mv}}.
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
#' @author Nana-adjoa Kwarteng
#'   \email{nana-adjoa.kwarteng@@uniklinik-freiburg.de},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{netmeta}}
#'
#' @references
#' Kwarteng N, Evrenoglou T, Mueller J, Elsaesser M, Schramm E, Schwarzer G,
#' Nikolakopoulou A (2026):
#' Illustrating the assumptions of meta-regression in treatment networks.
#' Preprint available at \emph{Research Square},
#' \doi{10.21203/rs.3.rs-8235913/v1}
#' 
#' Viechtbauer W (2010):
#' Conducting meta-analyses in R with the metafor package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#'
#' @keywords models regression
#' 
#' @examples
#' #
#' # 1) Smoking cessation example
#' #
#' data(smokingcessation)
#' 
#' # Add variable with (fictitious) risk of bias values
#' # with 1 = "low risk" and 2 = "high risk"
#' #
#' smokingcessation$rob <- rep(1:2, 12)
#' 
#' # Only consider first five studies (to reduce runtime of example)
#' #
#' smoke5 <- head(smokingcessation, n = 5)
#' 
#' # Transform data from arm-based format to contrast-based format
#' #
#' pw1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smoke5, sm = "OR")
#' 
#' # Run network meta-analysis (NMA)
#' nma1 <- netmeta(pw1, common = FALSE, ref = "A")
#' 
#' # Network meta-regression with continuous covariate and assumption of
#' # independent slopes
#' nr1 <- netmetareg(nma1, rob)
#' nr1
#' 
#' \donttest{
#' # 2) Pain Prevention of Propofol Injection Example (Jalota2011)
#' #
#' # Create pairwise object
#' pw2 <- pairwise(treat = trt, event =  pain, n = n,
#'   studlab = id, data = Jalota2011, sm = "OR")
#' 
#' # Run network meta-analysis (NMA)
#' nma2 <- netmeta(pw2, common = FALSE, ref = "Hand vein")
#' 
#' # NMR with independent and consistent assumption (default)
#' nmr2_i <- netmetareg(nma2, covar = seTE, consistency = TRUE,
#'   assumption = "i")
#' # NMR with common and consistent assumption
#' nmr2_c <- netmetareg(nma2, covar = seTE, consistency = TRUE,
#'   assumption = "c")
#'   
#' # NMR with exchangeable and consistent assumption
#' nmr2_e <- netmetareg(nma2, covar = seTE, consistency = TRUE,
#'   assumption = "e")
#'   
#' # NMR with independent UMIE model (no consistency in interactions) using
#' # default treatment order
#' nmr2_i_umie <- netmetareg(nma2, covar = seTE, consistency = FALSE,
#'   assumption = "i")
#'
#' # 3) Physical therapy example (Hong 2015)
#' #
#' dat3 <- Hong2015
#' # Externally create a customizable treatment order indicator variable,
#' # which will later be modified
#' dat3$trtnum <- as.numeric(as.factor(dat3$trt))
#' dat3 <- do.call(rbind, lapply(split(dat3, dat3$id),
#'   function(x) {
#'     x$direction <- ifelse(x$trtnum == min(x$trtnum), 0, 1)
#'   x}))
#' 
#' pw3 <- pairwise(treat = trt, mean = mean_pain, n = n, sd = sd_pain,
#'   studlab=id, data = dat3, sm = "MD")
#' # Create difference variable for covariate
#' pw3$disability_diff <- pw3$mean_disability1 - pw3$mean_disability2
#' 
#' # Update directionality parameters so that they will exclude imputed values
#' sel.NA <- is.na(pw3$disability_diff)
#' pw3$direction1[sel.NA] <- 0
#' pw3$direction2[sel.NA] <- 0
#' 
#' # Mean imputations for missing covariate values
#' pw3$disability_diff[sel.NA] <- mean(pw3$disability_diff, na.rm = TRUE)
#' 
#' # NMA
#' nma3 <- netmeta(pw3, common = FALSE, ref = "No treatment")
#' 
#' # NMR with independent and consistent assumption (default)
#' nmr3_i <- netmetareg(nma3, covar = disability_diff, consistency = TRUE,
#'   assumption = "i")
#' # NMR with common and consistent assumption
#' nmr3_c <- netmetareg(nma3, covar = disability_diff, consistency = TRUE,
#'   assumption = "c")
#' # NMR with exchangeable and consistent assumption
#' nmr3_e <- netmetareg(nma3, covar = disability_diff, consistency = TRUE,
#'   assumption = "e")
#' # NMR with independent UMIE model (no consistency in interactions) using
#' # default treatment order variable
#' nmr3_i_umie <- netmetareg(nma3, covar = disability_diff, consistency = FALSE,
#'   assumption = "i")
#' # NMR with independent UMIE model (no consistency in interactions) using
#' # custom treatment order variable which excludes the mean imputations from
#' # interaction estimation
#' nmr3_i_umie_noimp <- netmetareg(nma3, covar = disability_diff,
#'   consistency = FALSE, assumption = "i",
#'    direction1 = "direction1", direction2 = "direction2")
#' # NMR with common UMIE model (no consistency in interactions) using default
#' # treatment order variable
#' nmr3_c_umie <- netmetareg(nma3, covar = disability_diff, consistency = FALSE,
#'   assumption = "c")
#' # NMR with common UMIE model (no consistency in interactions) using custom
#' # treatment order variable which excludes the mean imputations from
#' # interaction estimation
#' nmr3_c_umie_noimp <- netmetareg(nma3, covar = disability_diff,
#'   consistency = FALSE, assumption = "c",
#'   direction1 = "direction1", direction2 = "direction2")
#'   #' # NMR with exchangeable UMIE model (no consistency in interactions) using default
#' # treatment order variable
#' nmr3_e_umie <- netmetareg(nma3, covar = disability_diff, consistency = FALSE,
#'   assumption = "e")
#' # NMR with exchangeable UMIE model (no consistency in interactions) using custom
#' # treatment order variable which excludes the mean imputations from
#' # interaction estimation
#' nmr3_e_umie_noimp <- netmetareg(nma3, covar = disability_diff,
#'   consistency = FALSE, assumption = "e",
#'   direction1 = "direction1", direction2 = "direction2")
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
                               direction1 = NULL, direction2 = NULL,
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
  assumption <- setchar(assumption, c("independent", "common", "exchangeable"))
  #
  reference.group <- setref(reference.group, x$trts)
  method.tau <- setchar(method.tau, c("REML", "ML", "FE"))
  chklevel(level)
  #
  sm <- x$sm
  #
  args <- list(...)
  # Check whether first argument is a list. In this case only use
  # this list as input.
  if (length(args) > 0 && is.list(args[[1]]))
    args <- args[[1]]
  #
  max.ia <- replaceNULL(args[["max.ia"]], FALSE)
  chklogical(max.ia)
  #
  trts <- x$trts
  trts.abbr <- treats(x$trts, nchar.trts)
  #
  reference.group <- trts.abbr[trts == reference.group]
  
  
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
  covar.name <- deparse(substitute(covar))
  covar.name.internal <- make.names(covar.name) # R/matrix friendly covariate name
  #
  covar <- catch("covar", mc, data, sfsp)
  #
  if (!is.null(x$data)) {
    sel.data <- rep(TRUE, nrow(data))
    if (isCol(data, ".subset"))
      sel.data <- sel.data & data$.subset %in% TRUE
    if (isCol(data, ".excl"))
      sel.data <- sel.data & data$.excl %in% FALSE
    #
    if (length(covar) > 1L)
      covar <- covar[sel.data]
    data <- data[sel.data, , drop = FALSE]
    #
    x$data <- data
  }
  #
  trts <- x$trts
  trts.abbr <- treats(x$trts, nchar.trts)
  trts.abbr <- make.names(trts.abbr)# R/matrix friendly reference treat names
  #
  reference.group <- trts.abbr[trts == reference.group]
  reference.group.internal <- make.names(reference.group)# R/matrix friendly reference treat names
  
  dat <- data.frame(studlab = data$.studlab,
                    treat1 = data$.treat1,
                    treat2 = data$.treat2,
                    TE = data$.TE,
                    seTE = data$.seTE)
  #
  dat$treat1 <-
    as.character(factor(dat$treat1, levels = trts, labels = trts.abbr))
  dat$treat2 <-
    as.character(factor(dat$treat2, levels = trts, labels = trts.abbr))
  #
  # dat$treat1.internal <- make.names(dat$treat1)# R/matrix friendly treat names
  # dat$treat2.internal <- make.names(dat$treat2)# R/matrix friendly treat names
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
  dat[[covar.name.internal]] <- covar
  
  available.n <- available.events <- available.times <- available.sds <- FALSE
  if (!is.null(data$.n1) && !is.null(data$.n2)) {
    dat$n1 <- data$.n1
    dat$n2 <- data$.n2
    available.n <- TRUE
  }
  #
  if (!is.null(data$.mean1) && !is.null(data$.mean2)) {
    dat$mean1 <- data$.mean1
    dat$mean2 <- data$.mean2
  }
  #
  if (!is.null(data$.sd1) && !is.null(data$.sd2)) {
    dat$sd1 <- data$.sd1
    dat$sd2 <- data$.sd2
    available.sds <- TRUE
  }
  #
  if (!is.null(data$.event1) && !is.null(data$.event2)) {
    dat$event1 <- data$.event1
    dat$event2 <- data$.event2
    available.events <- TRUE
  }
  #
  if (!is.null(data$.time1) && !is.null(data$.time2)) {
    dat$time1 <- data$.time1
    dat$time2 <- data$.time2
    available.times <- TRUE
  }
  #
  if (!is.null(data$.incr1) && !is.null(data$.incr2)) {
    dat$incr1 <- data$.incr1
    dat$incr2 <- data$.incr2
  }
  
  # import or create directionality parameters
  if (!consistency) {
    # Handle directionality parameters for inconsistency model
    has_treat1_dir <- !is.null(substitute(direction1))
    has_treat2_dir <- !is.null(substitute(direction2))
    
    if (has_treat1_dir & has_treat2_dir) {
      # Use provided directionality parameters
      dat$Ival1 <-
        ifelse(data$treat1 == data$.treat1,
               data[[substitute(direction1)]],
               data[[substitute(direction2)]])
      #
      dat$Ival2 <-
        ifelse(data$treat2 == data$.treat2,
               data[[substitute(direction2)]],
               data[[substitute(direction1)]])
    }
    else if (!has_treat1_dir & !has_treat2_dir) {
      # Default to treatment order when no directionality specified
      dat <- do.call(rbind, lapply(split(dat, dat$studlab), mkIval))
      rownames(dat) <- NULL
      warning("No directionality parameters specified. ",
              "Defaulting to treatment order.")
    }
    # Validate directionality parameters
    check_Ival(dat)
  }
  #
  dat <- dat[order(dat$studlab, dat$treat1, dat$treat2), , drop = FALSE]
  #
  keep <- logical(0)
  wo <- logical(0)
  #
  if (consistency) {
    # Check for constant covariate in reference treatment
    if (length(unique(dat[dat$treat1 == reference.group.internal |
                          dat$treat2 == reference.group.internal,
                          covar.name.internal])) == 1)
      stop("Invalid reference treatment for interaction. ",
           "Insufficient variation in observed covariate values for ",
           "reference treatment.")
    #
    for (i in unique(dat$studlab)) {
      d.i <- dat[dat$studlab == i, , drop = FALSE]
      trts.i <- unique(sort(c(d.i$treat1, d.i$treat2)))
      if (reference.group.internal %in% trts.i)
        ref.i <- reference.group.internal
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
  }
  else if (!consistency) {
    # calcV() is sensitive to the choice of reference treatment, but if we
    # have study specific references we have to update the keep and wo
    # statements
    #
    dat$Ival_diff <- dat$Ival1 - dat$Ival2
    
    
    # Basic Grouping info, identify multi arm and directionality
    is_multi <- ave(seq_along(dat$studlab), dat$studlab, FUN = function(x) length(x) > 1)
    has_diff <- ave(seq_along(dat$Ival_diff), dat$studlab, FUN = function(x) any(dat$Ival_diff[x] != 0))
    
    # Equivalent !(d.i$treat1 != ref.i & d.i$treat2 != ref.i) from netmeta
    get_ref_name <- function(stud_idx) {
      # Subset treatments for just this study
      trts <- unique(sort(c(dat$treat1[stud_idx], dat$treat2[stud_idx])))
      if (reference.group.internal %in% trts) {
        return(reference.group.internal)
      } else {
        return(rev(trts)[1])
      }
    }
    
    # Create a vector of the reference treatment name for every row
    ref_names <- ave(dat$studlab, dat$studlab, FUN = function(x) {
      idx <- which(dat$studlab == x[1]) 
      get_ref_name(idx)
    })
    
    # Logic for the fallback: Keep if either treatment matches the named reference
    keep_named <- (dat$treat1 == ref_names | dat$treat2 == ref_names)
    
    # Condition A: 2-Arm study -> Always keep, even if Ival_diff is 0 so that the treatment effect remains
    # Condition B: Multi-Arm AND has a difference -> Keep if Ival_diff != 0, this was the original behavior
    # Condition C: Multi-Arm AND all are 0 -> Use the Named Reference Fallback o that so that the treatment effect remains
    keep.i <- !is_multi | 
      (is_multi & has_diff & dat$Ival_diff != 0) | 
      (is_multi & !has_diff & keep_named)
    
    # Equivalent to d.i$treat1 == ref.i
    # If Multi-Arm and all are 0, use the named reference for the 'wo' flag
    # Otherwise, use Ival1 == 0
    wo.i <- ifelse(is_multi & !has_diff, 
                   dat$treat1 == ref_names, # use the name if the reference isn't defined by Ival
                   dat$Ival1 == 0) # old logic
    
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
    t2.tmp <- dat$treat2
    e2.tmp <- dat$event2
    n2.tmp <- dat$n2
    mean2.tmp <- dat$mean2
    sd2.tmp <- dat$sd2
    time2.tmp <- dat$time2
    incr2.tmp <- dat$incr2
    Ival2.tmp <- dat$Ival2
    #
    dat$treat2[wo] <- dat$treat1[wo]
    dat$event2[wo] <- dat$event1[wo]
    dat$n2[wo] <- dat$n1[wo]
    dat$mean2[wo] <- dat$mean1[wo]
    dat$sd2[wo] <- dat$sd1[wo]
    dat$time2[wo] <- dat$time1[wo]
    dat$incr2[wo] <- dat$incr1[wo]
    dat$Ival2[wo] <- dat$Ival1[wo]
    #
    dat$treat1[wo] <- t2.tmp[wo]
    dat$event1[wo] <- e2.tmp[wo]
    dat$n1[wo] <- n2.tmp[wo]
    dat$mean1[wo] <- mean2.tmp[wo]
    dat$sd1[wo] <- sd2.tmp[wo]
    dat$time1[wo] <- time2.tmp[wo]
    dat$incr1[wo] <- incr2.tmp[wo]
    dat$Ival1[wo] <- Ival2.tmp[wo]
  }
  #
  ncols1 <- ncol(dat)
  dat <- contrmat(dat, grp1 = "treat1", grp2 = "treat2")
  dat <- dat[order(dat$studlab), ]
  trts.all <- names(dat)[(ncols1 + 1):ncol(dat)]
  #
  # Ensure that the reference treatment is correctly assigned
  trts <- trts.all[trts.all !=reference.group.internal]
  
  
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
                                   paste(paste0(trts, ":", covar.name.internal),
                                         collapse = " + "))))
    }
    else {
      #
      # The interaction is all non-reference treatments vs reference,
      # so we can define a variable indicating whether it is reference
      # or non-reference.
      #
      error <-
        try(dat$nonref <- as.numeric(dat[[reference.group.internal]] != 0),
            silent = TRUE)
      #
      # Necessary, if reference treatment contains a white space
      #
      if (inherits(error, "try-error"))
        dat$nonref <- as.numeric(dat[[reference.group.internal]] != 0)
      #
      # Then get the interaction using the colon(:). The output has a
      # colon which is in the same format as the independent. Helpful
      # for later extraction
      #
      formula.nmr_default <- # renamed for construction of manual matrix below
        as.formula(paste("~ 0 + ",
                         paste(trts, collapse = " + "),
                         if (!is.null(covar))
                           paste0( " + ",
                                   paste(paste0("nonref:", covar.name.internal),
                                         collapse = " + "))))
    } 
	if(assumption=="exchangeable"){
      formula.independent_matrix <- as.formula(paste("~ 0 ", if (!is.null(covar)) 
                                                paste0(" + ", paste(paste0(trts, ":", covar.name.internal), 
                                                                    collapse = " + "))))
    }
  }
  else {
    dat$Ival <- dat$Ival1 - dat$Ival2
      dat$.comp_ <- paste(pmin(dat$treat1, dat$treat2),
                          pmax(dat$treat1, dat$treat2),
                          sep = "_vs_")
      count_comp <- table(dat$.comp_)
      # for dev use
      if (!max.ia) {
        # Requires at least two observations
        dat$.comp_ <-
          ifelse(dat$.comp_ %in% names(count_comp[count_comp >= 2]),
                 dat$.comp_, "insufficient_data")
      }
    if(assumption=="independent"){
      formula.nmr_default <-
        as.formula(paste("~ 0 + ", paste(trts, collapse = " + "),
                         if (!is.null(covar))
                           paste0("+ .comp_:Ival:", covar.name.internal))
        )
    }
    else {
      formula.nmr_default <-
        as.formula(paste("~ 0 + ", paste(trts, collapse = " + "),
                         if (!is.null(covar))
                           paste0("+ Ival:", covar.name.internal))
        )
    }
    if(assumption=="exchangeable"){
      stop("exchangeable UMIE not yet implemented")
      formula.independent_matrix<- as.formula(paste("~ 0 ", if (!is.null(covar)) 
        paste0("+ .comp_:Ival:", covar.name.internal)))
    }
  }
  
  # irrespective of consistency
  if(assumption=="exchangeable") {
    Z_slope_matrix <- model.matrix(formula.independent_matrix, data = dat)
    Z_slope_matrix <- Z_slope_matrix[, grepl("insufficient_data", colnames(Z_slope_matrix)) == FALSE]
    if (is.factor(covar)) 
      Z_slope_matrix <- Z_slope_matrix[, !grepl(paste0(levels(covar)[1], "$"), colnames(Z_slope_matrix)), 
                                       drop = FALSE]
    else if (is.character(covar)) 
      Z_slope_matrix <- Z_slope_matrix[, !grepl(paste0(min(covar, na.rm = TRUE), "$"), 
                                                colnames(Z_slope_matrix)), drop = FALSE]
    else if (is.logical(covar)) 
      Z_slope_matrix <- Z_slope_matrix[, !grepl(paste0(as.logical(min(covar)), "$"), 
                                                colnames(Z_slope_matrix)), drop = FALSE]
  }
  
  #
  # Checks for non-numeric covariate
  #
  mm <- model.matrix(formula.nmr_default, data = dat) # default model matrix
  mm <- mm[, !grepl("insufficient_data", colnames(mm))]
  #
  # Drop interactions which contain the reference covariate level if covariate
  # is of mode factor, character, or logical
  #
  if (is.factor(covar))
    mm <-
    mm[, !grepl(paste0(levels(covar)[1], "$"), colnames(mm)), drop = FALSE]
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
  treat1 <- treat2 <- comparison <- NULL
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
  if (available.n & (available.events | available.times | (available.sds)))
    V <- bldiag(lapply(split(dat, dat$studlab), calcV, sm = sm))
  else
    V <- dat$seTE^2
  
  if(assumption %in% c("independent", "common")){
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
    res$tau2B<-0
  } else{
    suppressWarnings(res <- runNN(rma.mv.exch, list(yi = dat.TE,
                                                    V = V,
                                                    data = dat,
                                                    mods = formula.nmr,
                                                    random = as.call(~factor(comparison) |studlab),
                                                    Z_slope=Z_slope_matrix,
                                                    rho = 0.5, method = method.tau, level = 100 *level,
                                                    ...)))
  }
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
  res$full_results <- nmr_full_results_new(res)
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
                             print.tau=TRUE,
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
    paste0(round(100 * x$.netmeta$level, 1), "% CI")
  #
  dat$z <- formatN(dat$z, digits = digits.stat)
  dat$pval <- formatPT(dat$pval, digits = digits.pval)
  names(dat)[names(dat) == "pval"] <- "p-value"
  #
  if(x$.netmeta$assumption%in%c("common","exchangeable") & x$.netmeta$consistency==TRUE){
    rownames(dat)<-gsub("nonref:",paste0("Non-",x$.netmeta$reference.group,":"),rownames(dat))
  }
  if( x$.netmeta$consistency==FALSE){
    rownames(dat)<-gsub("Ival:",paste0("Directionality:"),rownames(dat))
    rownames(dat)<-gsub("\\.comp_","",rownames(dat))
  }
  
  prmatrix(dat, quote = FALSE, right = TRUE, ...)
  #
  # Get rid of warning 'Undefined global functions or variables'
  comparison <- NULL
  # add heterogeneity estimators
  if(print.tau){
    text.hetstat<-paste0("\nQuantifying heterogeneity:\n")
    text.hetstat<-paste0(text.hetstat,
                         formatPT(x=x$tau2,
                                  lab = TRUE, labval = "tau^2",
                                  digits = 4,
                                  lab.NA = "NA"),
                         # add tau
                         "; ", formatPT(x=sqrt(x$tau2),
                                        lab = TRUE, labval = "tau",
                                        digits = 4,
                                        lab.NA = "NA")
    )
    if(x$.netmeta$assumption%in%c("exchangeable")){
      text.hetstat<-paste0(text.hetstat,
                           "; ",
                           formatPT(x=x$tau2B,
                                    lab = TRUE, labval = "tau.B^2",
                                    digits = 4,
                                    lab.NA = "NA"),
                           # add tau
                           "; ", formatPT(x=sqrt(x$tau2B),
                                          lab = TRUE, labval = "tau.B",
                                          digits = 4,
                                          lab.NA = "NA")
      )
    }
    cat(text.hetstat,"\n")
  }
  if (!is.null(covar) & details.methods) {
    cat("\nDetails on network meta-regression methods:\n")
    #
    if (x$.netmeta$method.tau == "FE"){
      cat("- Common effects model\n")
    }else{
      cat("- Random effects model\n")
      cat("- Restricted maximum-likelihood estimator for tau^2\n")
    }
    if (consistency)
      cat("- Consistency model\n")
    else
      cat("- Unrelated mean interaction effect (UMIE) model\n")
    #
    if (assumption == "independent")
      cat("- Independent slopes\n")
    else if (assumption == "common")
      cat("- Common slope\n")
    else if (assumption == "exchangeable")
      cat("- Exchangeable slope\n")
    #
    cat(paste0("- Reference group: ", x$.netmeta$reference.group, "\n"))
  }
  else if (details.methods) {
    cat("\nDetails on network meta-analysis methods:\n")
    #
    if (x$.netmeta$method.tau == "FE"){
      cat("- Common effects model\n")
    }else{
      cat("- Random effects model\n")
      cat("- Restricted maximum-likelihood estimator for tau^2\n")
    }
  }
  
  showlegend<-sum(x$.netmeta$trts!=x$.netmeta$trts.abbr)>0# if any treatment names are modified (shortened or special characters)
  legendabbr(unique(x$.netmeta$trts),
             make.names(treats(x$.netmeta$trts, x$.netmeta$nchar.trts)),
             condition=showlegend)
  
  
  
  #
  invisible(NULL)
}
