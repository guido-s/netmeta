forest.netsplit <- function(x,
                            pooled = ifelse(x$comb.random, "random", "fixed"),
                            showall = FALSE,
                            overall = TRUE,
                            direct = TRUE,
                            indirect = TRUE,
                            prediction = x$prediction,
                            subgroup = "comparison",
                            leftcols,
                            leftlabs,
                            rightcols = c("effect", "ci"),
                            rightlabs = NULL,
                            digits = gs("digits.forest"),
                            backtransf = x$backtransf,
                            lab.NA = "",
                            smlab,
                            ...) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  meta:::chkclass(x, "netsplit")
  ## x <- upgradenetmeta(x)
  ##
  pooled <- meta:::setchar(pooled, c("fixed", "random"))
  ##
  meta:::chklogical(showall)
  meta:::chklogical(overall)
  meta:::chklogical(direct)
  meta:::chklogical(indirect)
  meta:::chklogical(prediction)
  ##
  subgroup <- meta:::setchar(subgroup, c("comparison", "estimate"))
  ##
  meta:::chknumeric(digits, min = 0, single = TRUE)
  meta:::chklogical(backtransf)
  ##
  meta:::chkchar(lab.NA)
  ##
  if (pooled == "fixed")
    prediction <- FALSE
  ##
  if (!any(c(overall, direct, indirect)))
    stop("At least, one of the following estimates must be included in forest plot:\n- network estimates (argument 'overall')\n- direct estimates (argument 'direct')\n- indirect estimates (argument 'indirect')")
  ##
  if (missing(leftcols))
    if (direct)
      leftcols <- c("studlab", "k", "prop")
    else
      leftcols <- "studlab"
  ##
  if (missing(leftlabs))
    if (direct)
      leftlabs <- c("Comparison", "Number of\nStudies", "Direct\nEvidence")
    else
      leftlabs <- "Comparison"
  
  
  ##
  ##
  ## (2) Extract results for fixed effect and random effects model
  ##
  ##
  if (pooled == "fixed") {
    dat.direct <- x$direct.fixed
    dat.indirect <- x$indirect.fixed
    dat.overall <- x$fixed
    ##
    dat.direct$prop <- x$prop.fixed
    dat.indirect$prop <- NA
    dat.overall$prop <- NA
    ##
    if (missing(smlab))
      smlab <- "Fixed effect model"
  }
  else {
    dat.direct <- x$direct.random
    dat.indirect <- x$indirect.random
    dat.overall <- x$random
    ##
    dat.direct$prop <- x$prop.random
    dat.indirect$prop <- NA
    dat.overall$prop <- NA
    ##
    if (missing(smlab))
      smlab <- "Random effects model"
  }
  ##
  dat.predict <- data.frame(TE = NA, seTE = NA,
                            lower = x$predict$lower,
                            upper = x$predict$upper,
                            z = NA, p = NA, prop = NA)
  ##
  dat.direct$comps <- dat.indirect$comps <-
    dat.overall$comps <- dat.predict$comps <- x$comparison
  ##
  dat.direct$k <- x$k
  dat.indirect$k <- dat.overall$k <- dat.predict$k <- NA
  ##
  dat.direct$evidence   <- " Direct estimate"
  dat.indirect$evidence <- " Indirect estimate"
  dat.overall$evidence  <- " Network estimate"
  dat.predict$evidence  <- " Prediction interval"
  ##
  dat.direct$type.study <- "square"
  dat.indirect$type.study <- "square"
  dat.overall$type.study <- "diamond"
  dat.predict$type.study <- "diamond" # prediction
  ##
  dat.predict$TE <- dat.overall$TE
  dat.predict$seTE <- sqrt(dat.overall$seTE^2 + x$tau^2)
  ##
  if (!showall) {
    dat.direct <- dat.direct[x$k > 0, ]
    dat.indirect <- dat.indirect[x$k > 0, ]
    dat.overall <- dat.overall[x$k > 0, ]
    dat.predict <- dat.predict[x$k > 0, ]
  }
  
  
  if (subgroup == "comparison") {
    dat <- rbind(
      if (direct) dat.direct,
      if (indirect) dat.indirect,
      if (overall) dat.overall,
      if (prediction) dat.predict
    )
    ##
    n.rows <- direct + indirect + overall + prediction
    ##
    if (n.rows > 1)
      m <- metagen(dat$TE, dat$seTE, studlab = dat$evidence, data = dat,
                   byvar = dat$comps, print.byvar = FALSE)
    else
      m <- metagen(dat$TE, dat$seTE, studlab = dat$comps, data = dat)
    ##
    forest(m,
           digits = digits,
           comb.fixed = FALSE, comb.random = FALSE,
           hetstat = FALSE,
           leftcols = leftcols,
           leftlabs = leftlabs,
           rightcols = rightcols,
           rightlabs = rightlabs,
           lab.NA = lab.NA,
           smlab = smlab,
           backtransf = backtransf,
           type.study = dat$type.study,
           ...)
  }
  else {
    dat <- rbind(
      if (direct) dat.direct,
      if (indirect) dat.indirect,
      if (overall) dat.overall,
      if (prediction) dat.predict
    )
    ##
    n.rows <- direct + indirect + overall + prediction
    ##
    if (n.rows > 1)
      m <- metagen(dat$TE, dat$seTE, studlab = dat$comps, data = dat,
                   byvar = dat$evidence, print.byvar = FALSE)
    else
      m <- metagen(dat$TE, dat$seTE, studlab = dat$comps, data = dat)
    ##
    forest(m,
           digits = digits,
           comb.fixed = FALSE, comb.random = FALSE,
           hetstat = FALSE,
           leftcols = leftcols,
           leftlabs = leftlabs,
           rightcols = rightcols,
           rightlabs = rightlabs,
           lab.NA = lab.NA,
           backtransf = backtransf,
           smlab = smlab,
           type.study = dat$type.study,
           ...)
  }
  
  
  invisible(NULL)
}
