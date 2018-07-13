forest.netsplit <- function(x,
                            pooled = ifelse(x$comb.random, "random", "fixed"),
                            show = "both",
                            ##
                            subgroup = "comparison",
                            ##
                            overall = TRUE,
                            direct = TRUE,
                            indirect = TRUE,
                            prediction = x$prediction,
                            ##
                            text.overall = "Network estimate",
                            text.direct = "Direct estimate",
                            text.indirect = "Indirect estimate",
                            text.predict = "Prediction interval",
                            ##
                            type.overall,
                            type.direct,
                            type.indirect,
                            ##
                            col.square = "gray",
                            col.square.lines = col.square,
                            col.inside = "white",
                            col.diamond = "gray",
                            col.diamond.lines = "black",
                            col.predict = "red",
                            col.predict.lines = "black",
                            ##
                            equal.size = FALSE,
                            ##
                            leftcols,
                            leftlabs,
                            rightcols = c("effect", "ci"),
                            rightlabs = NULL,
                            ##
                            digits = gs("digits.forest"),
                            digits.prop = max(gs("digits.pval") - 2, 2),
                            ##
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
  ##
  chkchar <- meta:::chkchar
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  formatPT <- meta:::formatPT
  setchar <- meta:::setchar
  ##
  pooled <- setchar(pooled, c("fixed", "random"))
  ##
  subgroup <- setchar(subgroup, c("comparison", "estimate"))
  ##
  chklogical(overall)
  chklogical(direct)
  chklogical(indirect)
  chklogical(prediction)
  ##
  chkchar(text.overall)
  chkchar(text.direct)
  chkchar(text.indirect)
  chkchar(text.predict)
  ##
  missing.type.overall <- missing(type.overall)
  if (missing.type.overall)
    type.overall <- "diamond"
  else
    type.overall <- setchar(type.overall, c("diamond", "square"))
  ##
  if (missing(type.direct))
    type.direct <- "square"
  else
    type.direct <- setchar(type.direct, c("diamond", "square"))
  if (missing(type.indirect))
    type.indirect <- "square"
  else
    type.indirect <- setchar(type.indirect, c("diamond", "square"))
  ##
  chkchar(col.square)
  chkchar(col.square.lines)
  chkchar(col.inside)
  chkchar(col.diamond)
  chkchar(col.diamond.lines)
  chkchar(col.predict)
  chkchar(col.predict.lines)
  ##
  chklogical(equal.size)
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.prop, min = 0, single = TRUE)
  chklogical(backtransf)
  ##
  chkchar(lab.NA)
  ##
  if (pooled == "fixed") {
    if (!(missing(prediction)) & prediction)
      warning("Prediction intervals not shown for estimates from fixed effect model.")
    prediction <- FALSE
  }
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
  if (missing(leftlabs)) {
    leftlabs <- rep(NA, length(leftcols))
    leftlabs[leftcols == "studlab"] <- "Comparison"
    leftlabs[leftcols == "k"] <- "Number of\nStudies"
    leftlabs[leftcols == "prop"] <- "Direct\nEvidence"
  }
  ##
  n.subgroup <- direct + indirect + overall + prediction
  missing.smlab <- missing(smlab)
  ##
  if (n.subgroup == 1 & overall & missing.type.overall)
    type.overall <- "square"
  ##
  if (missing(text.predict))
    if (!(length(x$level.predict) == 0) &&
        x$level.comb != x$level.predict)
      text.predict <- paste(text.predict, " (",
                            round(x$level.predict * 100), "%-PI)", sep = "")
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  ## Check whether first argument is a list. In this case only use
  ## this list as input.
  if (length(args) > 0 && is.list(args[[1]]))
    args <- args[[1]]
  ##
  additional.arguments <- names(args)
  ##
  if (length(additional.arguments) > 0) {
    if (!is.na(charmatch("showa", additional.arguments)))
      if (!missing(show))
        warning("Deprecated argument 'showall' ignored as argument 'show' is also provided.")
      else {
        warning("Deprecated argument 'showall' has been replaced by argument 'show'.")
        show <- args[[charmatch("showa", additional.arguments)]]
        if (show)
          show <- "all"
        else
          show <- "both"
      }
  }
  ##
  show <- setchar(show, c("all", "both", "with.direct", "direct.only", "indirect.only"))


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
    dat.direct$prop <- formatPT(x$prop.fixed, digits = digits.prop)
    dat.indirect$prop <- NA
    dat.overall$prop <- NA
    ##
    if (missing.smlab)
      smlab <- "Fixed effect model"
  }
  else {
    dat.direct <- x$direct.random
    dat.indirect <- x$indirect.random
    dat.overall <- x$random
    ##
    dat.direct$prop <- formatPT(x$prop.random, digits = digits.prop)
    dat.indirect$prop <- NA
    dat.overall$prop <- NA
    ##
    if (missing.smlab)
      smlab <- "Random effects model"
  }
  ##
  if (missing.smlab & n.subgroup == 1)
    smlab <- paste(if (direct)
                     paste(text.direct, "\n", sep = ""),
                   if (indirect)
                     paste(text.indirect, "\n", sep = ""),
                   if (overall)
                     paste(text.overall, "\n", sep = ""),
                   "(",
                   tolower(smlab),
                   ")",
                   sep = "")
  ##
  dat.predict <- x$predict
  dat.predict$TE <- dat.predict$seTE <-
    dat.predict$z <- dat.predict$p <- dat.predict$prop <- NA
  dat.predict <- dat.predict[, c("comparison", "TE", "seTE",
                                 "lower", "upper", "z", "p", "prop")]
  ##
  dat.direct$comps <- dat.indirect$comps <-
    dat.overall$comps <- dat.predict$comps <- x$comparison
  ##
  dat.direct$k <- x$k
  dat.indirect$k <- dat.overall$k <- dat.predict$k <- NA
  ##
  dat.direct$evidence   <- text.direct
  dat.indirect$evidence <- text.indirect
  dat.overall$evidence  <- text.overall
  dat.predict$evidence  <- text.predict
  ##
  dat.direct$type.study <- type.direct
  dat.indirect$type.study <- type.indirect
  dat.overall$type.study <- type.overall
  dat.predict$type.study <- "predict"
  ##
  dat.direct$col.estimate <- if (type.direct == "square")
                               col.square
                             else
                               col.diamond
  dat.indirect$col.estimate <- if (type.indirect == "square")
                                 col.square
                               else
                                 col.diamond
  dat.overall$col.estimate <- if (type.overall == "square")
                                col.square
                              else
                                col.diamond
  ##
  dat.direct$col.lines <- if (type.direct == "square")
                            col.square.lines
                          else
                            col.diamond.lines
  dat.indirect$col.lines <- if (type.indirect == "square")
                              col.square.lines
                            else
                              col.diamond
  dat.overall$col.lines <- if (type.overall == "square")
                             col.square.lines
                           else
                             col.diamond.lines
  ##
  dat.predict$col.estimate <- col.predict
  dat.predict$col.lines <- col.predict.lines
  ##
  ## col.square.lines = col.square,
  ## col.inside = "white",
  ## col.diamond.lines = "black",
  ## col.predict.lines = "black",
  ##
  dat.predict$TE <- dat.overall$TE
  dat.predict$seTE <- sqrt(dat.overall$seTE^2 + x$tau^2)


  ##
  ##
  ## (3) Select treatment comparisons to show in forest plot
  ##
  ##
  if (show == "all")
    sel <- rep_len(TRUE, length(x$direct.fixed$TE))
  else if (show == "with.direct")
    sel <- (!is.na(x$direct.fixed$TE) & !is.na(x$direct.random$TE))
  else if (show == "both")
    sel <- (!is.na(x$direct.fixed$TE)  & !is.na(x$indirect.fixed$TE) &
            !is.na(x$direct.random$TE) & !is.na(x$indirect.random$TE))
  else if (show == "direct.only")
    sel <- (!is.na(x$direct.fixed$TE)  & is.na(x$indirect.fixed$TE) &
            !is.na(x$direct.random$TE) & is.na(x$indirect.random$TE))
  else if (show == "indirect.only")
    sel <- (is.na(x$direct.fixed$TE)  & !is.na(x$indirect.fixed$TE) &
            is.na(x$direct.random$TE) & !is.na(x$indirect.random$TE))
  ##
  dat.direct <- dat.direct[sel, ]
  dat.indirect <- dat.indirect[sel, ]
  dat.overall <- dat.overall[sel, ]
  dat.predict <- dat.predict[sel, ]


  ##
  ##
  ## (4) Forest plot
  ##
  ##
  if (subgroup == "comparison") {
    dat <- rbind(
      if (direct) dat.direct,
      if (indirect) dat.indirect,
      if (overall) dat.overall,
      if (prediction) dat.predict
    )
    ##
    if (n.subgroup > 1)
      m <- metagen(dat$TE, dat$seTE, studlab = dat$evidence, data = dat,
                   sm = x$sm, byvar = dat$comps, print.byvar = FALSE)
    else
      m <- metagen(dat$TE, dat$seTE, studlab = dat$comps, data = dat, sm = x$sm)
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
           col.square = dat$col.estimate,
           col.square.lines = dat$col.lines,
           weight.study = if (equal.size) "same" else "fixed",
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
    if (n.subgroup > 1)
      m <- metagen(dat$TE, dat$seTE, studlab = dat$comps, data = dat,
                   sm = x$sm, byvar = dat$evidence, print.byvar = FALSE)
    else
      m <- metagen(dat$TE, dat$seTE, studlab = dat$comps, data = dat, sm = x$sm)
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
           col.square = dat$col.estimate,
           col.square.lines = dat$col.lines,
           weight.study = if (equal.size) "same" else "fixed",
           ...)
  }


  invisible(NULL)
}
