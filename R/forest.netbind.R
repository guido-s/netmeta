forest.netbind <- function(x,
                           pooled = ifelse(x$comb.random, "random", "fixed"),
                           ##
                           equal.size = FALSE,
                           ##
                           leftcols = "studlab",
                           leftlabs = "Treatment",
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
  meta:::chkclass(x, "netbind")
  ##
  chkchar <- meta:::chkchar
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  formatPT <- meta:::formatPT
  setchar <- meta:::setchar
  ##
  pooled <- setchar(pooled, c("fixed", "random"))
  ##
  chklogical(equal.size)
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.prop, min = 0, single = TRUE)
  chklogical(backtransf)
  ##
  chkchar(lab.NA)
  
  
  ##
  ##
  ## (2) Extract results for fixed effect and random effects model
  ##
  ##
  sel <- x$fixed$treat != x$reference.group
  ##
  if (pooled == "fixed") {
    m <- metagen(x$fixed$TE, x$fixed$seTE, studlab = x$fixed$name,
                 sm = x$sm, comb.fixed = FALSE, comb.random = FALSE,
                 byvar = x$fixed$treat, print.byvar = FALSE,
                 subset = x$fixed$treat != x$reference.group)
    ##
    m$TE <- x$fixed$TE[sel]
    m$seTE <- x$fixed$seTE[sel]
    m$lower <- x$fixed$lower[sel]
    m$upper <- x$fixed$upper[sel]
    m$zval <- x$fixed$zval[sel]
    m$pval <- x$fixed$pval[sel]
    ##
    m$col.study <- x$fixed$col.study[sel]
    m$col.square <- x$fixed$col.square[sel]
    m$col.square.lines <- x$fixed$col.square.lines[sel]
    m$col.inside <- x$fixed$col.inside[sel]
    ##
    text.pooled <- "Fixed Effects Model"
  }
  else {
    m <- metagen(x$random$TE, x$random$seTE, studlab = x$random$name,
                 sm = x$sm, comb.fixed = FALSE, comb.random = FALSE,
                 byvar = x$random$treat, print.byvar = FALSE,
                 subset = x$random$treat != x$reference.group)
    ##
    m$TE <- x$random$TE[sel]
    m$seTE <- x$random$seTE[sel]
    m$lower <- x$random$lower[sel]
    m$upper <- x$random$upper[sel]
    m$zval <- x$random$zval[sel]
    m$pval <- x$random$pval[sel]
    ##
    m$col.study <- x$random$col.study[sel]
    m$col.square <- x$random$col.square[sel]
    m$col.square.lines <- x$random$col.square.lines[sel]
    m$col.inside <- x$random$col.inside[sel]
    ##
    text.pooled <- "Random Effects Model"
  }
  ##
  if (missing(smlab))
    if (x$baseline.reference)
      smlab <- paste0("Comparison: other vs '",
                      x$reference.group, "'\n(",
                      text.pooled,
                      ")")
    else
      smlab <- paste0("Comparison: '",
                      x$reference.group, "' vs other\n(",
                      text.pooled,
                      ")")
  
  
  ##
  ##
  ## (3) Forest plot
  ##
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
         ##
         col.study = m$col.study,
         col.square = m$col.square,
         col.square.lines = m$col.square.lines,
         col.inside = m$col.inside,
         col.inside.fixed = "black",
         col.inside.random = "black",
         ##
         weight.study = if (equal.size) "same" else pooled,
         calcwidth.subgroup = TRUE,
         ...)
  
  
  invisible(NULL)
}
