netsplit <- function(x, lower = FALSE,
                     reference.group = x$reference.group,
                     baseline.reference = x$baseline.reference,
                     sep.trts = x$sep.trts, quote = "") {
  
  
  meta:::chkclass(x, "netmeta")
  meta:::chklogical(lower)
  meta:::chklogical(baseline.reference)
  
  
  seq.comps <- rownames(x$Cov.fixed)
  x.sep.trts <- x$sep.trts
  if (x.sep.trts == ".")
    x.sep.trts <- "\\."
  ##
  treats <- matrix(unlist(strsplit(seq.comps, x.sep.trts)),
                   ncol = 2, byrow = TRUE)
  treats <- as.data.frame(treats, stringsAsFactors = FALSE)
  names(treats) <- c("treat1", "treat2")
  ##
  if (lower) {
    ##
    ## Comparison names are column:row (and must be switched)
    ##
    t1 <- treats$treat1
    treats$treat1 <- treats$treat2
    treats$treat2 <- t1
  }
  
  
  ##
  ## Number of studies with direct comparisons
  ##
  k <- lowertri(x$A.matrix)
  
  
  ##
  ## Change treatment order if
  ## - reference group is specified, i.e., unequal to ""
  ## - reference group is first treatment (argument 'baseline.reference' is TRUE)
  ## - reference group is second treatment (argument 'baseline.reference' is FALSE)
  ##
  wo <- rep_len(FALSE, length(seq.comps))
  ##
  if (reference.group != "") {
    reference.group <- setref(reference.group, colnames(x$TE.fixed))
    ##
    if (baseline.reference)
      wo[grep(reference.group, treats$treat1)] <- TRUE
    else
      wo[grep(reference.group, treats$treat2)] <- TRUE
  }
  ##
  if (any(wo)) {
    t1.wo <- treats$treat1[wo]
    treats$treat1[wo] <- treats$treat2[wo]
    treats$treat2[wo] <- t1.wo
  }
  ##
  comparison <- as.character(interaction(paste(quote, treats$treat1,
                                               quote, sep = ""),
                                         paste(quote, treats$treat2,
                                               quote, sep = ""),
                                         sep = sep.trts))
  
  
  ##
  ## Fixed effect model
  ##
  fixed.low <- data.frame(comparison,
                          TE = lowertri(x$TE.fixed),
                          seTE = lowertri(x$seTE.fixed),
                          lower = lowertri(x$lower.fixed),
                          upper = lowertri(x$upper.fixed),
                          z = lowertri(x$zval.fixed),
                          p = lowertri(x$pval.fixed))
  ##
  direct.fixed.low <- data.frame(comparison,
                                 TE = lowertri(x$TE.direct.fixed),
                                 seTE = lowertri(x$seTE.direct.fixed),
                                 lower = lowertri(x$lower.direct.fixed),
                                 upper = lowertri(x$upper.direct.fixed),
                                 z = lowertri(x$zval.direct.fixed),
                                 p = lowertri(x$pval.direct.fixed))
  ##
  indirect.fixed.low <- data.frame(comparison,
                                   TE = lowertri(x$TE.indirect.fixed),
                                   seTE = lowertri(x$seTE.indirect.fixed),
                                   lower = lowertri(x$lower.indirect.fixed),
                                   upper = lowertri(x$upper.indirect.fixed),
                                   z = lowertri(x$zval.indirect.fixed),
                                   p = lowertri(x$pval.indirect.fixed))
  ##
  fixed.upp <- data.frame(comparison,
                          TE = uppertri(x$TE.fixed),
                          seTE = uppertri(x$seTE.fixed),
                          lower = uppertri(x$lower.fixed),
                          upper = uppertri(x$upper.fixed),
                          z = uppertri(x$zval.fixed),
                          p = uppertri(x$pval.fixed))
  ##
  direct.fixed.upp <- data.frame(comparison,
                                 TE = uppertri(x$TE.direct.fixed),
                                 seTE = uppertri(x$seTE.direct.fixed),
                                 lower = uppertri(x$lower.direct.fixed),
                                 upper = uppertri(x$upper.direct.fixed),
                                 z = uppertri(x$zval.direct.fixed),
                                 p = uppertri(x$pval.direct.fixed))
  ##
  indirect.fixed.upp <- data.frame(comparison,
                                   TE = uppertri(x$TE.indirect.fixed),
                                   seTE = uppertri(x$seTE.indirect.fixed),
                                   lower = uppertri(x$lower.indirect.fixed),
                                   upper = uppertri(x$upper.indirect.fixed),
                                   z = uppertri(x$zval.indirect.fixed),
                                   p = uppertri(x$pval.indirect.fixed))
  ##
  if (lower) {
    fixed <- fixed.low
    direct.fixed <- direct.fixed.low
    indirect.fixed <- indirect.fixed.low
    ##
    if (any(wo)) {
      fixed[wo, ] <- fixed.upp[wo, ]
      direct.fixed[wo, ] <- direct.fixed.upp[wo, ]
      indirect.fixed[wo, ] <- indirect.fixed.upp[wo, ]
    }
  }
  else {
    fixed <- fixed.upp
    direct.fixed <- direct.fixed.upp
    indirect.fixed <- indirect.fixed.upp
    ##
    if (any(wo)) {
      fixed[wo, ] <- fixed.low[wo, ]
      direct.fixed[wo, ] <- direct.fixed.low[wo, ]
      indirect.fixed[wo, ] <- indirect.fixed.low[wo, ]
    }
  }
  ##
  m.fixed <- metagen(direct.fixed$TE - indirect.fixed$TE,
                     sqrt(direct.fixed$seTE^2 + indirect.fixed$seTE^2),
                     level = x$level.comb)
  ##
  compare.fixed <- data.frame(comparison,
                              TE = m.fixed$TE,
                              seTE = m.fixed$seTE,
                              lower = m.fixed$lower,
                              upper = m.fixed$upper,
                              z = m.fixed$zval,
                              p = m.fixed$pval)
  
  
  ##
  ## Random effects model
  ##
  random.low <- data.frame(comparison,
                           TE = lowertri(x$TE.random),
                           seTE = lowertri(x$seTE.random),
                           lower = lowertri(x$lower.random),
                           upper = lowertri(x$upper.random),
                           z = lowertri(x$zval.random),
                           p = lowertri(x$pval.random))
  ##
  direct.random.low <- data.frame(comparison,
                                  TE = lowertri(x$TE.direct.random),
                                  seTE = lowertri(x$seTE.direct.random),
                                  lower = lowertri(x$lower.direct.random),
                                  upper = lowertri(x$upper.direct.random),
                                  z = lowertri(x$zval.direct.random),
                                  p = lowertri(x$pval.direct.random))
  ##
  indirect.random.low <- data.frame(comparison,
                                    TE = lowertri(x$TE.indirect.random),
                                    seTE = lowertri(x$seTE.indirect.random),
                                    lower = lowertri(x$lower.indirect.random),
                                    upper = lowertri(x$upper.indirect.random),
                                    z = lowertri(x$zval.indirect.random),
                                    p = lowertri(x$pval.indirect.random))
  ##
  predict.low <- data.frame(comparison,
                            lower = lowertri(x$lower.predict),
                            upper = lowertri(x$upper.predict))
  ##
  random.upp <- data.frame(comparison,
                           TE = uppertri(x$TE.random),
                           seTE = uppertri(x$seTE.random),
                           lower = uppertri(x$lower.random),
                           upper = uppertri(x$upper.random),
                           z = uppertri(x$zval.random),
                           p = uppertri(x$pval.random))
  ##
  direct.random.upp <- data.frame(comparison,
                                  TE = uppertri(x$TE.direct.random),
                                  seTE = uppertri(x$seTE.direct.random),
                                  lower = uppertri(x$lower.direct.random),
                                  upper = uppertri(x$upper.direct.random),
                                  z = uppertri(x$zval.direct.random),
                                  p = uppertri(x$pval.direct.random))
  ##
  indirect.random.upp <- data.frame(comparison,
                                    TE = uppertri(x$TE.indirect.random),
                                    seTE = uppertri(x$seTE.indirect.random),
                                    lower = uppertri(x$lower.indirect.random),
                                    upper = uppertri(x$upper.indirect.random),
                                    z = uppertri(x$zval.indirect.random),
                                    p = uppertri(x$pval.indirect.random))
  ##
  predict.upp <- data.frame(comparison,
                            lower = uppertri(x$lower.predict),
                            upper = uppertri(x$upper.predict))
  ##
  if (lower) {
    random <- random.low
    direct.random <- direct.random.low
    indirect.random <- indirect.random.low
    predict <- predict.low
    ##
    if (any(wo)) {
      random[wo, ] <- random.upp[wo, ]
      direct.random[wo, ] <- direct.random.upp[wo, ]
      indirect.random[wo, ] <- indirect.random.upp[wo, ]
      predict[wo, ] <- predict.upp[wo, ]
    }
  }
  else {
    random <- random.upp
    direct.random <- direct.random.upp
    indirect.random <- indirect.random.upp
    predict <- predict.upp
    ##
    if (any(wo)) {
      random[wo, ] <- random.low[wo, ]
      direct.random[wo, ] <- direct.random.low[wo, ]
      indirect.random[wo, ] <- indirect.random.low[wo, ]
      predict[wo, ] <- predict.low[wo, ]
    }
  }
  ##
  m.random <- metagen(direct.random$TE - indirect.random$TE,
                      sqrt(direct.random$seTE^2 + indirect.random$seTE^2),
                      level = x$level.comb)
  ##
  compare.random <- data.frame(comparison,
                               TE = m.random$TE,
                               seTE = m.random$seTE,
                               lower = m.random$lower,
                               upper = m.random$upper,
                               z = m.random$zval,
                               p = m.random$pval)
  
  
  res <- list(comparison = comparison,
              ##
              k = k,
              ##
              prop.fixed = x$prop.direct.fixed[seq.comps],
              fixed = fixed,
              direct.fixed = direct.fixed,
              indirect.fixed = indirect.fixed,
              compare.fixed = compare.fixed,
              ##
              prop.random = x$prop.direct.random[seq.comps],
              random = random,
              direct.random = direct.random,
              indirect.random = indirect.random,
              compare.random = compare.random,
              predict = predict,
              ##
              sm = x$sm,
              level.comb = x$level.comb,
              comb.fixed = x$comb.fixed,
              comb.random = x$comb.random,
              ##
              prediction = x$prediction,
              level.predict = x$level.predict,
              tau = x$tau,
              ##
              reference.group = reference.group,
              baseline.reference = baseline.reference,
              sep.trts = sep.trts,
              quote = quote,
              ##
              backtransf = x$backtransf,
              ##
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- "netsplit"
  
  res
}
