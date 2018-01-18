netsplit <- function(x, upper = TRUE,
                     reference.group = x$reference.group,
                     baseline.reference = x$baseline.reference,
                     sep.trts = x$sep.trts, quote.trts = "",
                     tol.direct = 0.0005) {
  
  
  meta:::chkclass(x, "netmeta")
  meta:::chklogical(upper)
  meta:::chklogical(baseline.reference)
  
  
  seq.comps <- rownames(x$Cov.fixed)
  ##
  trts <- matrix(unlist(compsplit(seq.comps, x$sep.trts)),
                 ncol = 2, byrow = TRUE)
  trts <- as.data.frame(trts, stringsAsFactors = FALSE)
  names(trts) <- c("treat1", "treat2")
  ##
  if (!upper) {
    ##
    ## Comparison names are column:row (and must be switched)
    ##
    t1 <- trts$treat1
    trts$treat1 <- trts$treat2
    trts$treat2 <- t1
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
      wo[grep(reference.group, trts$treat1)] <- TRUE
    else
      wo[grep(reference.group, trts$treat2)] <- TRUE
  }
  ##
  if (any(wo)) {
    t1.wo <- trts$treat1[wo]
    trts$treat1[wo] <- trts$treat2[wo]
    trts$treat2[wo] <- t1.wo
  }
  ##
  comparison <- as.character(interaction(paste(quote.trts, trts$treat1,
                                               quote.trts, sep = ""),
                                         paste(quote.trts, trts$treat2,
                                               quote.trts, sep = ""),
                                         sep = sep.trts))
  
  
  ##
  ## Set direct evidence estimates to NA if only indirect evidence is
  ## available
  ##
  sel.zero.fixed <- abs(x$P.fixed) < tol.direct
  ##
  TE.direct.fixed <- x$TE.direct.fixed
  seTE.direct.fixed <- x$seTE.direct.fixed
  lower.direct.fixed <- x$lower.direct.fixed
  upper.direct.fixed <- x$upper.direct.fixed
  zval.direct.fixed <- x$zval.direct.fixed
  pval.direct.fixed <- x$pval.direct.fixed
  ##
  TE.direct.fixed[sel.zero.fixed] <- NA
  seTE.direct.fixed[sel.zero.fixed] <- NA
  lower.direct.fixed[sel.zero.fixed] <- NA
  upper.direct.fixed[sel.zero.fixed] <- NA
  zval.direct.fixed[sel.zero.fixed] <- NA
  pval.direct.fixed[sel.zero.fixed] <- NA
  ##
  sel.zero.random <- abs(x$P.random) < tol.direct
  ##
  TE.direct.random <- x$TE.direct.random
  seTE.direct.random <- x$seTE.direct.random
  lower.direct.random <- x$lower.direct.random
  upper.direct.random <- x$upper.direct.random
  zval.direct.random <- x$zval.direct.random
  pval.direct.random <- x$pval.direct.random
  ##
  TE.direct.random[sel.zero.random] <- NA
  seTE.direct.random[sel.zero.random] <- NA
  lower.direct.random[sel.zero.random] <- NA
  upper.direct.random[sel.zero.random] <- NA
  zval.direct.random[sel.zero.random] <- NA
  pval.direct.random[sel.zero.random] <- NA
  
  
  ##
  ## Indirect estimate is NA if only direct evidence is available
  ##
  sel.one.fixed <- abs(x$P.fixed - 1) < tol.direct
  ##
  TE.indirect.fixed <- x$TE.indirect.fixed
  seTE.indirect.fixed <- x$seTE.indirect.fixed
  lower.indirect.fixed <- x$lower.indirect.fixed
  upper.indirect.fixed <- x$upper.indirect.fixed
  zval.indirect.fixed <- x$zval.indirect.fixed
  pval.indirect.fixed <- x$pval.indirect.fixed
  ##
  TE.indirect.fixed[sel.one.fixed] <- NA
  seTE.indirect.fixed[sel.one.fixed] <- NA
  lower.indirect.fixed[sel.one.fixed] <- NA
  upper.indirect.fixed[sel.one.fixed] <- NA
  zval.indirect.fixed[sel.one.fixed] <- NA
  pval.indirect.fixed[sel.one.fixed] <- NA
  ##
  sel.one.random <- abs(x$P.random - 1) < tol.direct
  ##
  TE.indirect.random <- x$TE.indirect.random
  seTE.indirect.random <- x$seTE.indirect.random
  lower.indirect.random <- x$lower.indirect.random
  upper.indirect.random <- x$upper.indirect.random
  zval.indirect.random <- x$zval.indirect.random
  pval.indirect.random <- x$pval.indirect.random
  ##
  TE.indirect.random[sel.one.random] <- NA
  seTE.indirect.random[sel.one.random] <- NA
  lower.indirect.random[sel.one.random] <- NA
  upper.indirect.random[sel.one.random] <- NA
  zval.indirect.random[sel.one.random] <- NA
  pval.indirect.random[sel.one.random] <- NA
  
  
  ##
  ## Fixed effect model
  ##
  fixed.low <- data.frame(comparison,
                          TE = lowertri(x$TE.fixed),
                          seTE = lowertri(x$seTE.fixed),
                          lower = lowertri(x$lower.fixed),
                          upper = lowertri(x$upper.fixed),
                          z = lowertri(x$zval.fixed),
                          p = lowertri(x$pval.fixed),
                          stringsAsFactors = FALSE)
  ##
  direct.fixed.low <- data.frame(comparison,
                                 TE = lowertri(TE.direct.fixed),
                                 seTE = lowertri(seTE.direct.fixed),
                                 lower = lowertri(lower.direct.fixed),
                                 upper = lowertri(upper.direct.fixed),
                                 z = lowertri(zval.direct.fixed),
                                 p = lowertri(pval.direct.fixed),
                                 stringsAsFactors = FALSE)
  ##
  indirect.fixed.low <- data.frame(comparison,
                                   TE = lowertri(TE.indirect.fixed),
                                   seTE = lowertri(seTE.indirect.fixed),
                                   lower = lowertri(lower.indirect.fixed),
                                   upper = lowertri(upper.indirect.fixed),
                                   z = lowertri(zval.indirect.fixed),
                                   p = lowertri(pval.indirect.fixed),
                                   stringsAsFactors = FALSE)
  ##
  fixed.upp <- data.frame(comparison,
                          TE = uppertri(x$TE.fixed),
                          seTE = uppertri(x$seTE.fixed),
                          lower = uppertri(x$lower.fixed),
                          upper = uppertri(x$upper.fixed),
                          z = uppertri(x$zval.fixed),
                          p = uppertri(x$pval.fixed),
                          stringsAsFactors = FALSE)
  ##
  direct.fixed.upp <- data.frame(comparison,
                                 TE = uppertri(TE.direct.fixed),
                                 seTE = uppertri(seTE.direct.fixed),
                                 lower = uppertri(lower.direct.fixed),
                                 upper = uppertri(upper.direct.fixed),
                                 z = uppertri(zval.direct.fixed),
                                 p = uppertri(pval.direct.fixed),
                                 stringsAsFactors = FALSE)
  ##
  indirect.fixed.upp <- data.frame(comparison,
                                   TE = uppertri(TE.indirect.fixed),
                                   seTE = uppertri(seTE.indirect.fixed),
                                   lower = uppertri(lower.indirect.fixed),
                                   upper = uppertri(upper.indirect.fixed),
                                   z = uppertri(zval.indirect.fixed),
                                   p = uppertri(pval.indirect.fixed),
                                   stringsAsFactors = FALSE)
  ##
  if (!upper) {
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
                              p = m.fixed$pval,
                              stringsAsFactors = FALSE)
  
  
  ##
  ## Random effects model
  ##
  random.low <- data.frame(comparison,
                           TE = lowertri(x$TE.random),
                           seTE = lowertri(x$seTE.random),
                           lower = lowertri(x$lower.random),
                           upper = lowertri(x$upper.random),
                           z = lowertri(x$zval.random),
                           p = lowertri(x$pval.random),
                           stringsAsFactors = FALSE)
  ##
  direct.random.low <- data.frame(comparison,
                                  TE = lowertri(TE.direct.random),
                                  seTE = lowertri(seTE.direct.random),
                                  lower = lowertri(lower.direct.random),
                                  upper = lowertri(upper.direct.random),
                                  z = lowertri(zval.direct.random),
                                  p = lowertri(pval.direct.random),
                                  stringsAsFactors = FALSE)
  ##
  indirect.random.low <- data.frame(comparison,
                                    TE = lowertri(TE.indirect.random),
                                    seTE = lowertri(seTE.indirect.random),
                                    lower = lowertri(lower.indirect.random),
                                    upper = lowertri(upper.indirect.random),
                                    z = lowertri(zval.indirect.random),
                                    p = lowertri(pval.indirect.random),
                                    stringsAsFactors = FALSE)
  ##
  predict.low <- data.frame(comparison,
                            lower = lowertri(x$lower.predict),
                            upper = lowertri(x$upper.predict),
                            stringsAsFactors = FALSE)
  ##
  random.upp <- data.frame(comparison,
                           TE = uppertri(x$TE.random),
                           seTE = uppertri(x$seTE.random),
                           lower = uppertri(x$lower.random),
                           upper = uppertri(x$upper.random),
                           z = uppertri(x$zval.random),
                           p = uppertri(x$pval.random),
                           stringsAsFactors = FALSE)
  ##
  direct.random.upp <- data.frame(comparison,
                                  TE = uppertri(TE.direct.random),
                                  seTE = uppertri(seTE.direct.random),
                                  lower = uppertri(lower.direct.random),
                                  upper = uppertri(upper.direct.random),
                                  z = uppertri(zval.direct.random),
                                  p = uppertri(pval.direct.random),
                                  stringsAsFactors = FALSE)
  ##
  indirect.random.upp <- data.frame(comparison,
                                    TE = uppertri(TE.indirect.random),
                                    seTE = uppertri(seTE.indirect.random),
                                    lower = uppertri(lower.indirect.random),
                                    upper = uppertri(upper.indirect.random),
                                    z = uppertri(zval.indirect.random),
                                    p = uppertri(pval.indirect.random),
                                    stringsAsFactors = FALSE)
  ##
  predict.upp <- data.frame(comparison,
                            lower = uppertri(x$lower.predict),
                            upper = uppertri(x$upper.predict),
                            stringsAsFactors = FALSE)
  ##
  if (!upper) {
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
                               p = m.random$pval,
                               stringsAsFactors = FALSE)
  
  
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
              quote.trts = quote.trts,
              ##
              tol.direct = tol.direct,
              backtransf = x$backtransf,
              ##
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- "netsplit"
  
  res
}
