netsplit <- function(x) {
  ##
  lowertri <- function(x)
    x[lower.tri(x)]
  ##
  comparison <- names(x$prop.direct.fixed)
  ##
  ## Fixed effect model
  ##
  ##
  prop.fixed <- x$prop.direct.fixed
  ##
  fixed <- list(TE = lowertri(x$TE.fixed),
                seTE = lowertri(x$seTE.fixed),
                lower = lowertri(x$lower.fixed),
                upper = lowertri(x$upper.fixed),
                z = lowertri(x$zval.fixed),
                p = lowertri(x$pval.fixed))
  ##
  direct.fixed <- list(TE = lowertri(x$TE.direct.fixed),
                       seTE = lowertri(x$seTE.direct.fixed),
                       lower = lowertri(x$lower.direct.fixed),
                       upper = lowertri(x$upper.direct.fixed),
                       z = lowertri(x$zval.direct.fixed),
                       p = lowertri(x$pval.direct.fixed))
  ##
  indirect.fixed <- list(TE = lowertri(x$TE.indirect.fixed),
                         seTE = lowertri(x$seTE.indirect.fixed),
                         lower = lowertri(x$lower.indirect.fixed),
                         upper = lowertri(x$upper.indirect.fixed),
                         z = lowertri(x$zval.indirect.fixed),
                         p = lowertri(x$pval.indirect.fixed))
  ##
  m.fixed <- metagen(direct.fixed$TE - indirect.fixed$TE,
                     sqrt(direct.fixed$seTE^2 + indirect.fixed$seTE^2),
                     level = x$level.comb)
  ##
  compare.fixed <- list(TE = m.fixed$TE,
                        seTE = m.fixed$seTE,
                        lower = m.fixed$lower,
                        upper = m.fixed$upper,
                        z = m.fixed$zval,
                        p = m.fixed$pval)
  ##
  ## Random effects model
  ##
  prop.random <- x$prop.direct.random
  ##
  random <- list(TE = lowertri(x$TE.random),
                 seTE = lowertri(x$seTE.random),
                 lower = lowertri(x$lower.random),
                 upper = lowertri(x$upper.random),
                 z = lowertri(x$zval.random),
                 p = lowertri(x$pval.random))
  ##
  direct.random <- list(TE = lowertri(x$TE.direct.random),
                        seTE = lowertri(x$seTE.direct.random),
                        lower = lowertri(x$lower.direct.random),
                        upper = lowertri(x$upper.direct.random),
                        z = lowertri(x$zval.direct.random),
                        p = lowertri(x$pval.direct.random))
  ##
  indirect.random <- list(TE = lowertri(x$TE.indirect.random),
                          seTE = lowertri(x$seTE.indirect.random),
                          lower = lowertri(x$lower.indirect.random),
                          upper = lowertri(x$upper.indirect.random),
                          z = lowertri(x$zval.indirect.random),
                          p = lowertri(x$pval.indirect.random))
  ##
  m.random <- metagen(direct.random$TE - indirect.random$TE,
                      sqrt(direct.random$seTE^2 + indirect.random$seTE^2),
                      level = x$level.comb)
  ##
  compare.random <- list(TE = m.random$TE,
                         seTE = m.random$seTE,
                         lower = m.random$lower,
                         upper = m.random$upper,
                         z = m.random$zval,
                         p = m.random$pval)
  ##
  res <- list(comparison = comparison,
              ##
              prop.fixed = prop.fixed,
              fixed = fixed,
              direct.fixed = direct.fixed,
              indirect.fixed = indirect.fixed,
              compare.fixed = compare.fixed,
              ##
              prop.random = prop.random,
              random = random,
              direct.random = direct.random,
              indirect.random = indirect.random,
              compare.random = compare.random,
              ##
              comb.fixed = x$comb.fixed,
              comb.random = x$comb.random,
              sm = x$sm,
              level.comb = x$level.comb,
              version = packageDescription("netmeta")$Version
              )
  ##
  class(res) <- "netsplit"
  ##
  res
}
