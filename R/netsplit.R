netsplit <- function(x) {
  ##
  lowertri <- function(x)
    x[lower.tri(x)]
  ##
  comparison <- names(x$prop.direct.fixed)
  ##
  ## Fixed effect model
  ##
  TE.fixed <- lowertri(x$TE.fixed)
  seTE.fixed <- lowertri(x$seTE.fixed)
  TE.direct.fixed <- lowertri(x$TE.direct.fixed)
  seTE.direct.fixed <- lowertri(x$seTE.direct.fixed)
  TE.indirect.fixed <- lowertri(x$TE.indirect.fixed)
  seTE.indirect.fixed <- lowertri(x$seTE.indirect.fixed)
  ##
  prop.fixed <- x$prop.direct.fixed
  ##
  m.fixed <- metagen(TE.direct.fixed - TE.indirect.fixed,
                     sqrt(seTE.direct.fixed^2 + seTE.indirect.fixed^2))
  ##
  ## Random effects model
  ##
  TE.direct.random <- lowertri(x$TE.direct.random)
  seTE.direct.random <- lowertri(x$seTE.direct.random)
  TE.indirect.random <- lowertri(x$TE.indirect.random)
  seTE.indirect.random <- lowertri(x$seTE.indirect.random)
  TE.random <- lowertri(x$TE.random)
  seTE.random <- lowertri(x$seTE.random)
  ##
  prop.random <- x$prop.direct.random
  ##
  m.random <- metagen(TE.direct.random - TE.indirect.random,
                      sqrt(seTE.direct.random^2 + seTE.indirect.random^2))
  ##
  res <- list(comparison = comparison,
              ##
              prop.fixed = prop.fixed,
              TE.fixed = TE.fixed,
              seTE.fixed = seTE.fixed,
              TE.direct.fixed = TE.direct.fixed,
              seTE.direct.fixed = seTE.direct.fixed,
              TE.indirect.fixed = TE.indirect.fixed,
              seTE.indirect.fixed = seTE.indirect.fixed,
              ##
              TE.diff.fixed = m.fixed$TE,
              seTE.diff.fixed = m.fixed$seTE,
              lower.diff.fixed = m.fixed$lower,
              upper.diff.fixed = m.fixed$upper,
              zval.diff.fixed = m.fixed$zval,
              pval.diff.fixed = m.fixed$pval,
              ##
              prop.random = prop.random,
              TE.random = TE.random,
              seTE.random = seTE.random,
              TE.direct.random = TE.direct.random,
              seTE.direct.random = seTE.direct.random,
              TE.indirect.random = TE.indirect.random,
              seTE.indirect.random = seTE.indirect.random,
              ##
              TE.diff.random = m.random$TE,
              seTE.diff.random = m.random$seTE,
              lower.diff.random = m.random$lower,
              upper.diff.random = m.random$upper,
              zval.diff.random = m.random$zval,
              pval.diff.random = m.random$pval,
              ##
              comb.fixed = x$comb.fixed,
              comb.random = x$comb.random,
              sm = x$sm,
              level.comb = x$level.comb
              )
  ##
  class(res) <- "netsplit"
  ##
  res
}
