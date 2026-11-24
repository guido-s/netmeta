fit_nma_check_warnings <- function(pw) {
  expect_no_warning(fit <- netmeta(pw))
  #
  fit
}

check_saturated <- function(x) {
  expect_equal(x$df.Q, 0, tolerance = 0)
  expect_true(is.na(x$tau))
  expect_true(is.na(x$tau2))
  expect_true(is.na(x$pval.Q))
  expect_equal(x$TE.random, x$TE.common)
  expect_equal(x$seTE.random, x$seTE.common)
}
