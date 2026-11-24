#
# Setup
#

dat <- data.frame(
  effect = 1:6, se.effect = 1:6,
  trt1 = c("G", "B", "B", "D", "A", "F"),
  trt2 = c("B", "C", "E", "E", "H", "A"),
  study = 1:6)
#
nma <- netmeta(effect, se.effect, trt1, trt2, studlab = study,
               data = dat,
               subset =
                 (trt1 %in% c("A", "F", "H") & trt2 %in% c("A", "F", "H")))

#
# Tests
#

test_that("heterogeneity is non-estimable in the selected saturated network", {
  expect_true(is.na(nma$tau))
  expect_true(is.na(nma$tau2))
})

test_that("direct random-effects results fall back to common effects when tau is NA", {
  expect_equal(nma$TE.direct.random, nma$TE.direct.common)
  expect_equal(nma$seTE.direct.random, nma$seTE.direct.common)
  #
  comparisons <- cbind(c("A", "F", "A", "H"), c("F", "A", "H", "A"))
  expect_equal(unname(nma$TE.direct.random[comparisons]), c(-6, 6, 5, -5))
  expect_equal(unname(nma$seTE.direct.random[comparisons]), c(6, 6, 5, 5))
})

test_that("network random-effects results fall back to common effects when tau is NA", {
  expect_equal(nma$TE.random, nma$TE.common)
  expect_equal(nma$seTE.random, nma$seTE.common)
})
