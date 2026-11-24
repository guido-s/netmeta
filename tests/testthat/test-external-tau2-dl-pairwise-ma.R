#
# Setup
#

dat <- data.frame(
  MD = c(-0.2, -0.8, -1.4),
  seMD = c(0.1, 0.2, 0.3),
  trt1 = "A",
  trt2 = "B",
  study = paste0("study", 1:3))

#
# Tests
#

test_that("pairwise networks retain the DerSimonian-Laird estimate", {
  nma <- netmeta(MD, seMD, trt1, trt2, study, data = dat, sm = "MD")
  #
  expect_identical(nma$method.tau, "DL")
  #
  ma <- metagen(MD, seMD, data = dat, sm = "MD", method.tau = "DL")
  #
  mf <- metafor::rma(yi = MD, sei = seMD, data = dat,
                      measure = "MD", method = "DL")
  #
  expect_equal(nma$df.Q, 2, tolerance = 0)
  expect_equal(nma$Q, ma$Q)
  expect_equal(nma$Q, mf$QE)
  expect_equal(nma$tau2, ma$tau2)
  expect_equal(nma$tau2, mf$tau2)
})
