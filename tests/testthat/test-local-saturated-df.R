#
# Setup
#

#
# 1) Two three-arm studies sharing one treatment
#

dat1 <- data.frame(
  trt = c("A", "B", "C", "C", "D", "E"),
  ev = c(20, 40, 30, 30, 20, 40),
  n = 100,
  study = rep(1:2, each = 3))
#
# Relabeling treatments
#
oldlabs <- c("E", "D", "C", "B", "A")
newlabs <- c("B", "A", "D", "C", "E")
#
dat2 <- dat1
dat2$trt <- as.character(factor(dat2$trt, levels = oldlabs, labels = newlabs))

pw1 <- pairwise(treat = trt, event = ev, n = n,
                studlab = study, data = dat1, sm = "OR")
#
pw2 <- pairwise(treat = trt, event = ev, n = n,
                studlab = study, data = dat2, sm = "OR")

#
# 2) Single five-arm study
#

dat3 <- data.frame(
  trt = LETTERS[1:5],
  ev = c(20, 25, 30, 35, 40),
  n = 100,
  study = "Five-arm study")
#
pw3 <- pairwise(treat = trt, event = ev, n = n,
                studlab = study, data = dat3, sm = "OR")

#
# 3) Mixture of one two-arm, one three-arm and one four-arm study
#

dat5 <- data.frame(
  trt = c("A", "B", "B", "C", "D", "D", "E", "F", "G", "A", "G"),
  ev = c(20, 30, 30, 20, 40, 40, 25, 35, 45, 20, 45),
  n = 100,
  study = c(rep("two", 2), rep("three", 3), rep("four", 4), rep("loop", 2)))
#
dat4 <- subset(dat5, study != "loop")

pw4 <- pairwise(treat = trt, event = ev, n = n,
                studlab = study, data = dat4, sm = "OR")
#
pw5 <- pairwise(treat = trt, event = ev, n = n,
                studlab = study, data = dat5, sm = "OR")

#
# 4) Data from pairwise meta-analysis
#

dat6 <- data.frame(
  MD = c(-0.2, -0.8, -1.4),
  seMD = c(0.1, 0.2, 0.3),
  trt1 = "A",
  trt2 = "B",
  study = paste0("study", 1:3))

#
# Tests
#

#
# 1) Two three-arm studies sharing one treatment
#

test_that("two three-arm studies form a saturated network", {
  check_saturated(fit_nma_check_warnings(pw1))
  check_saturated(fit_nma_check_warnings(pw2))
})

test_that("relabeling treatments preserves random-effects results", {
  fit1 <- fit_nma_check_warnings(pw1)
  fit2 <- fit_nma_check_warnings(pw2)
  #
  expect_equal(
    unname(fit1$TE.random[oldlabs, oldlabs]),
    unname(fit2$TE.random[newlabs, newlabs])
  )
  expect_equal(
    unname(fit1$seTE.random[oldlabs, oldlabs]),
    unname(fit2$seTE.random[newlabs, newlabs])
  )
})

#
# 2) Single five-arm study
#

test_that("a single five-arm study forms a saturated network", {
  check_saturated(fit_nma_check_warnings(pw3))
})

#
# 3) Mixture of one two-arm, one three-arm and one four-arm study
#

test_that("mixed-arm studies form a saturated network without a loop", {
  check_saturated(fit_nma_check_warnings(pw4))
})

test_that("closing one loop adds one residual degree of freedom", {
  fit5 <- fit_nma_check_warnings(pw5)
  #
  expect_equal(fit5$df.Q, 1, tolerance = 0)
  expect_true(is.finite(fit5$tau))
})

#
# 4) DerSimonian-Laird estimate in network and pairwise meta-analysis
#

test_that("pairwise networks retain the DerSimonian-Laird estimate", {
  fit6 <- netmeta(MD, seMD, trt1, trt2, study, data = dat6, sm = "MD",
                  method.tau = "DL")
  #
  weights <- 1 / dat6$seMD^2
  pooled <- weighted.mean(dat6$MD, weights)
  #
  Q <- sum(weights * (dat6$MD - pooled)^2)
  C <- sum(weights) - sum(weights^2) / sum(weights)
  tau2 <- max(0, (Q - 2) / C)
  #
  expect_equal(fit6$df.Q, 2, tolerance = 0)
  expect_equal(fit6$Q, Q)
  expect_equal(fit6$tau2, tau2)
})
