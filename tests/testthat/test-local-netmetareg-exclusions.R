#
# Setup
#

dat <- data.frame(
  study = rep(paste0("S", 1:7), each = 2),
  treatment = c("A", "B", "A", "C", "B", "C", "A", "B",
                "A", "C", "B", "C", "A", "D"),
  event = c(10, 20, 15, 25, 18, 22, 14, 23, 20, 31, 16, 28, 0, 0),
  n = 100,
  year = rep(c(2000, 2001, 2002, 2010, 2011, 2012, 2020), each = 2))
#
pw <- pairwise(treat = treatment, event = event, n = n, studlab = study,
                data = dat, sm = "OR", allstudies = FALSE)

#
# Tests
#

test_that("meta-regression excludes comparisons omitted from the network", {
  expect_warning(
    nma1 <- netmeta(pw, reference.group = "A"),
    "Comparison with missing TE / seTE or zero seTE")
  #
  nma2 <- netmeta(subset(pw, studlab != "S7"), reference.group = "A")
  #
  nmr1 <- netmetareg(nma1, year, method.tau = "FE")
  nmr2 <- netmetareg(nma2, year, method.tau = "FE")
  #
  expect_equal(coef(nmr1), coef(nmr2))
  expect_equal(vcov(nmr1), vcov(nmr2))
  expect_setequal(nmr1$.netmeta$x$data$.studlab, paste0("S", 1:6))
  expect_equal(nrow(nma1$data), 7)
})

test_that("meta-regression aligns external covariates with retained comparisons", {
  expect_warning(
    nma1 <- netmeta(pw, reference.group = "A"),
    "Comparison with missing TE / seTE or zero seTE")
  #
  moderator <- pw$year
  #
  nmr1 <- netmetareg(nma1, moderator, method.tau = "FE")
  nmr2 <- netmetareg(nma1, year, method.tau = "FE")
  #
  expect_equal(unname(coef(nmr1)), unname(coef(nmr2)))
  expect_equal(unname(vcov(nmr1)), unname(vcov(nmr2)))
})

test_that("meta-regression respects subset and exclusion flags together", {
  expect_warning(
    nma1 <- netmeta(pw, reference.group = "A", subset = studlab != "S6"),
    "Comparison with missing TE / seTE or zero seTE")
  #
  nma2 <-
    netmeta(subset(pw, !studlab %in% c("S6", "S7")), reference.group = "A")
  #
  nmr2 <- netmetareg(nma2, year, method.tau = "FE")
  #
  moderator <- pw$year
  #
  nmr1 <- netmetareg(nma1, year, method.tau = "FE")
  external <- netmetareg(nma1, moderator, method.tau = "FE")
  #
  expect_equal(coef(nmr1), coef(nmr2))
  expect_equal(vcov(nmr1), vcov(nmr2))
  expect_equal(unname(coef(external)), unname(coef(nmr2)))
  expect_equal(unname(vcov(external)), unname(vcov(nmr2)))
  expect_setequal(nmr1$.netmeta$x$data$.studlab, paste0("S", 1:5))
  expect_false(anyNA(nmr1$.netmeta$x$data$.studlab))
})

test_that("meta-regression respects subset without missing comparisons", {
  nma1 <- netmeta(subset(pw, studlab != "S7"), reference.group = "A",
                  subset = studlab != "S6")
  nma2 <-
    netmeta(subset(pw, !studlab %in% c("S6", "S7")), reference.group = "A")
  #
  nmr1 <- netmetareg(nma1, year, method.tau = "FE")
  nmr2 <- netmetareg(nma2, year, method.tau = "FE")
  #
  expect_equal(coef(nmr1), coef(nmr2))
  expect_equal(vcov(nmr1), vcov(nmr2))
})
