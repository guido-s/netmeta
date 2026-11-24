#
# Setup
#

treat1 <- c("acar", "acar", "acar", "benf", "metf", "metf", "metf", "metf",
            "migl", "piog", "piog", "plac", "plac", "plac", "rosi")
#
treat2 <- c("metf", "plac", "sulf", "plac", "piog", "plac", "rosi", "sulf",
            "plac", "plac", "rosi", "rosi", "sita", "vild", "sulf")
#
comps <- paste(treat1, treat2, sep = ":")
#
TE <- rep(0, 15)
#
seTE <- c(5.159883e-01, 1.308290e-01, 1.662456e-01, 1.000000e+03,
          5.053986e-02, 1.194157e-01, 9.216967e-02, 1.786125e-01,
          1.000000e+03, 2.467612e-01, 6.232506e-02, 1.688789e-02,
          1.000000e+03, 1.000000e+03, 9.397779e-02)
#
nma <- netmeta(TE, seTE, treat1, treat2, studlab = comps, random = FALSE)

#
# Tests
#

test_that("Ruecker hat matrices are extracted from the network analysis", {
  hm <- hatmatrix(nma, method = "Ruecker")
  #
  expect_identical(hm$common, nma$H.matrix.common)
  expect_identical(hm$random, nma$H.matrix.random)
})

test_that("deprecated inverse input is ignored with a warning", {
  for (method in c("Ruecker", "Davies")) {
    reference <- hatmatrix(nma, method = method)
    #
    expect_warning(
      hm <- hatmatrix(nma, method = method, func.inverse = invmat),
      "Argument 'func.inverse' is deprecated and ignored.", fixed = TRUE)
    #
    expect_identical(hm$common, reference$common)
    expect_identical(hm$random, reference$random)
  }
})

test_that("Davies calculates the hat matrix stored in the network analysis", {
  for (type in c("short", "long")) {
    hm <- hatmatrix(nma, method = "Davies", type = type)
    # Each observed comparison belongs to a different two-arm study.
    expect_equal(
      unname(hm$common[comps, comps]),
      unname(nma$H.matrix.common))
    expect_equal(nrow(hm$common), if (type == "short") 15 else 45)
    expect_equal(ncol(hm$common), 15)
  }
})
