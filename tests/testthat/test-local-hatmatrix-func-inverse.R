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
nma <- netmeta(TE, seTE, treat1, treat2, studlab = comps,
               func.inverse = MASS::ginv, random = FALSE)

#
# Tests
#

test_that("Ruecker hat matrices are extracted from the network analysis", {
  hm <- hatmatrix(nma, method = "Ruecker")
  #
  expect_identical(hm$common, nma$H.matrix.common)
  expect_identical(hm$random, nma$H.matrix.random)
  expect_identical(hm$func.inverse, nma$func.inverse)
})

test_that("Ruecker warns and ignores an explicitly supplied inverse", {
  expect_warning(
    hm <- hatmatrix(nma, method = "Ruecker", func.inverse = invmat),
    "Argument 'func.inverse' ignored", fixed = TRUE)
  #
  expect_identical(hm$common, nma$H.matrix.common)
  expect_identical(hm$random, nma$H.matrix.random)
  expect_identical(hm$func.inverse, nma$func.inverse)
})

test_that("Davies uses the inverse stored in the network analysis", {
  for (type in c("short", "long")) {
    hm <- hatmatrix(nma, method = "Davies", type = type)
    # Each observed comparison belongs to a different two-arm study.
    expect_equal(
      unname(hm$common[comps, comps]),
      unname(nma$H.matrix.common))
    expect_identical(hm$func.inverse, "MASS::ginv")
    expect_equal(nrow(hm$common), if (type == "short") 15 else 45)
    expect_equal(ncol(hm$common), 15)
  }
})

test_that("network analyses accept inverse functions and character names", {
  for (fun in list(MASS::ginv, "ginv", "MASS::ginv")) {
    fit <- netmeta(TE, seTE, treat1, treat2, studlab = comps,
                   func.inverse = fun, random = FALSE)
    #
    expect_identical(fit$func.inverse, "MASS::ginv")
    expect_equal(fit$H.matrix.common, nma$H.matrix.common)
    #
    for (type in c("short", "long")) {
      expect_equal(
        hatmatrix(fit, method = "Davies", type = type)$common,
        hatmatrix(nma, method = "Davies", type = type)$common)
    }
  }
})

test_that("Davies accepts explicit inverse functions and character names", {
  for (type in c("short", "long")) {
    reference <- hatmatrix(nma, method = "Davies", type = type)$common
    #
    for (fun in list(MASS::ginv, "ginv", "MASS::ginv")) {
      hm <- hatmatrix(nma, method = "Davies", type = type,
                       func.inverse = fun)
      #
      expect_equal(hm$common, reference)
      expect_identical(hm$func.inverse, "MASS::ginv")
    }
  }
})

test_that("Davies calls an explicitly supplied user-defined inverse", {
  calls <- 0L
  #
  inverse <- function(X) {
    calls <<- calls + 1L
    MASS::ginv(X)
  }
  #
  for (type in c("short", "long")) {
    calls <- 0L
    hm <- hatmatrix(nma, method = "Davies", type = type,
                     func.inverse = inverse)
    #
    expect_gt(calls, 0L)
    expect_identical(hm$func.inverse, inverse)
    expect_equal(
      hm$common,
      hatmatrix(nma, method = "Davies", type = type)$common)
  }
})
