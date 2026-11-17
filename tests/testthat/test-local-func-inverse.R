#
# Setup
#

#
# 1) netmeta()
#

data(Woods2010, package = "metabook")
#
pw1 <- pairwise(
  treat = treatment, event = r, n = N, studlab = author,
  data = Woods2010, sm = "OR")
#
nma1 <- netmeta(pw1)
nma1.ginv <- netmeta(pw1, func.inverse = ginv)

#
# 2) netmetabin()
#

data(Dogliotti2014, package = "metabook")
#
pw2 <- pairwise(
  treat = treatment, n = total, event = stroke, studlab = study,
  data = subset(Dogliotti2014, id %in% 1:3), sm = "OR")
#
nma2 <- netmeta(pw2)
nma2.ginv <- netmeta(pw2, func.inverse = ginv)
#
nmabin2 <- netmetabin(pw2, method = "Inverse")
nmabin2.ginv <- netmetabin(pw2, method = "Inverse", func.inverse = ginv)

#
# 3) netcomb() and discomb()
#

data(Linde2016, package = "metabook")
#
pw3 <- pairwise(
  treat = int, event = resp, n = n, studlab = paste(author, year),
  data = subset(Linde2016, id %in% c(16, 24, 49, 118)), sm = "OR")
#
nma3 <- netmeta(pw3)
nma3.ginv <- netmeta(pw3, func.inverse = ginv)
#
cnma3 <- netcomb(nma3, inactive = "placebo")
cnma3.ginv <- netcomb(nma3.ginv, inactive = "placebo")
cnma3.ginv2 <- netcomb(nma3, inactive = "placebo", func.inverse = ginv)
#
dnma3 <- discomb(pw3, inactive = "placebo")
dnma3.ginv <- discomb(pw3, inactive = "placebo", func.inverse = ginv)

fits <- list(
  netmeta = function(fun) netmeta(pw1, func.inverse = fun),
  netmetabin = function(fun)
    netmetabin(pw2, method = "Inverse", func.inverse = fun),
  netcomb = function(fun)
    netcomb(nma3, inactive = "placebo", func.inverse = fun),
  discomb = function(fun)
    discomb(pw3, inactive = "placebo", func.inverse = fun))
#
reference.default <- list(nma1, nmabin2, cnma3, dnma3)
reference.ginv <- list(nma1.ginv, nmabin2.ginv, cnma3.ginv2, dnma3.ginv)
#
results <- c("TE.common", "seTE.common", "TE.random", "seTE.random", "tau2")

#
# Tests
#

for (i in seq_along(fits)) {
  test_that(paste(names(fits)[i], "accepts inverse functions and character names"), {
    for (fun in list(MASS::ginv, "ginv", "MASS::ginv")) {
      fit <- fits[[i]](fun)
      #
      expect_identical(fit$func.inverse, "MASS::ginv")
      expect_equal(fit[results], reference.ginv[[i]][results])
    }
    #
    for (fun in list(invmat, "invmat", "netmeta::invmat")) {
      fit <- fits[[i]](fun)
      #
      expect_identical(fit$func.inverse, "netmeta::invmat")
      expect_equal(fit[results], reference.default[[i]][results])
    }
  })

  test_that(paste(names(fits)[i], "preserves and calls a user-defined inverse"), {
    calls <- 0L
    inverse <- function(X) {
      calls <<- calls + 1L
      MASS::ginv(X)
    }
    #
    fit <- fits[[i]](inverse)
    #
    expect_gt(calls, 0L)
    expect_identical(fit$func.inverse, inverse)
    expect_equal(fit[results], reference.ginv[[i]][results])
  })

  test_that(paste(names(fits)[i], "rejects invalid inverse functions"), {
    for (fun in list(NULL, 1, TRUE)) {
      expect_error(
        fits[[i]](fun),
        "Argument 'func.inverse' must specify a function.", fixed = TRUE)
    }
  })
}

test_that("netcomb inherits the inverse from the network analysis", {
  expect_identical(cnma3$func.inverse, nma3$func.inverse)
  expect_identical(cnma3.ginv$func.inverse, nma3.ginv$func.inverse)
  expect_equal(cnma3.ginv[results], cnma3.ginv2[results])
})

test_that("inverse-variance netmetabin results agree with netmeta", {
  expect_equal(nmabin2[results], nma2[results])
  expect_equal(nmabin2.ginv[results], nma2.ginv[results])
})

test_that("netmetabin forwards the inverse to every internal network analysis", {
  original.netmeta <- netmeta
  calls <- 0L
  inverse <- function(X) MASS::ginv(X)
  #
  local_mocked_bindings(
    netmeta = function(..., func.inverse = invmat) {
      calls <<- calls + 1L
      expect_identical(func.inverse, inverse)
      original.netmeta(..., func.inverse = func.inverse)
    },
    .package = "netmeta")
  #
  netmetabin(pw2, method = "Inverse", func.inverse = inverse)
  expect_gt(calls, 0L)
})

test_that("discomb forwards the inverse to its reference network analysis", {
  original.netmeta <- netmeta
  calls <- 0L
  inverse <- function(X) MASS::ginv(X)
  #
  local_mocked_bindings(
    netmeta = function(..., func.inverse = invmat) {
      calls <<- calls + 1L
      expect_identical(func.inverse, inverse)
      original.netmeta(..., func.inverse = func.inverse)
    },
    .package = "netmeta")
  #
  discomb(pw3, inactive = "placebo", func.inverse = inverse)
  expect_gt(calls, 0L)
})

test_that("Mantel-Haenszel netmetabin stores inverse functions and names", {
  reference <- netmetabin(pw2, method = "MH", func.inverse = MASS::ginv)
  #
  for (fun in list("ginv", "MASS::ginv")) {
    fit <- netmetabin(pw2, method = "MH", func.inverse = fun)
    #
    expect_identical(fit$func.inverse, "MASS::ginv")
    expect_equal(fit[results], reference[results])
  }
  #
  inverse <- function(X) MASS::ginv(X)
  fit <- netmetabin(pw2, method = "MH", func.inverse = inverse)
  #
  expect_identical(reference$func.inverse, "MASS::ginv")
  expect_identical(fit$func.inverse, inverse)
  expect_equal(fit[results], reference[results])
})
