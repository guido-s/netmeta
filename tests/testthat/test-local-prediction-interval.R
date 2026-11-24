#
# Setup
#

dat <- data.frame(
  TE = c(0.1, 0.4),
  seTE = c(0.1, 0.2),
  treat1 = c("A", "A"),
  treat2 = c("B", "B"),
  studlab = 1:2)

#
# Tests
#

test_that("prediction intervals are available with one residual degree of freedom", {
  for (method in c("classic", "t-dist")) {
    nma <- netmeta(TE, seTE, treat1, treat2, studlab, data = dat,
                    prediction = TRUE, method.random.ci = method)
    #
    expected.seTE <- sqrt(nma$seTE.random^2 + nma$tau^2)
    expected <- ci(nma$TE.random, expected.seTE,
                   level = nma$level.predict, df = 1)
    #
    diag(expected$lower) <- diag(expected$upper) <- 0
    #
    expect_equal(nma$df.Q, 1, tolerance = 0)
    expect_true(nma$prediction)
    expect_equal(nma$seTE.predict, expected.seTE)
    expect_equal(nma$lower.predict, expected$lower)
    expect_equal(nma$upper.predict, expected$upper)
    expect_true(all(is.finite(nma$lower.predict)))
    expect_true(all(is.finite(nma$upper.predict)))
    #
    output <-
      capture.output(print(nma, reference.group = "A", prediction = TRUE))
    #
    expect_true(any(grepl(paste0(100 * nma$level.predict, "% PI"),
                          output, fixed = TRUE)))
    #
    split <- netsplit(nma)
    index <- cbind("B", "A")
    #
    expect_equal(split$predict$lower, unname(nma$lower.predict[index]))
    expect_equal(split$predict$upper, unname(nma$upper.predict[index]))
  }
})

test_that("prediction intervals are unavailable with zero residual degrees of freedom", {
  for (method in c("classic", "t-dist")) {
    nma <- netmeta(TE, seTE, treat1, treat2, studlab, data = dat[1, ],
                    prediction = TRUE, method.random.ci = method)
    #
    expect_equal(nma$df.Q, 0, tolerance = 0)
    expect_false(nma$prediction)
    expect_true(is.na(nma$tau))
    expect_true(all(is.na(nma$seTE.predict)))
    expect_true(all(is.na(nma$lower.predict)))
    expect_true(all(is.na(nma$upper.predict)))
    #
    split <- netsplit(nma)
    #
    expect_false(split$prediction)
    expect_true(all(is.na(split$predict$lower)))
    expect_true(all(is.na(split$predict$upper)))
  }
})
