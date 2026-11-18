#
# Setup
#

#
# 1) Comparison with GLS
#

# A three-arm study with a singular full covariance block, plus two
# independent two-arm studies. Direct GLS uses only two independent
# contrasts from the three-arm study.
#
S <- diag(c(0.04, 0.09, 0.012, 1e6))
S[1, 2] <- S[2, 1] <- 0.01
X <- rbind(c(1, 0), c(0, 1), c(1, 0), c(0, 1))
y <- c(0.2, 0.5, 0.25, 0.3)
C <- solve(crossprod(X, solve(S, X)))
beta <- C %*% crossprod(X, solve(S, y))
#
dat <- data.frame(
  TE = c(0.2, 0.5, 0.3, 0.25, 0.3),
  se = sqrt(c(0.04, 0.09, 0.11, 0.012, 1e6)),
  t1 = c("A", "A", "B", "A", "A"),
  t2 = c("B", "C", "C", "B", "C"),
  study = c("multi", "multi", "multi", "precise", "ignored"),
  correlated = c(TRUE, TRUE, TRUE, FALSE, FALSE))

#
# 2) Impact of precision loss on netimpact()
#

se.base <- c(0.8, 0.4, 0.04, 0.1)


#
# Tests
#

#
# 1) Comparison with GLS
#

test_that("block covariance inversion agrees with GLS regardless of row order", {
  rows <- seq_len(nrow(dat))
  #
  for (rows.i in list(rows, rev(rows))) {
    nma <- netmeta(TE, se, t1, t2, study, correlated = correlated,
                    data = dat[rows.i, ], tau.preset = 0)
    #
    expect_lt(abs(nma$TE.common["A", "B"] - beta[1]), 1e-8)
    expect_lt(abs(nma$seTE.common["A", "B"] - sqrt(C[1, 1])), 1e-8)
    expect_lt(abs(nma$seTE.common["A", "C"] - sqrt(C[2, 2])), 1e-8)
  }
})

#
# 2) Impact of precision loss on netimpact()
#

test_that("precision loss for independent studies is invariant to scale", {
  # For independent studies of one pair, precision loss has a closed form.
  for (scale in c(0.1, 1, 10)) {
    se <- se.base * scale
    nma <- netmeta(rep(0.2, 4), se, rep("A", 4), rep("B", 4), studlab = 1:4,
                   tau.preset = 0)
    #
    impact <- netimpact(nma)$impact.common[, 1]
    w <- 1 / se^2
    expected <- (w - 1 / (1000 * max(se))^2) / sum(w)
    #
    expect_lt(max(abs(impact - expected)), 1e-7)
  }
})
