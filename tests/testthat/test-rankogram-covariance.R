#
# Setup
#

#
# 1) Three independent two-arm studies on a tree
#

dat1 <- data.frame(
  MD =  c(0.3, -0.2, 0.15), seMD = c(0.1, 0.5, 0.12),
  trt1 = c("A", "B", "C"), trt2 = c("B", "C", "D"),
  study = 1:3)

#
# 2) Analytical SUCRA from pairwise normal probabilities
#

pw2 <- pairwise(
  treat = c("A", "B", "B", "C", "C", "D"),
  event = c(62, 34, 60, 66, 237, 208),
  n = c(100, 100, 100, 100, 500, 500),
  studlab = rep(1:3, each = 2),
  sm = "OR")
#
theta <- c(0, -cumsum(pw2$TE))
V <- matrix(0, 4, 4)
for (i in 1:3)
  for (j in (i + 1):4)
    V[i, j] <- V[j, i] <- sum(pw2$seTE[seq.int(i, j - 1)]^2)
#
prob <- pnorm(-outer(theta, theta, "-") / sqrt(V))
diag(prob) <- NA_real_
#
expected.sucra <- rowMeans(prob, na.rm = TRUE)

#
# Tests
#

#
# 1) Three independent two-arm studies on a tree
#

test_that("sampling covariance agrees with sums of independent edge variances", {
  nma1 <- netmeta(MD, seMD, trt1, trt2, study, data = dat1,
                  sm = "MD", tau.preset = 0.2)
  #
  seen <- NULL
  #
  local_mocked_bindings(
    rmvnorm = function(n, mean, sigma) {
      seen <<- sigma
      matrix(rep(mean, each = n), nrow = n)
    },
    .package = "netmeta")
  #
  for (model in c("common", "random")) {
    variances <- dat1$seMD^2 + if (model == "random") 0.2^2 else 0
    expected <- matrix(0, 4, 4)
    #
    for (i in 2:4)
      for (j in 2:4)
        expected[i, j] <- sum(variances[seq_len(min(i, j) - 1)])
    #
    samples_netmeta(nma1, 2, model)
    expect_lt(max(abs(seen - expected)), 1e-10)
  }
})

#
# 2) Analytical SUCRA from pairwise normal probabilities
#

test_that("sampled SUCRA agrees with analytical pairwise probabilities", {
  nma2 <- netmeta(pw2, random = FALSE)
  #
  set.seed(284)
  #
  actual <- netrank(nma2, method = "SUCRA", nsim = 100000)$ranking.common
  #
  expect_lt(max(abs(actual - expected.sucra)), 0.006)
  expect_identical(names(which.max(actual)), "D")
})

#
# 3) Two-treatment network with singular anchored covariance
#

test_that("sampling supports a singular anchored covariance for two treatments", {
  nma3 <- netmeta(0.2, 0.1, "A", "B", studlab = 1, sm = "MD", random = FALSE)
  #
  seen <- NULL
  #
  local_mocked_bindings(
    rmvnorm = function(n, mean, sigma) {
      seen <<- sigma
      matrix(rep(mean, each = n), nrow = n)
    },
    .package = "netmeta")
  #
  samples_netmeta(nma3, 2, "common")
  expect_lt(max(abs(seen - diag(c(0, 0.01)))), 1e-12)
})
