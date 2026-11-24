#
# Setup
#

#
# 1) Data set with three three-arm studies
#

dat1 <- data.frame(
  TE = c(0.6, 0.6, 0, 0.6, 0.6, 0, -0.4, -0.4, 0),
  se = sqrt(c(0.091, 0.091, 0.002, 0.091, 0.091, 0.002,
              0.011, 0.011, 0.002)),
  t1 = rep(c("A", "A", "B"), 3),
  t2 = rep(c("B", "C", "C"), 3),
  study = rep(c(20, 100, 3), each = 3))

# Known shared-arm variances; independently profile the Gaussian likelihood
#
Y <- rbind(c(-0.6, -0.6), c(-0.6, -0.6), c(0.4, 0.4))
S <- lapply(c(0.09, 0.09, 0.01), function(v) diag(0.001, 2) + v)
T <- matrix(c(1, 0.5, 0.5, 1), 2)
#
profile <- function(x) {
  V <- lapply(S, function(s) s + x * T)
  #
  W <- lapply(V, solve)
  info <- Reduce("+", W)
  #
  beta <- solve(info, Reduce("+", lapply(1:3, function(i) W[[i]] %*% Y[i, ])))
  #
  logdet <- function(v) as.numeric(determinant(v, logarithm = TRUE)$modulus)
  #
  rss <- sum(vapply(1:3, function(i) {
    r <- Y[i, ] - beta
    as.numeric(t(r) %*% W[[i]] %*% r)
  }, numeric(1)))
  sum(vapply(V, logdet, numeric(1))) + rss + logdet(info)
}
#
optimum <- optimize(profile, c(0, 2), tol = 1e-11)
#
tau2.opt <- if (profile(0) < optimum$objective) 0 else optimum$minimum

#
# 2) Four correlated arms require unequal off-diagonal sampling covariances
#

S4 <- matrix(c(0.2, 0.03, 0.01, 0.03, 0.3, 0.025, 0.01, 0.025, 0.4), 3)
C4 <- matrix(0, 4, 4)
C4[-1, -1] <- S4
#
variance <- outer(diag(C4), diag(C4), "+") - 2 * C4
#
pairs <- t(combn(1:4, 2))
means <- rbind(c(0, 0.2, 0.6, 0.3), c(0, 0.8, 0.3, -0.1),
               c(0, -0.2, 0.1, 0.7))
#
dat2 <- do.call(rbind, lapply(1:3, function(s)
  data.frame(TE = means[s, pairs[, 1]] - means[s, pairs[, 2]],
             se = sqrt(variance[pairs]), t1 = LETTERS[pairs[, 1]],
             t2 = LETTERS[pairs[, 2]], study = s)))


#
# Tests
#

#
# 1) Comparison with the profiled restricted likelihood
#

test_that("REML agrees with the profiled likelihood across references and row orders", {
  for (ref in c("A", "B", "C")) {
    for (reverse.rows in c(FALSE, TRUE)) {
      if (reverse.rows)
        z <- dat1[rev(seq_len(nrow(dat1))), ]
      else
        z <- dat1
      #
      nma1.ref <- netmeta(TE, se, t1, t2, study, data = z, sm = "HR",
                          method.tau = "REML", reference.group = ref)
      #
      expect_lt(abs(nma1.ref$tau - sqrt(tau2.opt)), 2e-5)
    }
  }
})

#
# 2) Sampling covariance for maximum likelihood estimation
#

test_that("ML and REML treat numeric treatment labels as names", {
  dat.numeric <- dat1
  dat.numeric$t1 <- 10 * match(dat1$t1, c("A", "B", "C"))
  dat.numeric$t2 <- 10 * match(dat1$t2, c("A", "B", "C"))
  dat.character <- dat.numeric
  dat.character$t1 <- as.character(dat.character$t1)
  dat.character$t2 <- as.character(dat.character$t2)
  #
  for (method in c("ML", "REML")) {
    fit.numeric <- netmeta(
      TE, se, t1, t2, study, data = dat.numeric, sm = "HR",
      method.tau = method, reference.group = "10", keeprma = TRUE)
    fit.character <- netmeta(
      TE, se, t1, t2, study, data = dat.character, sm = "HR",
      method.tau = method, reference.group = "10", keeprma = TRUE)
    #
    expect_equal(fit.numeric$tau2, fit.character$tau2)
    expect_equal(fit.numeric$TE.random, fit.character$TE.random)
    expect_equal(fit.numeric$seTE.random, fit.character$seTE.random)
    expect_equal(fit.numeric$rma.tau$V, fit.character$rma.tau$V)
  }
})

test_that("ML uses the sampling covariance of the retained contrasts", {
  # Check the covariance, independently of the optimizer's choice of optimum.
  nma1.ml <- netmeta(TE, se, t1, t2, study, data = dat1, sm = "HR",
                     method.tau = "ML", reference.group = "A", keeprma = TRUE)
  #
  study.order <- match(unique(nma1.ml$rma.tau$data$studlab), unique(dat1$study))
  #
  expect_lt(
    max(abs(nma1.ml$rma.tau$V - metafor::bldiag(S[study.order]))),
    1e-12)
})

#
# 3) Network meta-analysis with four correlated arms
#

test_that("REML retains unequal sampling covariances for correlated arms", {
  nma2 <- netmeta(TE, se, t1, t2, study, data = dat2, correlated = TRUE,
                  sm = "MD", method.tau = "REML", reference.group = "A",
                  keeprma = TRUE)
  #
  expect_lt(max(abs(nma2$rma.tau$V - kronecker(diag(3), S4))), 1e-12)
})

#
# 4) Non-estimable heterogeneity for a single two-arm study
#

test_that("REML heterogeneity is non-estimable for a single two-arm study", {
  suppressWarnings(
    nma3 <- netmeta(0.2, 0.1, "A", "B", studlab = 1, sm = "MD",
                    method.tau = "REML"))
  #
  expect_true(is.na(nma3$tau))
})
