suppressPackageStartupMessages(library(netmeta))

d <- data.frame(TE = c(.6, .6, 0, .6, .6, 0, -.4, -.4, 0),
  se = sqrt(c(.091, .091, .002, .091, .091, .002, .011, .011, .002)),
  t1 = rep(c("A", "A", "B"), 3), t2 = rep(c("B", "C", "C"), 3),
  study = rep(c(20, 100, 3), each = 3))
Y <- rbind(c(-.6, -.6), c(-.6, -.6), c(.4, .4))
# Known shared-arm variances; independently profile the Gaussian likelihood.
S <- lapply(c(.09, .09, .01), function(v) diag(.001, 2) + v)
T <- matrix(c(1, .5, .5, 1), 2)
profile <- function(tau2) {
  V <- lapply(S, function(s) s + tau2 * T)
  W <- lapply(V, solve)
  info <- Reduce("+", W)
  beta <- solve(info, Reduce("+", lapply(1:3, function(i) W[[i]] %*% Y[i, ])))
  logdet <- function(v) as.numeric(determinant(v, logarithm = TRUE)$modulus)
  rss <- sum(vapply(1:3, function(i) {
    r <- Y[i, ] - beta
    as.numeric(t(r) %*% W[[i]] %*% r)
  }, numeric(1)))
  sum(vapply(V, logdet, numeric(1))) + rss + logdet(info)
}
optimum <- optimize(profile, c(0, 2), tol = 1e-11)
tau2 <- if (profile(0) < optimum$objective) 0 else optimum$minimum
for (ref in c("A", "B", "C")) for (reverse.rows in c(FALSE, TRUE)) {
  z <- if (reverse.rows) d[rev(seq_len(nrow(d))), ] else d
  fit <- netmeta(TE, se, t1, t2, study, data = z, sm = "HR",
                  method.tau = "REML", reference.group = ref)
  stopifnot(abs(fit$tau - sqrt(tau2)) < 2e-5)
}
# ML uses the same sampling covariance; do not conflate covariance repair
# with the optimizer's choice of a local versus boundary likelihood optimum.
ml <- netmeta(TE, se, t1, t2, study, data = d, sm = "HR", method.tau = "ML",
                reference.group = "A", keeprma = TRUE)
study.order <- match(unique(ml$rma.tau$data$studlab), unique(d$study))
stopifnot(max(abs(ml$rma.tau$V - metafor::bldiag(S[study.order]))) < 1e-12)
# Four correlated arms require unequal off-diagonal sampling covariances.
S4 <- matrix(c(.2, .03, .01, .03, .3, .025, .01, .025, .4), 3)
C4 <- matrix(0, 4, 4); C4[-1, -1] <- S4
variance <- outer(diag(C4), diag(C4), "+") - 2 * C4
pairs <- t(combn(1:4, 2))
means <- rbind(c(0, .2, .6, .3), c(0, .8, .3, -.1), c(0, -.2, .1, .7))
d4 <- do.call(rbind, lapply(1:3, function(s)
  data.frame(TE = means[s, pairs[, 1]] - means[s, pairs[, 2]],
             se = sqrt(variance[pairs]), t1 = LETTERS[pairs[, 1]],
             t2 = LETTERS[pairs[, 2]], study = s)))
fit <- netmeta(TE, se, t1, t2, study, data = d4, correlated = TRUE,
                sm = "MD", method.tau = "REML", reference.group = "A", keeprma = TRUE)
stopifnot(max(abs(fit$rma.tau$V - kronecker(diag(3), S4))) < 1e-12)
# A single two-arm study still returns non-estimable heterogeneity.
one <- suppressWarnings(netmeta(.2, .1, "A", "B", studlab = 1, sm = "MD",
                                 method.tau = "REML"))
stopifnot(is.na(one$tau))
