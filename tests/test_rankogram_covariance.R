suppressPackageStartupMessages(library(netmeta))

# Three independent two-arm studies on a tree; exact covariance follows
# from summing edge contrasts, not from a package ranking calculation.
TE <- c(.3, -.2, .15)
se <- c(.1, .5, .12)
x <- netmeta(TE, se, c("A", "B", "C"), c("B", "C", "D"),
             studlab = 1:3, sm = "MD", tau.preset = .2)
sampler <- getFromNamespace("samples_netmeta", "netmeta")
e <- new.env(parent = environment(sampler))
seen <- NULL
e$rmvnorm <- function(n, mean, sigma) {
  seen <<- sigma
  matrix(rep(mean, each = n), nrow = n)
}
environment(sampler) <- e
for (model in c("common", "random")) {
  variances <- se^2 + if (model == "random") .2^2 else 0
  expected <- matrix(0, 4, 4)
  for (i in 2:4) for (j in 2:4)
    expected[i, j] <- sum(variances[seq_len(min(i, j) - 1)])
  sampler(x, 2, model)
  stopifnot(max(abs(seen - expected)) < 1e-10)
}
# Public ranking API: limiting SUCRA is the mean of pairwise normal probabilities.
p <- pairwise(c("A", "B", "B", "C", "C", "D"),
               event = c(62, 34, 60, 66, 237, 208),
               n = c(100, 100, 100, 100, 500, 500),
               studlab = rep(1:3, each = 2), sm = "OR")
x <- netmeta(p, random = FALSE)
theta <- c(0, -cumsum(p$TE))
V <- matrix(0, 4, 4)
for (i in 1:3) for (j in (i + 1):4)
  V[i, j] <- V[j, i] <- sum(p$seTE[seq.int(i, j - 1)]^2)
prob <- pnorm(-outer(theta, theta, "-") / sqrt(V))
diag(prob) <- NA_real_
expected <- rowMeans(prob, na.rm = TRUE)
set.seed(284)
actual <- netrank(x, method = "SUCRA", nsim = 100000)$ranking.common
stopifnot(max(abs(actual - expected)) < .006,
          names(which.max(actual)) == "D")
# A two-treatment network remains supported (singular anchored covariance).
two <- netmeta(.2, .1, "A", "B", studlab = 1, sm = "MD", random = FALSE)
sampler(two, 2, "common")
stopifnot(max(abs(seen - diag(c(0, .01)))) < 1e-12)
