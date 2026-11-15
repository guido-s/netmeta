suppressPackageStartupMessages(library(netmeta))

# For independent studies of one pair, precision loss has a closed form.
for (scale in c(.1, 1, 10)) {
  se <- c(.8, .4, .04, .1) * scale
  x <- netmeta(rep(.2, 4), se, rep("A", 4), rep("B", 4),
                 studlab = 1:4, sm = "MD", tau.preset = 0)
  impact <- netimpact(x)$impact.common[, 1]
  w <- 1 / se^2
  expected <- (w - 1 / (1000 * max(se))^2) / sum(w)
  stopifnot(max(abs(impact - expected)) < 1e-7)
}
# A singular three-arm covariance block plus two independent two-arm studies.
# Direct GLS uses just two independent contrasts from the three-arm study.
S <- diag(c(.04, .09, .012, 1e6))
S[1, 2] <- S[2, 1] <- .01
X <- rbind(c(1, 0), c(0, 1), c(1, 0), c(0, 1))
y <- c(.2, .5, .25, .3)
C <- solve(crossprod(X, solve(S, X)))
beta <- C %*% crossprod(X, solve(S, y))
d <- data.frame(TE = c(.2, .5, .3, .25, .3),
  se = sqrt(c(.04, .09, .11, .012, 1e6)),
  t1 = c("A", "A", "B", "A", "A"), t2 = c("B", "C", "C", "B", "C"),
  study = c("multi", "multi", "multi", "precise", "ignored"),
  correlated = c(TRUE, TRUE, TRUE, FALSE, FALSE))
for (rows in list(seq_len(nrow(d)), rev(seq_len(nrow(d))))) {
  fit <- netmeta(TE, se, t1, t2, study, correlated = correlated,
                  data = d[rows, ], sm = "MD", tau.preset = 0)
  stopifnot(abs(fit$TE.common["A", "B"] - beta[1]) < 1e-8,
            abs(fit$seTE.common["A", "B"] - sqrt(C[1, 1])) < 1e-8,
            abs(fit$seTE.common["A", "C"] - sqrt(C[2, 2])) < 1e-8)
}
