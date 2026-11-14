library(netmeta)

# Two three-arm studies sharing one treatment have four independent
# contrasts and five treatments, hence no residual degrees of freedom.
pw <- pairwise(treat = c("A", "B", "C", "C", "D", "E"),
               event = c(20, 40, 30, 30, 20, 40), n = rep(100, 6),
               studlab = rep(1:2, each = 3), sm = "OR")

fit_network <- function(data) {
  warnings <- character()
  fit <- withCallingHandlers(
    netmeta(data),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  stopifnot(length(warnings) == 0L)
  fit
}

check_saturated <- function(fit) {
  stopifnot(fit$df.Q == 0, is.na(fit$tau), is.na(fit$tau2),
            is.na(fit$pval.Q),
            isTRUE(all.equal(fit$TE.random, fit$TE.common)),
            isTRUE(all.equal(fit$seTE.random, fit$seTE.common)))
}

fit <- fit_network(pw)
check_saturated(fit)

# Relabeling treatments must not change an estimate or its uncertainty.
labels <- list(c("E", "D", "C", "B", "A"),
               c("B", "A", "D", "C", "E"))
for (new_labels in labels) {
  renamed <- pw
  renamed$treat1 <- new_labels[match(pw$treat1, LETTERS[1:5])]
  renamed$treat2 <- new_labels[match(pw$treat2, LETTERS[1:5])]
  other <- fit_network(renamed)
  check_saturated(other)
  stopifnot(
    isTRUE(all.equal(unname(fit$TE.random),
                     unname(other$TE.random[new_labels, new_labels]))),
    isTRUE(all.equal(unname(fit$seTE.random),
                     unname(other$seTE.random[new_labels, new_labels]))))
}
check_saturated(fit_network(pw[rev(seq_len(nrow(pw))), ]))

# A single five-arm study is also saturated. Its ten comparisons still
# contribute only four independent contrasts.
five_arm <- pairwise(treat = LETTERS[1:5], event = c(20, 25, 30, 35, 40),
                     n = rep(100, 5), studlab = rep("five", 5), sm = "OR")
check_saturated(fit_network(five_arm))

# Mixed two-, three-, and four-arm studies contribute 1 + 2 + 3 contrasts.
# With seven treatments connected at single articulation points, df = 0.
mixed <- pairwise(treat = c("A", "B", "B", "C", "D", "D", "E", "F", "G"),
                  event = c(20, 30, 30, 20, 40, 40, 25, 35, 45),
                  n = rep(100, 9),
                  studlab = c("two", "two", rep("three", 3), rep("four", 4)),
                  sm = "OR")
check_saturated(fit_network(mixed))

# Exercise both calculation engines, including the legacy internal entry
# point, with the same saturated multi-arm input.
prepared <- netmeta:::prepare(pw$TE, pw$seTE, pw$treat1, pw$treat2,
                             pw$studlab, func.inverse = invmat)
d <- prepared$data
args <- list(TE = d$TE, W = as.matrix(prepared$W),
             seTE = sqrt(1 / d$weights),
             treat1 = d$treat1, treat2 = d$treat2,
             treat1.pos = d$treat1.pos, treat2.pos = d$treat2.pos,
             narms = d$narms, studlab = d$studlab,
             sm = "OR", level = 0.95, level.ma = 0.95,
             seTE.orig = d$seTE, func.inverse = invmat)
for (engine in c("nma.ruecker", "nma_ruecker")) {
  engine_args <- args
  if (engine == "nma_ruecker")
    engine_args$Cov0 <- prepared$Cov
  result <- do.call(getFromNamespace(engine, "netmeta"), engine_args)
  stopifnot(result$df == 0, is.na(result$tau), is.na(result$pval.Q))
}

# Positive residual df must retain the usual DerSimonian-Laird estimate.
# With two treatments this has an independent pairwise-analysis oracle.
effects <- c(-0.2, -0.8, -1.4)
ses <- c(0.1, 0.2, 0.3)
positive <- netmeta(effects, ses, rep("A", 3), rep("B", 3),
                    studlab = paste0("study", 1:3), sm = "MD")
weights <- 1 / ses^2
pooled <- weighted.mean(effects, weights)
Q <- sum(weights * (effects - pooled)^2)
C <- sum(weights) - sum(weights^2) / sum(weights)
tau2 <- max(0, (Q - 2) / C)
stopifnot(positive$df.Q == 2,
          isTRUE(all.equal(positive$Q, Q)),
          isTRUE(all.equal(positive$tau2, tau2)))

# Closing one loop in the mixed-arm network adds exactly one residual df.
loop <- pairwise(treat = c("A", "G"), event = c(20, 45), n = c(100, 100),
                 studlab = c("loop", "loop"), sm = "OR")
positive_mixed <- fit_network(rbind(mixed, loop))
stopifnot(positive_mixed$df.Q == 1, is.finite(positive_mixed$tau))
