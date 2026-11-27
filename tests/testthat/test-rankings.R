test_that("rankings preserves all rank columns when treatments take disjoint rank subsets", {
  set.seed(1)
  # Two well-separated pairs: A, B only take ranks 1-2; C, D only take ranks 3-4
  samples <- cbind(
    A = rnorm(1000, -10.1, 0.5),
    B = rnorm(1000, -9.9, 0.5),
    C = rnorm(1000, 9.9, 0.5),
    D = rnorm(1000, 10.1, 0.5)
  )

  rnk <- rankogram(samples, small.values = "desirable")

  expect_identical(dim(rnk$ranking.matrix.random), c(4L, 4L))
  expect_identical(colnames(rnk$ranking.matrix.random), as.character(1:4))
  expect_equal(unname(rowSums(rnk$ranking.matrix.random)), rep(1, 4), tolerance = 1e-12)
  expect_equal(unname(colSums(rnk$ranking.matrix.random)), rep(1, 4), tolerance = 1e-12)

  # C and D never reach rank 1 or 2
  expect_equal(unname(rnk$ranking.matrix.random[c("C", "D"), 1:2]), matrix(0, 2, 2))
  # A and B never reach rank 3 or 4
  expect_equal(unname(rnk$ranking.matrix.random[c("A", "B"), 3:4]), matrix(0, 2, 2))

  # SUCRA satisfies Salanti et al. (2011) identity: SUCRA_i = (a - mean_rank_i) / (a - 1)
  expected_sucra <- (4 - rnk$meanranks.random) / 3
  expect_equal(rnk$ranking.random, expected_sucra, tolerance = 1e-12)
  expect_gt(min(rnk$ranking.random[c("A", "B")]), 0.7)
  expect_lt(max(rnk$ranking.random[c("C", "D")]), 0.3)
})

test_that("rankogram handles deterministic rankings where each treatment achieves one rank", {
  set.seed(2)
  # Three well-separated treatments: A < B < C in every draw
  samples <- cbind(
    A = rnorm(200, -10, 0.1),
    B = rnorm(200, 0, 0.1),
    C = rnorm(200, 10, 0.1)
  )

  rnk <- rankogram(samples, small.values = "desirable")

  expect_identical(dim(rnk$ranking.matrix.random), c(3L, 3L))
  expect_equal(unname(rnk$ranking.matrix.random), diag(3))
  expect_equal(rnk$ranking.random, c(A = 1, B = 0.5, C = 0))
  expect_equal(rnk$meanranks.random, c(A = 1, B = 2, C = 3))
  expect_equal(rnk$medianranks.random, c(A = 1, B = 2, C = 3))
})

test_that("rankogram and netrank work for two-treatment networks and separated multi-treatment networks", {
  # Two-treatment network
  nma2 <- netmeta(0.2, 0.1, "A", "B", studlab = 1, sm = "MD", random = FALSE)
  set.seed(3)
  rnk2 <- rankogram(nma2, nsim = 500, small.values = "desirable")
  expect_identical(dim(rnk2$ranking.matrix.common), c(2L, 2L))
  expect_equal(rnk2$ranking.common, rnk2$ranking.matrix.common[, 1], tolerance = 1e-12)
  expect_equal(sum(rnk2$ranking.common), 1, tolerance = 1e-12)

  # Four-treatment network with two well-separated pairs (A,B vs C,D)
  dat4 <- data.frame(
    TE = c(-0.2, -20, -0.2),
    seTE = c(0.5, 0.5, 0.5),
    treat1 = c("A", "B", "C"),
    treat2 = c("B", "C", "D"),
    studlab = 1:3
  )
  nma4 <- netmeta(TE, seTE, treat1, treat2, studlab, data = dat4, sm = "MD", random = FALSE)
  set.seed(4)
  nr4 <- netrank(nma4, method = "SUCRA", nsim = 1000, small.values = "desirable")
  expect_gt(min(nr4$ranking.common[c("A", "B")]), 0.7)
  expect_lt(max(nr4$ranking.common[c("C", "D")]), 0.3)
})
