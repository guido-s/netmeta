#
# Setup
#

data(smokingcessation, package = "netmeta")
#
pw <-
  pairwise(
    list(treat1, treat2, treat3),
    event = list(event1, event2, event3),
    n = list(n1, n2, n3),
    data = smokingcessation, sm = "OR")
#
nma_list <-
  list(
    classic = netmeta(pw, method.random.ci = "classic"),
    "t-dist" = netmeta(pw, method.random.ci = "t-dist"))
#
split_list <- lapply(nma_list, netsplit)
#
nmab <- netmetabin(
  c(10, 12, 9), rep(100, 3), c(15, 14, 13), rep(100, 3),
  c("A", "A", "B"), c("B", "C", "C"), studlab = 1:3)

#
# Tests
#

test_that("stored inference respects the model and random-effects CI method", {
  expect_gt(nma_list$classic$df.Q, 0)
  #
  for (method in names(nma_list)) {
    nma <- nma_list[[method]]
    for (model in c("common", "random")) {
      df <- if (model == "random" && method == "t-dist") nma$df.Q else NULL
      #
      for (part in c("", "direct.", "indirect.")) {
        expected <- ci(
          nma[[paste0("TE.", part, model)]],
          nma[[paste0("seTE.", part, model)]],
          level = nma$level.ma, df = df)
        #
        expect_equal(nma[[paste0("lower.", part, model)]], expected$lower)
        expect_equal(nma[[paste0("upper.", part, model)]], expected$upper)
        expect_equal(nma[[paste0("pval.", part, model)]], expected$p)
      }
    }
  }
  #
  expect_equal(nma_list$classic$TE.random, nma_list[["t-dist"]]$TE.random)
  expect_equal(nma_list$classic$seTE.random, nma_list[["t-dist"]]$seTE.random)
  expect_gt(
    max(nma_list$classic$lower.random - nma_list[["t-dist"]]$lower.random,
        na.rm = TRUE),
    0)
})

test_that("pooled intervals use level.ma independently of level", {
  for (method in names(nma_list)) {
    nma <- netmeta(pw, method.random.ci = method, level = 0.90, level.ma = 0.95)
    #
    for (model in c("common", "random")) {
      df <- if (model == "random" && method == "t-dist") nma$df.Q else NULL
      #
      for (part in c("", "direct.", "indirect.")) {
        expected <- ci(
          nma[[paste0("TE.", part, model)]],
          nma[[paste0("seTE.", part, model)]],
          level = nma$level.ma, df = df)
        #
        expect_equal(nma[[paste0("lower.", part, model)]], expected$lower)
        expect_equal(nma[[paste0("upper.", part, model)]], expected$upper)
      }
    }
  }
})

test_that("split inference respects the model and random-effects CI method", {
  for (method in names(nma_list)) {
    nma <- nma_list[[method]]
    split <- split_list[[method]]
    #
    for (model in c("common", "random")) {
      df <- if (model == "random" && method == "t-dist") nma$df.Q else NULL
      #
      for (name in c(model, paste0(c("direct.", "indirect.", "compare."), model))) {
        dat <- split[[name]]
        expected <- ci(dat$TE, dat$seTE, level = nma$level.ma, df = df)
        #
        expect_equal(dat$lower, expected$lower)
        expect_equal(dat$upper, expected$upper)
        expect_equal(dat$p, expected$p)
      }
    }
  }
})

test_that("network and split print headings distinguish z and t statistics", {
  for (method in names(nma_list)) {
    for (model in c("common", "random")) {
      statistic <- if (model == "random" && method == "t-dist") "t" else "z"
      #
      output <- capture.output(print(
        nma_list[[method]], common = model == "common", random = model == "random",
        reference.group = "A"))
      header <- output[grepl("p-value", output, fixed = TRUE) &
                         grepl("95% CI", output, fixed = TRUE)]
      #
      expect_length(header, 1)
      expect_match(header, paste0("[[:space:]]", statistic, "[[:space:]]+p-value"))
      #
      output <- capture.output(print(
        split_list[[method]], common = model == "common", random = model == "random",
        ci = TRUE))
      header <- output[grepl("comparison", output, fixed = TRUE) &
                         grepl("p-value", output, fixed = TRUE)]
      #
      expect_length(header, 1)
      expect_match(header, paste0("[[:space:]]", statistic, "[[:space:]]+p-value"))
    }
  }
})

test_that("table disagreement p-values agree with split inference", {
  for (method in names(nma_list)) {
    nma <- nma_list[[method]]
    di <- direct.indirect(nma)
    trts <- comptrts(nma, TRUE, nma$reference.group, nma$baseline.reference,
                     NULL, nma$sep.trts, "")
    #
    for (model in c("common", "random")) {
      table <- mat2dat.table(
        di, model, trts, backtransf = FALSE,
        digits = 4, digits.I2 = 1, digits.pval = 8,
        scientific.pval = FALSE, zero.pval = TRUE, JAMA.pval = FALSE,
        big.mark = "", text.NA = ".", writexl = FALSE)
      split <- mat2dat.split(di, model, trts)
      #
      expected <- formatPT(
        split$compare$p, digits = 8, scientific = FALSE,
        zero = TRUE, JAMA = FALSE, lab.NA = ".")
      #
      expect_identical(table$Incoherence, expected)
      expect_equal(nrow(table), nrow(trts))
    }
  }
})

test_that("plot intervals retain the CI method at a different confidence level", {
  for (method in names(nma_list)) {
    nma <- nma_list[[method]]
    split <- split_list[[method]]
    #
    for (model in c("common", "random")) {
      dat <- plot(split, pooled = model, type = "Bland-Altman", level = 0.90)$data
      reference <- split[[paste0("compare.", model)]]
      index <- match(dat$comparison, reference$comparison)
      df <- if (model == "random" && method == "t-dist") nma$df.Q else NULL
      #
      expected <- ci(reference$TE, reference$seTE, level = 0.90, df = df)
      #
      expect_false(anyNA(index))
      expect_equal(dat$lower, expected$lower[index])
      expect_equal(dat$upper, expected$upper[index])
      expect_no_error(plot(split, pooled = model, type = "Wilson", level = 0.90))
    }
  }
})

test_that("binary network objects provide the classic CI default", {
  expect_identical(nmab$method.random.ci, "classic")
  expect_no_error(netsplit(nmab))
})

test_that("older objects acquire CI defaults without losing existing settings", {
  for (object in list(nma_list$classic, summary(nma_list$classic), split_list$classic, nmab)) {
    object$version <- "3.7-0"
    nested <- !inherits(object, "netmeta")
    #
    if (nested)
      object$x$method.random.ci <- NULL
    else
      object$method.random.ci <- NULL
    #
    updated <- updateversion(object)
    #
    expect_identical(
      if (nested) updated$x$method.random.ci else updated$method.random.ci,
      "classic")
    #
    if (nested)
      object$x$method.random.ci <- "t-dist"
    else
      object$method.random.ci <- "t-dist"
    updated <- updateversion(object)
    #
    expect_identical(
      if (nested) updated$x$method.random.ci else updated$method.random.ci,
      "t-dist")
  }
  #
  old.binary <- nmab
  old.binary$version <- "3.7-0"
  old.binary$method.random.ci <- NULL
  #
  expect_no_error(netsplit(old.binary))
  #
  old.split <- split_list$classic
  old.split$version <- "3.7-0"
  old.split$x$method.random.ci <- NULL
  #
  expect_no_error(capture.output(print(old.split)))
})

test_that("t-based inference falls back to normal inference with zero residual df", {
  nma <- netmeta(c(0.2, 0.3), c(0.1, 0.2), c("A", "B"), c("B", "C"),
                  studlab = 1:2, method.random.ci = "t-dist")
  #
  expect_equal(nma$df.Q, 0, tolerance = 0)
  expect_true(is.na(nma$tau))
  #
  for (part in c("", "direct.")) {
    for (field in c("TE.", "seTE.", "lower.", "upper.", "pval.")) {
      expect_equal(
        nma[[paste0(field, part, "random")]],
        nma[[paste0(field, part, "common")]])
    }
  }
  #
  expect_no_error(netsplit(nma))
})
