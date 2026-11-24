#
# Setup
#

nma <- netmeta(c(0, 1, 0), rep(1, 3), rep("A", 3),
                c("B+C", "A+C", "C+D"), studlab = 1:3, sm = "MD")
cnma <- suppressWarnings(netcomb(nma))
#
complex <- netcomplex(cnma, "A + B + C")
comparison <- netcomparison(cnma, "A", "B + D")

#
# Tests
#

test_that("netcomplex uses lab.NA for all unavailable result columns", {
  for (common in c(TRUE, FALSE)) {
    expect_identical(
      capture.output(print(complex, common = common, random = !common)),
      capture.output(print(complex, common = common, random = !common,
                           lab.NA = gs("lab.NA"))))
    #
    for (label in c("NA", ".", "missing")) {
      output <- capture.output(print(complex, common = common, random = !common,
                                     lab.NA = label))
      row <- grep("A + B + C", output, fixed = TRUE, value = TRUE)
      #
      expect_length(row, 1)
      values <- strsplit(trimws(row), "[[:space:]]+")[[1]]
      expect_identical(tail(values, 4), rep(label, 4))
    }
  }
})

test_that("netcomparison uses lab.NA for all unavailable result columns", {
  for (common in c(TRUE, FALSE)) {
    expect_identical(
      capture.output(print(comparison, common = common, random = !common)),
      capture.output(print(comparison, common = common, random = !common,
                           lab.NA = gs("lab.NA"))))
    #
    for (label in c("NA", ".", "missing")) {
      output <- capture.output(print(comparison, common = common, random = !common,
                                     lab.NA = label))
      row <- grep("A:B + D", output, fixed = TRUE, value = TRUE)
      #
      expect_length(row, 1)
      values <- strsplit(trimws(row), "[[:space:]]+")[[1]]
      expect_identical(tail(values, 4), rep(label, 4))
    }
  }
})
