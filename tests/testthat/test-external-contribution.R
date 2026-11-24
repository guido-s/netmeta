#
# Setup
#

dat <- read.csv(testthat::test_path("data_binary.csv"))

pw <- pairwise(treat = T, event = R, n = N,
               studlab = Study, data = dat, sm = "OR")
#
nma <- netmeta(pw, random = FALSE)
#
cm <- netcontrib(nma)
# The current version of CINeMA uses the old hat matrix
# cm <- netcontrib(nma, hatmatrix.F1000 = TRUE)
#
tomsTest <- BigDataStatMeth_results(dat, pw, cm)

#
# Tests
#

test_that("Study contributions are not given by default ", {
  cm <- netcontrib(nma)
  rscm <- cm$study.random
  expect_null(rscm)
})

test_that("Study contributions sum to 1 ",{
  cm <- netcontrib(nma, path = TRUE, study = TRUE)
  scm <- cm$study.common
  sc <- subset(scm, comparison == "A:B") |> dplyr::select(study, contribution)
  expect_equal(sum(sc$contribution), 1)
})

test_that("Tom's weights equal netmeta's weights", {
  o <- paste(nma$studlab, nma$treat1, nma$treat2, sep = "_")
  #
  expect_setequal(names(tomsTest$weights), o)
  #
  expect_equal(
    unname(tomsTest$weights[o]),
    unname(nma$w.common)
  )
})

test_that("Study contributions replicate Woods example with Tom's",{
  cm <- netcontrib(nma, path = TRUE, study = TRUE)
  scm <- cm$study.common
  tomsContrs <- tomsTest$contrs
  tc <- data.frame(names(tomsContrs),c(tomsContrs))
  sc <- subset(scm, comparison == "A:B") |> dplyr::select(study, contribution)
  #
  expect_equal(tc[, 2], sc[, 2])
})
