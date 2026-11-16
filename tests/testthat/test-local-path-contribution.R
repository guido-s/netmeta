#
# Setup
#

data(Woods2010, package = "metabook")

pw <- pairwise(treat = treatment, event = r, n = N,
                studlab = author, data = Woods2010, sm = "OR")
#
nma <- netmeta(pw)

#
# Tests
#

test_that("Path contributions are not given by default ", {
  cm <- netcontrib(nma)
  rpcm <- cm$path.random
  expect_null(rpcm)
})

test_that("Path contributions of Placebo:Salmeterol comparison of Woods2010 sum up to 1", {
  cm <- netcontrib(nma, path = TRUE)
  # Path contributions of common effects model
  cpcm <- cm$path.common
  expect_equal(
    subset(cpcm, comparison == "Placebo:Salmeterol")$contribution |> sum(),
    1)
})
