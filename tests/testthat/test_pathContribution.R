library(devtools)
library(testthat)
load_all()

data("Woods2010")

pw1 <- pairwise(treat = treatment, studlab = author,
                event = r, n = N,
                data = Woods2010, sm = "OR")

net1 <- netmeta(pw1)

test_that("Path contributions are not given by default ", {
  cm <- netcontrib(net1)
  rpcm <- cm$path.random
  expect_null(rpcm)
})

test_that("Path contributions of Placebo:Salmeterol comparison of Woods2010 sum up to 1", {
  cm <- netcontrib(net1, path = TRUE)
  # Path contributions of common effects model
  cpcm <- cm$path.common
  expect_equal(subset(cpcm, comparison == "Placebo:Salmeterol")$contribution %>%
                 sum(), 1)
})
