library(devtools)
load_all()
library(testthat)
library(dplyr)

data1 <- read.csv("NMA_data_binary_FE.csv")
source("test_H_and_B_matrices.R")

pw1 <- pairwise(treat = T, event = R, n = N,
   studlab = Study, data = data1, sm = "OR")
 
net1 <- netmeta(pw1)
cm <- netcontrib(net1)
tomsTest <- runTomsHandBTest(data1)

test_that("Study contributions are not given by default ", {
  cm <- netcontrib(net1)
  rscm <- cm$study.random
  expect_null(rscm)
})

test_that("Study contributions sum to 1 ",{
  cm <- netcontrib(net1, path = TRUE, study = TRUE)
  scm <- cm$study.common
  sc <- subset(scm, comparison == "A:B") %>% select(study, contribution)
  expect_equal(sum(sc$contribution), 1)
})

test_that("Tom's weights = netmeta's weights",{
  tryCatch(
   expect_true(
     all.equal(unname(tomsTest$weights), unname(net1$w.common))
   ),
    error = function(e) {
      fail(paste(tomsTest$weights, net1$w.common ))
    })
})

test_that("Study contributions replicate Woods example with Tom's",{
  cm <- netcontrib(net1, path = TRUE, study = TRUE)
  scm <- cm$study.common
  tomsContrs <- tomsTest$contrs
  tc <- data.frame(names(tomsContrs),c(tomsContrs))
  sc <- subset(scm, comparison == "A:B") %>% select(study, contribution)
  expect_equal(tc[, 2], sc[, 2])
})
