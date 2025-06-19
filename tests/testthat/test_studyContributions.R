#
library(devtools)
load_all()
library(testthat)
library(dplyr)

data1 <- read.csv("NMA_data_binary_FE.csv")
source("Test H and B matrices.R")

p1 <- pairwise(treat = T, event = R, n = N,
   studlab = Study, data = data1, sm = "OR")
 
net1 <- netmeta(p1)
cm <- netcontrib(net1)
tomsTest <- runTomsHandBTest(data1)

test_that("studyContibutions are not given by default ",{
  cm <- netcontrib(net1)
  rscm <- cm$random.scm
  expect_null(rscm)
})

test_that("studyContibutions sum to 1 ",{
  cm <- netcontrib(net1, pathContribution=T, studyContribution=T)
  scm <- cm$common.scm
  sc <- subset(scm, comparison=="A:B")[,c('study','x')]
  expect_equal(sum(sc$x),1)
})

test_that("studyContribution weights = netmeta's weights",{
  studyWeights <- readRDS("studyWeights.rds")
  tryCatch(
   expect_equal(studyWeights, net1$w.common),
    error = function(e) {
      fail(paste(studyWeights, net1$w.common ))
    })
})

test_that("Tom's weights = netmeta's weights",{
  tryCatch(
   expect_true(
     all.equal(unname(tomsTest$weights),unname(net1$w.common))
   ),
    error = function(e) {
      fail(paste(tomsTest$weights, net1$w.common ))
    })
})

test_that("studyContibutions replicate woods example with Tom's",{
  cm <- netcontrib(net1, pathContribution=T, studyContribution=T)
  scm <- cm$common.scm
  tomsContrs <- tomsTest$contrs
  tc <- data.frame(names(tomsContrs),c(tomsContrs))
  sc <- subset(scm, comparison=="A:B")[,c('study','x')]
  # # print(tc)
  # print(sc)
  expect_equal(tc[,2],sc[,2])
})
