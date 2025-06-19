#
library(devtools)
library(testthat)
load_all()

data("Woods2010")
data1 <- Woods2010

p1 <- pairwise(treat=treatment,
               studlab=author,
               event=r,
               n=N,
               data = data1, sm = "OR")
 
net1 <- netmeta(p1)


test_that("PathContibutions are not given by default ",{
  cm <- netcontrib(net1)
  rpcm <- cm$random.pcm
  expect_null(rpcm)
})

test_that("Path contributions of Placebo:Salmeterol comparison of Woods2010 sum up to 1",{
  cm <- netcontrib(net1, pathContribution=T)
  #Path contributions of common effects model
  cpcm <- cm$common.pcm
  print(cpcm)
  expect_equal(subset(cpcm,comparison=="Placebo:Salmeterol")$contribution %>% sum(),1)
})
