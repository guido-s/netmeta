test_that("single network treatments use sorted order", {
  nc <- netconnection(
    treat1 = c("B", "A", "A", "B", "D"),
    treat2 = c("E", "D", "C", "D", "E"),
    studlab = paste0("study", 1:5)
  )

  expect_equal(nc$n.subnets, 1)
  expect_identical(nc$seq, sort(nc$seq))
  expect_identical(rownames(nc$A.matrix), nc$seq)
})

test_that("disconnected networks retain subnetwork grouping", {
  nc <- netconnection(
    treat1 = c("A", "C"),
    treat2 = c("B", "D"),
    studlab = c("study1", "study2")
  )

  expect_equal(nc$n.subnets, 2)
  expect_equal(nc$seq, rownames(nc$A.matrix))
})
