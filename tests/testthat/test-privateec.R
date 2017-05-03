# test-privateec.R - Bill White - 4/26/17
#
# Test the privateEC privateEC function on simulated data.

library(privateEC)
context("Private EC")

test_that("privateEC returns sane results", {
  num.samples <- 100
  num.variables <- 100
  pb <- 0.1
  nbias <- pb * num.variables
  signals <- sprintf("gene%04d", 1:nbias)
  sim.data <- createSimulation(d=num.variables, n=num.samples,
                               type="sva", verbose=FALSE)
  pec.results <- privateEC(data.sets=sim.data, is.simulated=TRUE, n=num.samples,
                           signal.names=signals, verbose=FALSE)
  expect_equal(ncol(pec.results$plots.data), 5)
  expect_equal(ncol(pec.results$melted.data), 4)
  expect_equal(length(pec.results$correct), nrow(pec.results$plots.data))
})
