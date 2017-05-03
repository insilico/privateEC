# test-replicate-sims.R - Bill White - 4/27/17
#
# Test replicating the Bioinformatics paper.

library(privateEC)
context("Private EC Simulated Data Analysis")

test_that("run one workflow step", {
  num.samples <- 100
  num.variables <- 100
  pb <- 0.1
  bias <- 0.4
  one.step.result <- paperSimWorkflow(myrun="001",
                                      n=num.samples,
                                      d=num.variables,
                                      pb=pb,
                                      update.freq=50,
                                      verbose=FALSE)
  expect_equal(length(one.step.result$run.results), 3)
})
