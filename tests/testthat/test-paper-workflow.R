# test-paper-workflow.R - Bill White - 4/27/17
#
# Test replicating the Bioinformatics paper simulated and real data analysis.

library(privateEC)
context("Paper Workflows")

test_that("run one workflow of a simulation plus an analysis step", {
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
  expect_equal(length(one.step.result), 2)
  expect_equal(length(one.step.result$run.results), 3)
})

test_that("run one workflow for a real data analysis", {
  data(fullfMRI2)
  data(phenos)
  # only 100 variables for a test
  real.result <- paperRealWorkflow(corr.mat=fullfMRI2[, 2900:ncol(fullfMRI2)],
                                   phenos=phenos,
                                   update.freq=50,
                                   verbose=FALSE)
  expect_equal(length(real.result), 2)
  expect_equal(length(real.result$run.results), 3)
})
