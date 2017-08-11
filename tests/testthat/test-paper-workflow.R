# test-paper-workflow.R - Bill White - 4/27/17
#
# Test replicating the Bioinformatics paper simulated and real data analysis.

library(privateEC)
context("Paper Workflows")

test_that("run one workflow of a simulation plus an analysis step", {
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  upd.frq <- 0.1 * num.variables
  one.step.result <- paperSimWorkflow(n.samples = num.samples,
                                      n.variables = num.variables,
                                      pct.signals = pct.signals,
                                      update.freq = upd.frq,
                                      verbose = FALSE)
  expect_equal(length(one.step.result), 2)
  expect_equal(length(one.step.result$run.results), 3)
})

test_that("run one workflow for a real data analysis", {
  data(rsfMRIcorrMDD)
  # ~100 variables for a test
  test.mat <- rsfMRIcorrMDD[, 2900:ncol(rsfMRIcorrMDD)]
  real.result <- paperRealWorkflow(real.data = test.mat,
                                   label = "phenos",
                                   update.freq = 5,
                                   verbose = FALSE)
  expect_equal(length(real.result), 2)
  expect_equal(length(real.result$run.results), 3)
})
