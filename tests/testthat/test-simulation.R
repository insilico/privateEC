# test-simulation.R - Bill White - 4/30/17
#
# Test the creation of simulated data.

library(privateEC)
context("Simulation")

test_that("createSimulation mainEffect returns sane results", {
  sim.type <- "mainEffect"
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  bias <- 0.4
  sim.data <- createSimulation(num.samples=num.samples,
                               num.variables=num.variables,
                               pct.signals=pct.signals,
                               bias=bias,
                               sim.type=sim.type,
                               verbose=FALSE)
  expect_equal(nrow(sim.data$train) + nrow(sim.data$holdout) + nrow(sim.data$validation),
               num.samples * 3)
  expect_equal(ncol(sim.data$train), num.variables + 1)
})

test_that("createSimulation Erdos-Renyi network base returns sane results", {
  sim.type <- "interactionErdos"
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  bias <- 0.4
  sim.data <- createSimulation(num.variables=num.variables,
                               num.samples=num.samples,
                               pct.signals=pct.signals,
                               bias=bias,
                               sim.type=sim.type,
                               verbose=FALSE)
  expect_equal(nrow(sim.data$train) + nrow(sim.data$holdout) + nrow(sim.data$validation),
               num.samples * 3)
  expect_equal(ncol(sim.data$train), num.variables + 1)
})

test_that("createSimulation with interaction network returns sane results", {
  sim.type <- "interactionScalefree"
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  bias <- 0.4
  sim.data <- createSimulation(num.samples=num.samples,
                               num.variables=num.variables,
                               pct.signals=pct.signals,
                               bias=bias,
                               sim.type=sim.type,
                               verbose=FALSE)
  expect_equal(nrow(sim.data$train) + nrow(sim.data$holdout) + nrow(sim.data$validation),
               num.samples * 3)
  expect_equal(ncol(sim.data$train), num.variables + 1)
})

test_that("splitDataset makes the right sized splits", {
  data("fullfMRI2")
  n <- nrow(fullfMRI2)
  data.sets <- splitDataset(all.data=fullfMRI2,
                            pct.train=0.5,
                            pct.holdo=0.5,
                            pct.validation=0,
                            class.label="phenos")
  # splits
  expect_equal(nrow(data.sets$train), n * 0.5, tolerance=0.5)
  expect_equal(nrow(data.sets$holdout), n * 0.5, tolerance=0.5)
  # total
  if(is.null(nrow(data.sets$validation))) {
    expect_equal(nrow(data.sets$train) + nrow(data.sets$holdout), nrow(fullfMRI2))
  } else {
    expect_equal(nrow(data.sets$train) + nrow(data.sets$holdout) + nrow(data.sets$validation),
                 nrow(fullfMRI2))
  }
})
