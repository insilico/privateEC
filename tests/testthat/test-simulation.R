# test-simulation.R - Bill White - 4/30/17
# 
# modified by Saeid Parvandeh - July 2018
#
# Test the creation of simulated data.

library(privateEC)
context("Simulation")

test_that("createSimulation mainEffect returns sane results - with dichotomous outcome", {
  sim.type <- "mainEffect"
  label <- "class"
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  bias <- 0.4
  sim.data <- createSimulation(num.samples = num.samples,
                               num.variables = num.variables,
                               label = label,
                               pct.signals = pct.signals,
                               bias = bias,
                               sim.type = sim.type,
                               verbose = FALSE)
  expect_equal(nrow(sim.data$train) + nrow(sim.data$holdout) + nrow(sim.data$validation),
               num.samples)
  expect_equal(ncol(sim.data$train), num.variables + 1)
})

test_that("createSimulation Erdos-Renyi network base returns sane results - with dichotomous outcome", {
  sim.type <- "interactionErdos"
  label <- "class"
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  bias <- 0.4
  sim.data <- createSimulation(num.variables = num.variables,
                               num.samples = num.samples,
                               label = label,
                               pct.signals = pct.signals,
                               bias = bias,
                               sim.type = sim.type,
                               verbose = FALSE)
  expect_equal(nrow(sim.data$train) + nrow(sim.data$holdout) + nrow(sim.data$validation),
               num.samples)
  expect_equal(ncol(sim.data$train), num.variables + 1)
})

test_that("createSimulation with interaction network returns sane results - with dichotomous outcome", {
  sim.type <- "interactionScalefree"
  label <- "class"
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  bias <- 0.4
  sim.data <- createSimulation(num.samples = num.samples,
                               num.variables = num.variables,
                               label = label,
                               pct.signals = pct.signals,
                               bias = bias,
                               sim.type = sim.type,
                               verbose = FALSE)
  expect_equal(nrow(sim.data$train) + nrow(sim.data$holdout) + nrow(sim.data$validation),
               num.samples)
  expect_equal(ncol(sim.data$train), num.variables + 1)
})

test_that("createMixedSimulation with interaction network returns sane results - with dichotomous outcome", {
  mixed.type <- c("mainEffect","interactionScalefree")
  label <- "class"
  num.samples <- 100
  num.variables <- 100
  pct.mixed <- .5
  pct.signals <- 0.1
  bias <- 0.4
  sim.data <- createMixedSimulation(num.samples = num.samples,
                                    num.variables = num.variables,
                                    label = label,
                                    pct.signals = pct.signals,
                                    bias = bias,
                                    pct.mixed=pct.mixed,
                                    mixed.type=mixed.type,
                                    verbose = FALSE)
  expect_equal(nrow(sim.data$train) + nrow(sim.data$holdout) + nrow(sim.data$validation),
               num.samples)
  expect_equal(ncol(sim.data$train), num.variables + 1)
})

test_that("createSimulation mainEffect returns sane results - with quantitative outcome", {
  sim.type <- "mainEffect"
  label <- "qtrait"
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  bias <- 0.4
  sim.data <- createSimulation(num.samples = num.samples,
                               num.variables = num.variables,
                               label = label,
                               pct.signals = pct.signals,
                               bias = bias,
                               sim.type = sim.type,
                               verbose = FALSE)
  expect_equal(nrow(sim.data$train) + nrow(sim.data$holdout) + nrow(sim.data$validation),
               num.samples)
  expect_equal(ncol(sim.data$train), num.variables + 1)
})

test_that("createSimulation Erdos-Renyi network base returns sane results - with quantitative outcome", {
  sim.type <- "interactionErdos"
  label <- "qtrait"
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  bias <- 0.4
  sim.data <- createSimulation(num.variables = num.variables,
                               num.samples = num.samples,
                               label = label,
                               pct.signals = pct.signals,
                               bias = bias,
                               sim.type = sim.type,
                               verbose = FALSE)
  expect_equal(nrow(sim.data$train) + nrow(sim.data$holdout) + nrow(sim.data$validation),
               num.samples)
  expect_equal(ncol(sim.data$train), num.variables + 1)
})

test_that("createSimulation with interaction network returns sane results - with quantitative outcome", {
  sim.type <- "interactionScalefree"
  label <- "qtrait"
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  bias <- 0.4
  sim.data <- createSimulation(num.samples = num.samples,
                               num.variables = num.variables,
                               label = label,
                               pct.signals = pct.signals,
                               bias = bias,
                               sim.type = sim.type,
                               verbose = FALSE)
  expect_equal(nrow(sim.data$train) + nrow(sim.data$holdout) + nrow(sim.data$validation),
               num.samples)
  expect_equal(ncol(sim.data$train), num.variables + 1)
})

test_that("createMixedSimulation with interaction network returns sane results - with quantitative outcome", {
  mixed.type <- c("mainEffect","interactionScalefree")
  label <- "qtrait"
  num.samples <- 100
  num.variables <- 100
  pct.mixed <- .5
  pct.signals <- 0.1
  bias <- 0.4
  sim.data <- createMixedSimulation(num.samples = num.samples,
                                    num.variables = num.variables,
                                    label = label,
                                    pct.signals = pct.signals,
                                    bias = bias,
                                    pct.mixed=pct.mixed,
                                    mixed.type=mixed.type,
                                    verbose = FALSE)
  expect_equal(nrow(sim.data$train) + nrow(sim.data$holdout) + nrow(sim.data$validation),
               num.samples)
  expect_equal(ncol(sim.data$train), num.variables + 1)
})

test_that("splitDataset makes the correct sized train, holdout sized splits", {
  data("rsfMRIcorrMDD")
  n <- nrow(rsfMRIcorrMDD)
  data.sets <- splitDataset(all.data = rsfMRIcorrMDD,
                            pct.train = 0.5,
                            pct.holdout = 0.5,
                            pct.validation = 0,
                            label = "class")
  # splits
  expect_equal(nrow(data.sets$train), n * 0.5, tolerance = 0.5)
  expect_equal(nrow(data.sets$holdout), n * 0.5, tolerance = 0.5)
  expect_equal(nrow(data.sets$validation), 0)
  # total
  if (is.null(nrow(data.sets$validation))) {
    expect_equal(nrow(data.sets$train) + nrow(data.sets$holdout), nrow(rsfMRIcorrMDD))
  } else {
    expect_equal(nrow(data.sets$train) +
                   nrow(data.sets$holdout) +
                   nrow(data.sets$validation),
                 nrow(rsfMRIcorrMDD))
  }
})

test_that("splitDataset makes the correct sized train, holdout AND validation sized splits", {
  data("rsfMRIcorrMDD")
  n <- nrow(rsfMRIcorrMDD)
  data.sets <- splitDataset(all.data = rsfMRIcorrMDD,
                            pct.train = 1 / 3,
                            pct.holdout = 1 / 3,
                            pct.validation = 1 / 3,
                            label = "class")
  # splits
  a.third <- 1 / 3
  ntimes.a.third <- floor(n * a.third)
  expect_equal(nrow(data.sets$train), ntimes.a.third + 1)
  expect_equal(nrow(data.sets$holdout), ntimes.a.third)
  expect_equal(nrow(data.sets$validation), ntimes.a.third)
  # total
  if (is.null(nrow(data.sets$validation))) {
    expect_equal(nrow(data.sets$train) +
                   nrow(data.sets$holdout),
                 nrow(rsfMRIcorrMDD))
  } else {
    expect_equal(nrow(data.sets$train) +
                   nrow(data.sets$holdout) +
                   nrow(data.sets$validation),
                 nrow(rsfMRIcorrMDD))
  }
})
