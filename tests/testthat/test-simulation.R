# test-simulation.R - Bill White - 4/30/17
#
# Test the creation of simulated data.

library(privateEC)
context("Simulation")

test_that("createSimulation SVA returns sane results", {
  sim.type <- "sva"
  num.samples <- 100
  num.variables <- 100
  pb <- 0.1
  bias <- 0.4
  nbias <- pb * num.variables
  signals <- sprintf("gene%04d", 1:nbias)
  sim.data <- createSimulation(d=num.variables, n=num.samples, pb=pb,
                               bias=bias, type=sim.type, verbose=FALSE)
  expect_equal(nrow(sim.data$train) +
                 nrow(sim.data$holdout) +
                 nrow(sim.data$test), num.samples * 3)
  expect_equal(ncol(sim.data$train), num.variables)
})

test_that("createSimulation Erdos-Renyi network base returns sane results", {
  sim.type <- "er"
  num.samples <- 100
  num.variables <- 100
  pb <- 0.1
  bias <- 0.4
  nbias <- pb * num.variables
  signals <- sprintf("gene%04d", 1:nbias)
  sim.data <- createSimulation(d=num.variables, n=num.samples, pb=pb,
                               bias=bias, type=sim.type, verbose=FALSE)
  expect_equal(nrow(sim.data$train) +
                 nrow(sim.data$holdout) +
                 nrow(sim.data$test), num.samples * 3)
  expect_equal(ncol(sim.data$train), num.variables)
})

test_that("createSimulation with interaction network returns sane results", {
  sim.type <- "inte"
  num.samples <- 100
  num.variables <- 100
  pb <- 0.1
  bias <- 0.4
  nbias <- pb * num.variables
  signals <- sprintf("gene%04d", 1:nbias)
  sim.data <- createSimulation(d=num.variables, n=num.samples, pb=pb,
                               bias=bias, type=sim.type, verbose=FALSE)
  expect_equal(nrow(sim.data$train) +
                 nrow(sim.data$holdout) +
                 nrow(sim.data$test), num.samples * 3)
  expect_equal(ncol(sim.data$train), num.variables)
})
