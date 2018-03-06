# test-classification.R - Bill White - 4/26/17
#
# Test the privateEC classification functions on simulated data.

library(privateEC)
context("Classification")

test_that("privateEC returns sane results - Relief-F + randomForest", {
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  rf.mtry <- 5
  sim.data <- createSimulation(num.samples = num.samples,
                               num.variables = num.variables,
                               pct.signals = pct.signals,
                               pct.train = 1 / 3,
                               pct.holdout = 1 / 3,
                               pct.validation = 1 / 3,
                               sim.type = "mainEffect",
                               verbose = FALSE)
  pec.results <- privateEC(train.ds = sim.data$train,
                           holdout.ds = sim.data$holdout,
                           validation.ds = sim.data$validation,
                           label = sim.data$class.label,
                           importance.name = "relieff",
                           learner.name = "randomforest",
                           rf.mtry = rf.mtry,
                           is.simulated = TRUE,
                           signal.names = sim.data$signal.names,
                           verbose = FALSE)
  expect_equal(ncol(pec.results$algo.acc), 5)
  expect_equal(ncol(pec.results$ggplot.data), 4)
  expect_equal(length(pec.results$correct), nrow(pec.results$algo.acc))
})

test_that("privateEC returns sane results - Relief-F + xgboost fixed parameters", {
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  sim.data <- createSimulation(num.samples = num.samples,
                               num.variables = num.variables,
                               pct.signals = pct.signals,
                               pct.train = 1 / 3,
                               pct.holdout = 1 / 3,
                               pct.validation = 1 / 3,
                               sim.type = "mainEffect",
                               verbose = FALSE)
  # use of scalars reduces grid search to one dimension
  hyper.rounds <- 1
  hyper.depth <- 4
  hyper.shrink <- 1.0
  pec.results <- privateEC(train.ds = sim.data$train,
                           holdout.ds = sim.data$holdout,
                           validation.ds = sim.data$validation,
                           label = sim.data$class.label,
                           importance.name = "relieff",
                           importance.algorithm = "ReliefFBestK",
                           learner.name = "xgboost",
                           xgb.num.rounds = hyper.rounds,
                           xgb.max.depth = hyper.depth,
                           xgb.shrinkage = hyper.shrink,
                           is.simulated = TRUE,
                           signal.names = sim.data$signal.names,
                           verbose = FALSE)
  expect_equal(ncol(pec.results$algo.acc), 5)
  expect_equal(ncol(pec.results$ggplot.data), 4)
  expect_equal(length(pec.results$correct), nrow(pec.results$algo.acc))
})

test_that("privateEC returns sane results - Relief-F + xgboost with hyperparameter optimization", {
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  sim.data <- createSimulation(num.samples = num.samples,
                               num.variables = num.variables,
                               pct.signals = pct.signals,
                               pct.train = 1 / 3,
                               pct.holdout = 1 / 3,
                               pct.validation = 1 / 3,
                               sim.type = "mainEffect",
                               verbose = FALSE)
  # use vectors to perform grid search on the parameters
  hyper.rounds <- c(1, 2, 3)
  hyper.depth <- c(2, 4, 8)
  hyper.shrink <- c(0.1, 0.5, 1.0)
  pec.results <- privateEC(train.ds = sim.data$train,
                           holdout.ds = sim.data$holdout,
                           validation.ds = sim.data$validation,
                           label = sim.data$class.label,
                           importance.name = "relieff",
                           importance.algorithm = "ReliefFBestK",
                           learner.name = "xgboost",
                           xgb.num.rounds = hyper.rounds,
                           xgb.max.depth = hyper.depth,
                           xgb.shrinkage = hyper.shrink,
                           is.simulated = TRUE,
                           signal.names = sim.data$signal.names,
                           verbose = FALSE)
  expect_equal(ncol(pec.results$algo.acc), 5)
  expect_equal(ncol(pec.results$ggplot.data), 4)
  expect_equal(length(pec.results$correct), nrow(pec.results$algo.acc))
})

test_that("originalThresholdout returns sane results", {
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  temp.pec.file <- tempfile(pattern = "pEc_temp", tmpdir = tempdir())
  sim.data <- createSimulation(num.variables = num.variables,
                               num.samples = num.samples,
                               sim.type = "mainEffect",
                               pct.train = 1 / 3,
                               pct.holdout = 1 / 3,
                               pct.validation = 1 / 3,
                               verbose = FALSE)
  pec.results <- privateEC(train.ds = sim.data$train,
                           holdout.ds = sim.data$holdout,
                           validation.ds = sim.data$validation,
                           label = sim.data$class.label,
                           is.simulated = TRUE,
                           signal.names = sim.data$signal.names,
                           save.file = temp.pec.file,
                           verbose = FALSE)
  por.results <- originalThresholdout(train.ds = sim.data$train,
                                      holdout.ds = sim.data$holdout,
                                      validation.ds = sim.data$validation,
                                      label = sim.data$class.label,
                                      is.simulated = TRUE,
                                      signal.names = sim.data$signal.names,
                                      pec.file = temp.pec.file,
                                      verbose = FALSE)
  file.remove(temp.pec.file)

  expect_equal(ncol(pec.results$algo.acc), 5)
  expect_equal(ncol(pec.results$ggplot.data), 4)
  expect_equal(length(pec.results$correct), nrow(pec.results$algo.acc))

  expect_equal(ncol(por.results$algo.acc), 5)
  expect_equal(ncol(por.results$ggplot.data), 4)
  expect_equal(length(por.results$correct), nrow(por.results$algo.acc))
})

test_that("privateRF returns sane results", {
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  temp.pec.file <- tempfile(pattern = "pEc_temp", tmpdir = tempdir())
  sim.data <- createSimulation(num.variables = num.variables,
                               num.samples = num.samples,
                               sim.type = "mainEffect",
                               pct.train = 1 / 3,
                               pct.holdout = 1 / 3,
                               pct.validation = 1 / 3,
                               verbose = FALSE)
  pec.results <- privateEC(train.ds = sim.data$train,
                           holdout.ds = sim.data$holdout,
                           validation.ds = sim.data$validation,
                           label = sim.data$class.label,
                           is.simulated = TRUE,
                           signal.names = sim.data$signal.names,
                           save.file = temp.pec.file,
                           verbose = FALSE)
  prf.results <- privateRF(train.ds = sim.data$train,
                           holdout.ds = sim.data$holdout,
                           validation.ds = sim.data$validation,
                           label = sim.data$class.label,
                           is.simulated = TRUE,
                           signal.names = sim.data$signal.names,
                           pec.file = temp.pec.file,
                           verbose = FALSE)
  file.remove(temp.pec.file)

  expect_equal(ncol(pec.results$algo.acc), 5)
  expect_equal(ncol(pec.results$ggplot.data), 4)
  expect_equal(length(pec.results$correct), nrow(pec.results$algo.acc))

  expect_equal(ncol(prf.results$algo.acc), 5)
  expect_equal(ncol(prf.results$ggplot.data), 4)
  expect_equal(length(prf.results$correct), nrow(prf.results$algo.acc))
})

test_that("standard random forest returns sane results", {
  num.samples <- 100
  num.variables <- 100
  pct.signals <- 0.1
  sim.data <- createSimulation(num.variables = num.variables,
                               num.samples = num.samples,
                               sim.type = "mainEffect",
                               pct.train = 1 / 3,
                               pct.holdout = 1 / 3,
                               pct.validation = 1 / 3,
                               verbose = FALSE)
  rra.results <- standardRF(train.ds = sim.data$train,
                            holdout.ds = sim.data$holdout,
                            validation.ds = sim.data$validation,
                            label = sim.data$class.label,
                            is.simulated = TRUE,
                            signal.names = sim.data$signal.names,
                            verbose = FALSE)

  expect_equal(ncol(rra.results$algo.acc), 5)
  expect_equal(ncol(rra.results$ggplot.data), 4)
  expect_equal(length(rra.results$correct), nrow(rra.results$algo.acc))
})
