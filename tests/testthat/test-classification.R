# test-classification.R - Bill White - 4/26/17
#
# Test the privateEC classification functions on simulated data.

library(privateEC)
context("Classification")

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

test_that("originalPrivacy returns sane results", {
  num.samples <- 100
  num.variables <- 100
  pb <- 0.1
  nbias <- pb * num.variables
  signals <- sprintf("gene%04d", 1:nbias)
  temp.pec.file <- tempfile(pattern="pEc_temp", tmpdir=tempdir())

  data.sets <- createSimulation(d=num.variables,
                                n=num.samples,
                                type="sva",
                                verbose=FALSE)
  pec.results <- privateEC(data.sets=data.sets,
                           is.simulated=TRUE,
                           n=num.samples,
                           d=num.variables,
                           signal.names=signals,
                           save.file=temp.pec.file,
                           verbose=FALSE)
  por.results <- originalPrivacy(data.sets=data.sets,
                                 is.simulated=TRUE,
                                 n=num.samples,
                                 d=num.variables,
                                 signal.names=signals,
                                 pec.file=temp.pec.file,
                                 verbose=FALSE)
  file.remove(temp.pec.file)

  expect_equal(ncol(pec.results$plots.data), 5)
  expect_equal(ncol(pec.results$melted.data), 4)
  expect_equal(length(pec.results$correct), nrow(pec.results$plots.data))

  expect_equal(ncol(por.results$plots.data), 5)
  expect_equal(ncol(por.results$melted.data), 4)
  expect_equal(length(por.results$correct), nrow(por.results$plots.data))
})

test_that("privateRF returns sane results", {
  num.samples <- 100
  num.variables <- 100
  pb <- 0.1
  nbias <- pb * num.variables
  signals <- sprintf("gene%04d", 1:nbias)
  temp.pec.file <- tempfile(pattern="pEc_temp", tmpdir=tempdir())

  data.sets <- createSimulation(d=num.variables,
                                n=num.samples,
                                type="sva",
                                verbose=FALSE)
  pec.results <- privateEC(data.sets=data.sets,
                           is.simulated=TRUE,
                           n=num.samples,
                           d=num.variables,
                           signal.names=signals,
                           save.file=temp.pec.file,
                           verbose=FALSE)
  prf.results <- privateRF(data.sets=data.sets,
                           is.simulated=TRUE,
                           n=num.samples,
                           d=num.variables,
                           signal.names=signals,
                           pec.file=temp.pec.file,
                           verbose=FALSE)
  file.remove(temp.pec.file)

  expect_equal(ncol(pec.results$plots.data), 5)
  expect_equal(ncol(pec.results$melted.data), 4)
  expect_equal(length(pec.results$correct), nrow(pec.results$plots.data))

  expect_equal(ncol(prf.results$plots.data), 5)
  expect_equal(ncol(prf.results$melted.data), 4)
  expect_equal(length(prf.results$correct), nrow(prf.results$plots.data))
})

test_that("standard random forest returns sane results", {
  num.samples <- 100
  num.variables <- 100
  pb <- 0.1
  nbias <- pb * num.variables
  signals <- sprintf("gene%04d", 1:nbias)
  data.sets <- createSimulation(d=num.variables, n=num.samples,
                                type="sva", verbose=FALSE)
  rra.results <- standardRF(data.sets=data.sets,
                            is.simulated=TRUE,
                            type=type,
                            verbose=FALSE,
                            signal.names=signals)

  expect_equal(ncol(rra.results$plots.data), 5)
  expect_equal(ncol(rra.results$melted.data), 4)
  expect_equal(length(rra.results$correct), nrow(rra.results$plots.data))
})

test_that("privateECinbix returns sane results", {
  num.samples <- 100
  num.variables <- 100
  pb <- 0.1
  nbias <- pb * num.variables
  signals <- sprintf("gene%04d", 1:nbias)
  sim.data <- createSimulation(d=num.variables, n=num.samples,
                               type="sva", verbose=FALSE)
  pec.results <- privateECinbix(data.sets=sim.data,
                                n=num.samples,
                                verbose=FALSE)
  expect_equal(ncol(pec.results$plots.data), 5)
  expect_equal(ncol(pec.results$melted.data), 4)
  expect_equal(length(pec.results$correct), nrow(pec.results$plots.data))
})
