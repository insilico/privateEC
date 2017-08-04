# workflow.R - Bill White - May 2017
#
# Worflow algorithms replicating those in the Bioinformatics paper:
# Differential privacy-based Evaporative Cooling feature selection and
# classification with Relief-F and Random Forests
# https://doi.org/10.1093/bioinformatics/btx298

#' Workflow for running one simulation of the Bioinformatics paper workflow
#'
#' Creates one simulation of train/holdout/validation data sets, then runs the
#' four algorithms on that data. Returns a data frame of run results for each.
#'
#' @param n.samples An integer for the number of samples
#' @param n.variables An integer for the number of variables
#' @param pct.signals A numeric for the significant variable bias
#' @param update.freq A integer for the number of steps before update
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list containing:
#' \describe{
#'   \item{run.results}{data frame of run results of each sim type}
#'   \item{elapsed}{total elapsed time}
#' }
#' @note Default parameter values match those from the Bioinformatics paper.
#' @family workflows
#' @examples
#'   num.samples <- 100
#'   num.variables <- 100
#'   pct.signals <- 0.1
#'   upd.frq <- 0.1 * num.variables
#'   one.step.result <- paperSimWorkflow(n.samples=num.samples,
#'                                       n.variables=num.variables,
#'                                       pct.signals=pct.signals,
#'                                       update.freq=upd.frq,
#'                                       verbose=FALSE)
#' @seealso The workflow consists of the sequence:
#' \code{\link{createSimulation}}
#' \code{\link{privateEC}} (optionally \code{\link{privateECinbix}})
#' \code{\link{originalThresholdout}}
#' \code{\link{privateRF}}
#' \code{\link{standardRF}} and
#' \code{\link{compileResults}}. A comparison analysis with real data (fMRI)
#' is in \code{\link{paperRealWorkflow}}.
#' @export
paperSimWorkflow <- function(n.samples=100,
                             n.variables=100,
                             pct.signals=0.1,
                             update.freq=10,
                             verbose=FALSE) {
  ptm <- proc.time()
  types <- c("mainEffect", "interactionErdos", "interactionScalefree")
  biases <- c(0.4, 0.4, 0.4)
  alg.steps <- seq(n.variables + 1, 1, -update.freq)[-1]
  num.steps <- length(alg.steps)
  temp.pec.file <- tempfile(pattern="pEc_temp", tmpdir=tempdir())
  num.sigs <- pct.signals * n.variables
  signal.names <- sprintf("gene%04d", 1:num.sigs)

  all.run.results <- lapply(seq_along(types), FUN=function(simtype.num) {
    type <- types[[simtype.num]]
    bias <- biases[[simtype.num]]
    if(verbose) cat("begin type/sim/classification loop for type/bias",
                    simtype.num, bias, "\n")
    if(verbose) cat("running simulation with n num.vars pct.signals",
                    n, n.variables, pct.signals, "\n")
    if(verbose) cat("--------------------------------------------\n")
    data.sets <- createSimulation(num.samples=n.samples,
                                  num.variables=n.variables,
                                  pct.signals=pct.signals,
                                  bias=bias,
                                  sim.type=type,
                                  verbose=verbose)
    if(verbose) cat("running private algorithms\n")
    if(verbose) cat("\nrunning privateEC\n")
    pec.result <- privateEC(train.ds=data.sets$train,
                            holdout.ds=data.sets$holdout,
                            validation.ds=data.sets$validation,
                            label=data.sets$class.label,
                            is.simulated=TRUE,
                            update.freq=update.freq,
                            save.file=temp.pec.file,
                            signal.names=data.sets$signal.names,
                            verbose=verbose)
    if(verbose) cat("\nrunning originalThresholdout.R\n")
    por.result <- originalThresholdout(train.ds=data.sets$train,
                                       holdout.ds=data.sets$holdout,
                                       validation.ds=data.sets$validation,
                                       label=data.sets$class.label,
                                       is.simulated=TRUE,
                                       signal.names=signal.names,
                                       pec.file=temp.pec.file,
                                       verbose=verbose)
    if(verbose) cat("\nrunning privaterf.R\n")
    pra.result <- privateRF(train.ds=data.sets$train,
                            holdout.ds=data.sets$holdout,
                            validation.ds=data.sets$validation,
                            label=data.sets$class.label,
                            is.simulated=TRUE,
                            signal.names=signal.names,
                            pec.file=temp.pec.file,
                            verbose=verbose)
    if(verbose) cat("\nrunning regularRF\n")
    rra.result <- standardRF(train.ds=data.sets$train,
                             holdout.ds=data.sets$holdout,
                             validation.ds=data.sets$validation,
                             label=data.sets$class.label,
                             is.simulated=TRUE,
                             signal.names=signal.names,
                             verbose=verbose)

    if(!is.null(temp.pec.file)) {
      file.remove(temp.pec.file)
    }

    # compile and return results
    all.results <- list(pec=pec.result,
                        por=por.result,
                        pra=pra.result,
                        rra=rra.result)
    run.results <- compileResults(run.results=all.results, verbose=verbose)
    if(verbose) cat("end sim/classification loop for type/bias", type, bias, "\n")
    run.results
  }) # end simtype.num loop 1:3 for sva. er, inte

  elapsed <- (proc.time() - ptm)[3]
  if(verbose) cat("Total elapsed time:", elapsed, "\n")

  list(run.results=all.run.results, elapsed=elapsed)
}

#' Compile the results of a simulation + classifier methods run
#'
#' @param run.results A list of run results
#' @param save.file A character vector for results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list with:
#' \describe{
#'   \item{algo.acc}{data frame of results, a row for each update}
#'   \item{ggplot.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#' }
#' @family workflows
compileResults <- function(run.results=NULL,
                           save.file=NULL,
                           verbose=FALSE) {
  if(is.null(run.results)) {
    stop("compileResults: No results list provided as first argument")
  }
  if(verbose) cat("compiling results for plotting\n")
  acc.dfs <- lapply(run.results, FUN=function(method.results) {
    method.results$algo.acc
  })
  all.acc <- do.call(rbind, acc.dfs)
  ggplot.dfs <- lapply(run.results, FUN=function(method.results) {
    method.results$ggplot.data
  })
  all.ggplot <- do.call(rbind, ggplot.dfs)
  all.correct <- lapply(run.results, FUN=function(method.results) {
    method.results$correct
  })
  if(!is.null(save.file)) {
    if(verbose) cat("saving compiled results", save.file, "\n")
    save(all.acc, all.correct, all.ggplot, file=save.file)
  }

  list(algo.acc=all.acc,
       ggplot.data=all.ggplot,
       correct=all.correct)
}

#' Runs the four comparison algorithms on the passed correlation matrix and phenotypes.
#' Returns a data frame of run results for each.
#'
#' @param real.data A matrix subject by region-pair correlations augmented with class label
#' @param label A character vector class label from columns of corr,mat
#' @param update.freq A integer for the number of steps before update
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list containing:
#' \describe{
#'   \item{run.results}{data frame of run results of each sim type}
#'   \item{elapsed}{total elapsed time}
#' }
#' @examples
#'   data(fullfMRI2)
#'   real.result <- paperRealWorkflow(real.data=fullfMRI2[, 2900:ncol(fullfMRI2)],
#'                                    label="phenos",
#'                                    update.freq=50,
#'                                    verbose=FALSE)
#' @seealso The workflow consists of the sequence:
#' \code{\link{privateEC}} (optionally \code{\link{privateECinbix}})
#' \code{\link{originalThresholdout}}
#' \code{\link{privateRF}}
#' \code{\link{standardRF}} and
#' \code{\link{compileResults}}. A comparison analysis with simulated data
#' is in \code{\link{paperSimWorkflow}}.
#' @family workflows
#' @export
paperRealWorkflow <- function(real.data=NULL,
                              label=NULL,
                              update.freq=50,
                              verbose=FALSE) {
  if(is.null(real.data)) {
    stop("No subject by correlation matrix provided")
  }
  ptm <- proc.time()

  temp.pec.file <- tempfile(pattern="pEc_temp", tmpdir=tempdir())
  # transform the data to fit the workflow expectations
  n <- nrow(real.data)
  num.vars <- ncol(real.data) - 1
  real.data.sets <- splitDataset(all.data=real.data,
                                 pct.train=0.5,
                                 pct.holdout=0.5,
                                 pct.validation=0,
                                 class.label="phenos")
  if(verbose) cat("\nrunning privateEC\n")
  pec.result <- privateEC(train.ds=real.data.sets$train,
                          holdout.ds=real.data.sets$holdout,
                          validation.ds=NULL,
                          label=label,
                          is.simulated=FALSE,
                          update.freq=update.freq,
                          save.file=temp.pec.file,
                          verbose=verbose)
  if(verbose) cat("\nrunning originalThresholdout.R\n")
  por.result <- originalThresholdout(train.ds=real.data.sets$train,
                                     holdout.ds=real.data.sets$holdout,
                                     validation.ds=NULL,
                                     label=label,
                                     is.simulated=FALSE,
                                     pec.file=temp.pec.file,
                                     verbose=verbose)
  if(verbose) cat("\nrunning privaterf.R\n")
  pra.result <- privateRF(train.ds=real.data.sets$train,
                          holdout.ds=real.data.sets$holdout,
                          validation.ds=NULL,
                          label=label,
                          is.simulated=FALSE,
                          pec.file=temp.pec.file,
                          verbose=verbose)
  if(verbose) cat("\nrunning regularRF\n")
  rra.result <- standardRF(train.ds=real.data.sets$train,
                           holdout.ds=real.data.sets$holdout,
                           validation.ds=NULL,
                           label=label,
                           is.simulated=FALSE,
                           verbose=verbose)

  file.remove(temp.pec.file)

  # compile and return results
  all.results <- list(pec=pec.result, por=por.result, pra=pra.result, rra=rra.result)
  final.results <- compileResults(run.results=all.results, verbose=verbose)

  elapsed <- (proc.time() - ptm)[3]
  if(verbose) cat("Total elapsed time:", elapsed, "\n")

  list(run.results=final.results, elapsed=elapsed)
}
