#' Workflow for running one simulation of the Bioinformatics paper workflow
#'
#' @param myrun A character vector identifying the run
#' @param n An integer for the number of samples
#' @param d An integer for the number of variables
#' @param pb A numeric for the significant variable bias
#' @param update.freq A integer for the number of steps before update
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list containing:
#' \describe{
#'   \item{run.results}{data frame of run results of each sim type}
#'   \item{elapsed}{total elapsed time}
#' }
#' @export
privateECworkflow <- function(myrun="001",
                              n=100,
                              d=100,
                              pb=0.1,
                              update.freq=50,
                              verbose=FALSE) {
  ptm <- proc.time()
  cat("run ID:", myrun, "\n")
  types <- c("sva", "er", "inte")
  biases <- c(0.4, 0.4, 0.4)
  alg.steps <- seq(d + 1, 1, -update.freq)[-1]
  num.steps <- length(alg.steps)
  save.file <- "tmp"
  nbias <- pb * d
  num.sigs <- pb * d
  signal.names <- sprintf("gene%04d", 1:nbias)

  all.run.results <- lapply(1:length(types), FUN=function(simtype.num) {
    type <- types[[simtype.num]]
    bias <- biases[[simtype.num]]
    if(verbose) cat("begin type/sim/classification loop for type/bias", type,
                    bias, "\n")
    if(verbose) cat("running simulation with n d pb", n, d, pb, "\n")
    shortname <- paste(round(bias, digits=1), pb, d, n, myrun, sep="_")
    if(verbose) cat("--------------------------------------------\n")
    data.sets <- createSimulation(n, d, pb, bias, shortname, type, myrun,
                                  verbose, save.file=FALSE)
    if(verbose) cat("running private algorithms\n")
    accuracy.ls <- list()
    correct.vars.ls <- list()
    if(verbose) cat("\nrunning privateEC\n")
    pec.result <- privateEC(data.sets, n, shortname, bias, type,
                            myrun, update.freq, save.file=TRUE, verbose)
    if(verbose) cat("\nrunning originalprivacy.R\n")
    por.result <- originalPrivacy(data.sets, n, d, shortname, type, myrun,
                                  save.file=FALSE, verbose)
    if(verbose) cat("\nrunning privaterf.R\n")
    pra.result <- privateRF(data.sets, n, d, shortname, type,
                            save.file=FALSE, verbose)
    if(verbose) cat("\nrunning regularRF\n")
    rra.result <- standardRF(data.sets, shortname, type, save.file=FALSE,
                             verbose)

    # compile and return results
    file.remove(save.file)
    all.results = list(pec=pec.result,
                       por=por.result,
                       pra=pra.result,
                       rra=rra.result)
    run.results <- compileResults(all.results, save.file=FALSE, verbose=verbose)
    if(verbose) cat("end sim/classification loop for type/bias", type, bias, "\n")
    run.results
  }) # end simtype.num loop 1:3 for sva. er, inte

  elapsed <- (proc.time() - ptm)[3]
  cat("Total elapsed time:", elapsed, "\n")

  list(run.results=all.run.results, elapsed=elapsed)
}

#' Compile the results of a simulation + classifier methods run
#'
#' @param run.results A list of run results
#' @param save.file A character vector for results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list with:
#' \describe{
#'   \item{plots.data}{data frame of results, a row for each update}
#'   \item{melted.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#' }
#' @export
compileResults <- function(run.results=NULL,
                           save.file=NULL,
                           verbose=FALSE) {
  if(is.null(run.results)) {
    stop("compileResults: No results list provided as first argument")
  }
  if(verbose) cat("compiling results for plotting\n")
  plots.dfs <- lapply(run.results, FUN=function(method.results) {
    method.results$plots.data
  })
  fp.plots <- do.call(rbind, plots.dfs)
  melted.dfs <- lapply(run.results, FUN=function(method.results) {
    method.results$melted.data
  })
  fp.melted <- do.call(rbind, melted.dfs)
  correct.dfs <- lapply(run.results, FUN=function(method.results) {
    method.results$correct
  })
  # ---------------------------------------------------------------------------
  if(!is.null(save.file)) {
    if(verbose) cat("saving compiled results", save.file, "\n")
    save(fp.plots, plots.dfs, fp.melted, file=save.file)
  }
  list(plots.data=fp.plots,
       melted.data=fp.melted,
       correct=correct.dfs)
}
