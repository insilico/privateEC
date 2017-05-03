#' Workflow for running one simulation of the Bioinformatics paper workflow
#'
#' Creates one simulation of train/holdout/test data sets, then runs the
#' four algorithms on that data. Returns a data frame of run results for each.
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
#' @note Default parameter values match those from the Bioinformatics paper.
#' @export
paperSimWorkflow <- function(myrun="001",
                             n=100,
                             d=100,
                             pb=0.1,
                             update.freq=50,
                             verbose=FALSE) {
  ptm <- proc.time()
  if(verbose) cat("run ID:", myrun, "\n")
  types <- c("sva", "er", "inte")
  biases <- c(0.4, 0.4, 0.4)
  alg.steps <- seq(d + 1, 1, -update.freq)[-1]
  num.steps <- length(alg.steps)
  temp.pec.file <- tempfile(pattern="pEc_temp", tmpdir=tempdir())
  nbias <- pb * d
  num.sigs <- pb * d
  signal.names <- sprintf("gene%04d", 1:nbias)

  all.run.results <- lapply(1:length(types), FUN=function(simtype.num) {
    type <- types[[simtype.num]]
    bias <- biases[[simtype.num]]
    if(verbose) cat("begin type/sim/classification loop for type/bias",
                    type, bias, "\n")
    if(verbose) cat("running simulation with n d pb", n, d, pb, "\n")
    shortname <- paste(round(bias, digits=1), pb, d, n, myrun, sep="_")
    if(verbose) cat("--------------------------------------------\n")
    data.sets <- createSimulation(n=n,
                                  d=d,
                                  pb=pb,
                                  bias=bias,
                                  shortname=shortname,
                                  type=type,
                                  myrun=myrun,
                                  verbose=verbose,
                                  save.file=FALSE)
    if(verbose) cat("running private algorithms\n")
    if(verbose) cat("\nrunning privateEC\n")
    pec.result <- privateEC(data.sets=data.sets,
                            is.simulated=TRUE,
                            n=n,
                            shortname=shortname,
                            bias=bias,
                            type=type,
                            myrun=myrun,
                            update.freq=update.freq,
                            save.file=temp.pec.file,
                            verbose=verbose,
                            signal.names=signal.names)
    if(verbose) cat("\nrunning originalprivacy.R\n")
    por.result <- originalPrivacy(data.sets=data.sets,
                                  is.simulated=TRUE,
                                  n=n,
                                  d=d,
                                  shortname=shortname,
                                  type=type,
                                  myrun=myrun,
                                  verbose=verbose,
                                  signal.names=signal.names,
                                  pec.file=temp.pec.file)
    if(verbose) cat("\nrunning privaterf.R\n")
    pra.result <- privateRF(data.sets=data.sets,
                            n=n,
                            d=d,
                            shortname=shortname,
                            type=type,
                            verbose=verbose,
                            signal.names=signal.names,
                            pec.file=temp.pec.file)
    if(verbose) cat("\nrunning regularRF\n")
    rra.result <- standardRF(data.sets=data.sets,
                             is.simulated=TRUE,
                             shortname=shortname,
                             type=type,
                             verbose=verbose,
                             signal.names=signal.names)

    if(!is.null(temp.pec.file)) {
      file.remove(temp.pec.file)
    }

    # compile and return results
    all.results = list(pec=pec.result,
                       por=por.result,
                       pra=pra.result,
                       rra=rra.result)
    run.results <- compileResults(run.results=all.results, verbose=verbose)
    if(verbose) cat("end sim/classification loop for type/bias",
                    type, bias, "\n")
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

#' Runs the four comparison algorithms on passed data set.
#' Returns a data frame of run results for each.
#'
#' @param corr.mat A matrix subject by region-pair correlations
#' @param phenos A vector of factors representing phenotypes for subjects
#' @param update.freq A integer for the number of steps before update
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list containing:
#' \describe{
#'   \item{run.results}{data frame of run results of each sim type}
#'   \item{elapsed}{total elapsed time}
#' }
#' @export
paperRealWorkflow <- function(corr.mat=NULL,
                              phenos=NULL,
                              update.freq=50,
                              verbose=FALSE) {
  if(is.null(corr.mat)) {
    stop("paperRealWorkflow: No subject by correlation matrix provided")
  }
  if(is.null(phenos)) {
    stop("paperRealWorkflow: No phenos provided")
  }
  ptm <- proc.time()

  temp.pec.file <- tempfile(pattern="pEc_temp", tmpdir=tempdir())

  # transform the data to fit the workflow expectations
  data <- data.frame(corr.mat, pheno=phenos)
  n <- nrow(data)
  d <- ncol(data)
  ind <- sample(2, n, replace=T)
  ind.case <- sample(2, n, replace=T)[1:floor(n / 2)]
  ind.ctrl <- sample(ind.case, floor(n / 2))
  ind.case <- c(ind.case, 2)
  data <- data.frame(data)
  data[,d] <- factor(data[, d])
  levels(data[,d]) <- c(-1, 1)
  data.case <- data[data[d]==1, ]
  data.ctrl <- data[data[d]==-1, ]
  X_train <- rbind(data.case[ind.case == 1, ], data.ctrl[ind.ctrl == 1, ])
  X_holdo <- rbind(data.case[ind.case == 2, ], data.ctrl[ind.ctrl == 2, ])
  data.sets <- list(train=X_train, holdout=X_holdo)

  if(verbose) cat("\nrunning privateEC\n")
  pec.result <- privateEC(data.sets=data.sets,
                          is.simulated=FALSE,
                          n=n,
                          d=d,
                          shortname="fmri",
                          type="REAL",
                          myrun="000",
                          update.freq=update.freq,
                          save.file=temp.pec.file,
                          verbose=verbose)
  if(verbose) cat("\nrunning originalprivacy.R\n")
  por.result <- originalPrivacy(data.sets=data.sets,
                                is.simulated=FALSE,
                                n=n,
                                d=d,
                                shortname="fmri",
                                type="REAL",
                                myrun="000",
                                verbose=verbose,
                                pec.file=temp.pec.file)
  if(verbose) cat("\nrunning privaterf.R\n")
  pra.result <- privateRF(data.sets=data.sets,
                          n=n,
                          d=d,
                          shortname="fmri",
                          type="REAL",
                          is.simulated=FALSE,
                          pec.file=temp.pec.file,
                          verbose=verbose)
  if(verbose) cat("\nrunning regularRF\n")
  rra.result <- standardRF(data.sets=data.sets,
                           is.simulated=FALSE,
                           shortname="fmri",
                           type="REAL",
                           verbose=verbose)

  file.remove(temp.pec.file)

  # compile and return results
  all.results = list(pec=pec.result,
                     por=por.result,
                     pra=pra.result,
                     rra=rra.result)
  run.results <- compileResults(run.results=all.results,
                                verbose=verbose)

  elapsed <- (proc.time() - ptm)[3]
  if(verbose) cat("Total elapsed time:", elapsed, "\n")

  list(run.results=run.results, elapsed=elapsed)
}
