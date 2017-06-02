# classification.R - Trang Le and Bill White - Fall 2016/Spring 2017
#
# Classification algorithms used in the Bioinformatics paper:
# Differential privacy-based Evaporative Cooling feature selection and
# classification with Relief-F and Random Forests
# https://doi.org/10.1093/bioinformatics/btx298

#' Compute and return importance scores (Relief-F scores)
#'
#' @param train.set A training data frame with last column as class
#' @param holdout.set A holdout data frame with last column as class
#' @param label A character vector of the class variable column name
#' @param imp.estimator A character vector CORElearn attribute importance estimator
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list with two data frames representing the importance scores
#' (Relief-F scores) for the train and holdout data sets.
getImportanceScores <- function(train.set=NULL,
                                holdout.set=NULL,
                                label="phenos",
                                imp.estimator="ReliefFbestK",
                                verbose=FALSE) {
  if(is.null(train.set)) {
    stop("getImportanceScores: No training data set provided")
  }
  if(is.null(holdout.set)) {
    stop("getImportanceScores: No holdout data set provided")
  }
  if(verbose) cat("\ttrain\n")
  train.relief <- CORElearn::attrEval(label,
                                      data=train.set,
                                      estimator=imp.estimator)
  if(verbose) cat("\tholdout\n")
  holdout.relief <- CORElearn::attrEval(label,
                                        data=holdout.set,
                                        estimator=imp.estimator)

  list(data.frame(train.relief), data.frame(holdout.relief))
}

#' Private Evaporative Cooling feature selection and classification
#'
#' @param train.ds A data frame with training data and class labels
#' @param holdout.ds A data frame with holdout data and class labels
#' @param validation.ds A data frame with validation data and class labels
#' @param label A character vector of the class variable column name
#' @param is.simulated Is the data simulated (or real?)
#' @param bias A numeric for effect size in simulated signal variables
#' @param update.freq An integer the number of steps before update
#' @param corelearn.estimator CORElearn Relief-F estimator
#' @param rf.ntree An integer the number of trees in the random forest
#' @param rf.mtry An integer the number of variables sampled at each random forest node split
#' @param start.temp A numeric EC starting temperature
#' @param final.temp A numeric EC final temperature
#' @param tau.param A numeric tau to control temperature reduction schedule
#' @param threshold A numeric, default 4 / sqrt(n) suggested in the
#'  thresholdout’s supplementary material (Dwork, et al.,2015)
#' @param tolerance A numeric, default 1 / sqrt(n) suggested in the
#'  thresholdout’s supplementary material (Dwork, et al.,2015)
#' @param signal.names A character vector of signal names in simulated data
#' @param save.file A character vector for results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list with:
#' \describe{
#'   \item{algo.acc}{data frame of results, a row for each update}
#'   \item{ggplot.data}{melted results data frame for plotting with ggplot}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @examples
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' sim.data <- createSimulation(n=num.samples,
#'                              num.vars=num.variables,
#'                              pct.signals=pct.signals,
#'                              sim.type="mainEffect",
#'                              verbose=FALSE)
#' pec.results <- privateEC(train.ds=sim.data$train,
#'                          holdout.ds=sim.data$holdout,
#'                          validation.ds=sim.data$validation,
#'                          label=sim.data$class.label,
#'                          is.simulated=TRUE,
#'                          signal.names=sim.data$signal.names,
#'                          verbose=FALSE)
#' @note Within thresholdout, we choose a threshold of 4 / sqrt(n) and
#' tolerance of 1 / sqrt(n) as suggested in the thresholdout’s supplementary
#' material (Dwork, et al., 2015).
#' @references
#' Trang Le, W. K. Simmons, M. Misaki, B.C. White, J. Savitz, J. Bodurka,
#' and B. A. McKinney. “Differential privacy-based Evaporative Cooling feature selection
#' and classification with Relief-F and Random Forests,”
#' Bioinformatics. Accepted. https://doi.org/10.1093/bioinformatics/btx298. 2017
#' @family classification
#' @export
privateEC <- function(train.ds=NULL, holdout.ds=NULL, validation.ds=NULL,
                      label="phenos",
                      is.simulated=TRUE,
                      bias=0.4,
                      update.freq=50,
                      corelearn.estimator="ReliefFbestK",
                      rf.ntree=500,
                      rf.mtry=NULL,
                      start.temp=0.1,
                      final.temp=10 ^ (-5),
                      tau.param=100,
                      threshold=4 / sqrt(nrow(holdout.ds)),
                      tolerance=1 / sqrt(nrow(holdout.ds)),
                      signal.names=NULL,
                      save.file=NULL,
                      verbose=FALSE) {
  if(is.null(train.ds) | is.null(holdout.ds)) {
    stop("At least train and holdout data sets must be provided")
  }
  if(is.simulated & is.null(signal.names)) {
    warning("No signal names provided")
  }
  n <- nrow(train.ds)
  d <- ncol(train.ds) - 1
  param.mtry <- rf.mtry
  if(is.null(rf.mtry)) {
    param.mtry <- floor(sqrt(d))
  } else {
    if((param.mtry < 1) | (param.mtry > d)) {
      stop(paste("mtry parameter", param.mtry, "out of range 1 < mtry < d"))
    }
  }
  ptm <- proc.time()
  if(verbose) cat("running Relief-F importance for training and holdout sets\n")
  important.scores <- getImportanceScores(train.set=train.ds,
                                          holdout.set=holdout.ds,
                                          label=label,
                                          imp.estimator=corelearn.estimator,
                                          verbose=verbose)
  if(verbose) cat("private EC importance:", corelearn.estimator,
                  "elapsed time:", (proc.time() - ptm)[3], "\n")

  q1.scores <- important.scores[[1]]
  q2.scores <- important.scores[[2]]
  att.names <- rownames(q1.scores)
  diff.scores <- abs(q1.scores - q2.scores)
  delta.q <- max(diff.scores)
  q1.scores.plot <- q1.scores

  T0 <- start.temp
  Tmin <- final.temp
  tau <- tau.param

  i <- 1
  myT <- T0

  fholds <- 0.5
  fvalidations <- 0.5
  ftrains <- 0.5
  correct.detect.ec <- vector(mode="numeric")
  oldAccuracy <- 0.5
  cur.vars.remain <- length(att.names)

  vars.remain <- c(0, cur.vars.remain)
  kept.atts <- att.names
  var.names <- list()

  if(verbose) cat("private EC optimization loop\n")
  num.updates <- 0
  while ((myT > Tmin) && (utils::tail(vars.remain, 1) > 2) &&
         (utils::tail(vars.remain, 2)[1] != utils::tail(vars.remain, 2)[2])) {

    diff <- diff.scores * (diff.scores > 10^(-3)) + 10^(-3) * (diff.scores < 10^(-3))
    PAs <- exp(-q1.scores / (2 * diff * myT))
    PAs <- PAs[kept.atts, ]
    sumPAs <- sum(PAs)
    scaled.PAs <- PAs / sum(PAs)
    cum.scaled.PAs <- cumsum(scaled.PAs)
    num.remv <- 1 # only remove 1 attribute
    prob.rands <- sort(stats::runif(num.remv, min=0, max=1))
    remv.atts <- kept.atts[prob.rands < cum.scaled.PAs][1]
    if (num.remv >= length(kept.atts)){
      kept.atts <- NULL
      break
    } else {
      kept.atts <- setdiff(kept.atts, remv.atts)
    }

    # Compute train and holdout accuracy for new S_A attributes:
    att.names <- kept.atts
    if((i %% update.freq) == 1) {
      num.updates <- num.updates + 1
      if(verbose) cat("step", i, "update", num.updates, "myT > Tmin",
                      myT, Tmin, "?\n")
      if(verbose) cat("\trecomputing scores with evaporated attributes removed\n")
      new.X_train <- train.ds[, c(kept.atts, label), drop=F]
      new.X_holdout <- holdout.ds[, c(kept.atts, label), drop=F]
      new.X_validation <- validation.ds[, c(kept.atts, label), drop=F]
      new.scores <- getImportanceScores(new.X_train, new.X_holdout, label=label)
      q1.scores <- new.scores[[1]]
      q2.scores <- new.scores[[2]]
      diff.scores <- abs(q1.scores - q2.scores)
      delta.q <- max(diff.scores)
      if(verbose) cat("\trunning randomForest\n")
      if(param.mtry >= ncol(new.X_train) - 1) {
        kept.atts <- NULL
        break
        # param.mtry <- floor(sqrt(ncol(new.X_train) - 1))
        #attempted.mtry <- param.mtry
        # param.mtry <-  floor((ncol(new.X_train) - 1) / 2)
        # cat("random forest mtry setting", attempted.mtry,
        #     ">= ", ncol(new.X_train) - 1,
        #     ", so setting to half the number of variables",
        #     param.mtry, "\n")
      }
      model.formula <- stats::as.formula(paste(label, "~.", sep=""))
      result.rf <- randomForest::randomForest(formula=model.formula,
                                              data=new.X_train,
                                              ntree=rf.ntree,
                                              mtry=param.mtry)
      ftrain <- 1 - mean(result.rf$confusion[,"class.error"])
      if(verbose) cat("\tpredict\n")
      holdout.pred <- stats::predict(result.rf, newdata=new.X_holdout)
      fholdout <- mean(holdout.pred == holdout.ds[, label])
      if(is.simulated) {
        validation.pred <- stats::predict(result.rf, newdata=new.X_validation)
        fvalidation <- mean(validation.pred == validation.ds[, label])
      } else {
        fvalidation <- 0
      }
      if(abs(ftrain - fholdout) < (threshold + stats::rnorm(1, 0, tolerance))) {
        fholdout <- ftrain
      } else {
        if(verbose) cat("\tadjust holdout with stats::rnorm\n")
        fholdout <- fholdout + stats::rnorm(1, 0, tolerance)
      }

      ftrains <- c(ftrains, ftrain)
      fholds <- c(fholds, fholdout)
      fvalidations <- c(fvalidations, fvalidation)

      if(verbose) cat("\tadjusting temperature\n")
      myT <- myT * exp(-1 / tau) # dropping T
      if(verbose) cat("\t", i - 1, ': ', myT, '\n')

      if(verbose) cat("\tcollecting results\n")
      cur.vars.remain <- length(att.names)
      vars.remain <- c(vars.remain, cur.vars.remain)
      var.names[[num.updates]] <- kept.atts
      if(is.simulated) {
        correct.detect.ec <- c(correct.detect.ec,
                               sum(var.names[[num.updates]] %in% signal.names))
      }
    }
    i <- i + 1
  }
  if(verbose) cat("private EC optimization loop elapsed time:",
                  (proc.time() - ptm)[3], "\n")

  vars.remain <- vars.remain[-1] # remove the first value 0
  fplots <- data.frame(vars.remain,
                       train.acc=ftrains,
                       holdout.acc=fholds,
                       validation.acc=fvalidations,
                       alg=1)
  fplots <- fplots[-1, ] # remove the first row
  melted.fs <- reshape2::melt(fplots, id=c("vars.remain", "alg"))
  if(!is.null(save.file)) {
    if(verbose) {
      cat("saving results to ", save.file, "\n")
    }
    save(fplots, melted.fs, correct.detect.ec, n, d, signal.names,
         threshold, tolerance, bias, file=save.file)
  }

  if(verbose) cat("privateEC elapsed time:", (proc.time() - ptm)[3], "\n")

  list(algo.acc=fplots,
       ggplot.data=melted.fs,
       correct=correct.detect.ec,
       elasped=(proc.time() - ptm)[3])
}

#' Original Thresholdout algorithm
#'
#' Original Thresholdout with Dwork’s linear classifier
#' (Dwork, et al., 2015)
#'
#' @param train.ds A data frame with training data and class labels
#' @param holdout.ds A data frame with holdout data and class labels
#' @param validation.ds A data frame with validation data and class labels
#' @param label A character vector of the class variable column name
#' @param is.simulated Is the data simulated (or real?)
#' @param update.freq A integer for the number of steps before update
#' @param pec.file A character vector filename of privateEC results
#' @param threshold A numeric, default 4 / sqrt(n) suggested in the
#'  thresholdout’s supplementary material (Dwork, et al.,2015)
#' @param tolerance A numeric, default 1 / sqrt(n) suggested in the
#'  thresholdout’s supplementary material (Dwork, et al.,2015)
#' @param signal.names A character vector of signal names in simulated data
#' @param save.file A character vector for results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list containing:
#' \describe{
#'   \item{algo.acc}{data frame of results, a row for each update}
#'   \item{ggplot.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @examples
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' temp.pec.file <- tempfile(pattern="pEc_temp", tmpdir=tempdir())
#'
#' sim.data <- createSimulation(num.vars=num.variables,
#'                              n=num.samples,
#'                              sim.type="mainEffect",
#'                              verbose=FALSE)
#' pec.results <- privateEC(train.ds=sim.data$train,
#'                          holdout.ds=sim.data$holdout,
#'                          validation.ds=sim.data$validation,
#'                          label=sim.data$class.label,
#'                          is.simulated=TRUE,
#'                          signal.names=sim.data$signal.names,
#'                          save.file=temp.pec.file,
#'                          verbose=FALSE)
#' por.results <- originalThresholdout(train.ds=sim.data$train,
#'                                    holdout.ds=sim.data$holdout,
#'                                    validation.ds=sim.data$validation,
#'                                    label=sim.data$class.label,
#'                                    is.simulated=TRUE,
#'                                    signal.names=sim.data$signal.names,
#'                                    pec.file=temp.pec.file,
#'                                    verbose=FALSE)
#' file.remove(temp.pec.file)
#' @family classification
#' @export
originalThresholdout <- function(train.ds=NULL, holdout.ds=NULL, validation.ds=NULL,
                                label="phenos",
                                is.simulated=TRUE,
                                update.freq=50,
                                pec.file=NULL,
                                threshold=4 / sqrt(nrow(holdout.ds)),
                                tolerance=1 / sqrt(nrow(holdout.ds)),
                                signal.names=NULL,
                                save.file=NULL,
                                verbose=FALSE) {
  if(is.null(train.ds) | is.null(holdout.ds)) {
    stop("At least train and holdout data sets must be provided")
  }
  if(is.simulated & is.null(signal.names)) {
    stop("No signal names provided")
  }
  if(is.null(pec.file)) {
    stop("No previous privateEC results file in pec.file argument")
  }
  if(!file.exists(pec.file)) {
    stop("privateEC results file expected in pec.file argument:", pec.file)
  }

  ptm <- proc.time()

  n <- nrow(train.ds)
  d <- ncol(train.ds) - 1

  # compare with original privacy::
  if(verbose) cat("loading resuls from privacy EC", pec.file, "\n")
  load(pec.file)

  predictors.train <- as.matrix(train.ds[, 1:d])
  predictors.holdout <- as.matrix(holdout.ds[, 1:d])
  if(is.simulated) {
    predictors.validation <- as.matrix(validation.ds[, 1:d])
  }
  train.pheno <- train.ds[, label]
  holdout.pheno <- holdout.ds[, label]
  if(is.simulated) {
    validation.pheno <- validation.ds[, label]
    validation.pheno <- as.numeric(levels(validation.pheno))[validation.pheno]
  }
  train.pheno <- as.numeric(levels(train.pheno))[train.pheno]
  holdout.pheno <- as.numeric(levels(holdout.pheno))[holdout.pheno]

  trainanswers <- (t(predictors.train) %*% train.pheno) / nrow(train.ds)
  holdoutanswers <- (t(predictors.holdout) %*% holdout.pheno) / nrow(holdout.ds)
  diffs <- abs(trainanswers - holdoutanswers)
  noise <- stats::rnorm(d, 0, tolerance)
  abovethr <- diffs > threshold + noise
  holdoutanswers[!abovethr] <- trainanswers[!abovethr]
  holdoutanswers[abovethr] <- (holdoutanswers +
                               stats::rnorm(d, 0, tolerance))[abovethr]
  trainpos <- trainanswers > 1 / sqrt(n)
  holdoutpos <- holdoutanswers > 1 / sqrt(n)
  trainneg <- trainanswers < -1 / sqrt(n)
  holdoutneg <- holdoutanswers < -1 / sqrt(n)

  selected <- (trainpos & holdoutpos) | (trainneg & holdoutneg)
  trainanswers[!selected] <- 0
  sortanswers <- order(abs(trainanswers))
  # vars.remain <- c(0,10,20,30,45,70,100,150,200,250,300,400,500)
  alg.steps <- seq(d + 1, 1, -update.freq)[-1]
  vars.remain <- alg.steps
  numks <- length(vars.remain)
  noisy_vals <- matrix(-1, numks, 3)
  correct.detect.ori <- vector(mode="numeric")
  var.names <- list()

  for (i in 1:numks){
    k <- vars.remain[i]
    topk <- utils::tail(sortanswers, k)
    weights <- matrix(0, d, 1)
    weights[topk] <- sign(trainanswers[topk])
    ftrain <- mean((sign(predictors.train %*% weights)) == train.pheno)
    fholdout <- mean((sign(predictors.holdout %*% weights)) == holdout.pheno)
    if(abs(ftrain - fholdout) < threshold + stats::rnorm(1, 0, tolerance)){
      fholdout <- ftrain
    } else {
      fholdout <- fholdout + stats::rnorm(1, 0, tolerance)
    }
    if(is.simulated) {
      fvalidation <- mean((sign(predictors.validation %*% weights)) == validation.pheno)
    } else {
      fvalidation <- 0
    }
    if(k == 0) {
      ftrain <- 0.5
      fholdout <- 0.5
      fvalidation <- 0.5
    }
    noisy_vals[i,] <- c(ftrain, fholdout, fvalidation)
    var.names[[i]] <- colnames(train.ds)[topk]
    if(is.simulated) {
      correct.detect.ori <- c(correct.detect.ori,
                              sum(var.names[[i]] %in% signal.names))
    }
  }

  colnames(noisy_vals) <- c("train.acc", "holdout.acc", "validation.acc")
  pplots <- data.frame(vars.remain, noisy_vals, alg=2)
  melted.ps <- reshape2::melt(pplots, id=c("vars.remain", "alg"))

  if(!is.null(save.file)) {
    if(verbose) cat("saving to", save.file, "\n")
    save(pplots, melted.ps, correct.detect.ori, file=save.file)
  }

  if(verbose) cat("originalThresholout elapsed time:", (proc.time() - ptm)[3], "\n")

  list(algo.acc=pplots,
       ggplot.data=melted.ps,
       correct=correct.detect.ori,
       elasped=(proc.time() - ptm)[3])
}

#' Private random forests algorithm
#'
#' Random Forest Thresholdout, which is TO with the feature selection
#' and classifier replaced with Random Forest.
#'
#' @param train.ds A data frame with training data and class labels
#' @param holdout.ds A data frame with holdout data and class labels
#' @param validation.ds A data frame with validation data and class labels
#' @param label A character vector of the class variable column name
#' @param is.simulated Is the data simulated (or real?)
#' @param rf.importance.measure A character vector for the random forest importance measure
#' @param rf.ntree An integer the number of trees in the random forest
#' @param rf.mtry An integer the number of variables sampled at each random forest node split
#' @param pec.file A character vector filename of privateEC results
#' @param update.freq A integer for the number of steps before update
#' @param threshold A numeric, default 4 / sqrt(n) suggested in the
#'  thresholdout’s supplementary material (Dwork, et al.,2015)
#' @param tolerance A numeric, default 1 / sqrt(n) suggested in the
#'  thresholdout’s supplementary material (Dwork, et al.,2015)
#' @param signal.names A character vector of signal names in simulated data
#' @param save.file A character vector for results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list containing:
#' \describe{
#'   \item{algo.acc}{data frame of results, a row for each update}
#'   \item{ggplot.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @examples
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' temp.pec.file <- tempfile(pattern="pEc_temp", tmpdir=tempdir())
#'
#' sim.data <- createSimulation(num.vars=num.variables,
#'                              n=num.samples,
#'                              sim.type="mainEffect",
#'                              verbose=FALSE)
#' pec.results <- privateEC(train.ds=sim.data$train,
#'                          holdout.ds=sim.data$holdout,
#'                          validation.ds=sim.data$validation,
#'                          label=sim.data$class.label,
#'                          is.simulated=TRUE,
#'                          signal.names=sim.data$signal.names,
#'                          save.file=temp.pec.file,
#'                          verbose=FALSE)
#' prf.results <- privateRF(train.ds=sim.data$train,
#'                          holdout.ds=sim.data$holdout,
#'                          validation.ds=sim.data$validation,
#'                          label=sim.data$class.label,
#'                          is.simulated=TRUE,
#'                          signal.names=sim.data$signal.names,
#'                          pec.file=temp.pec.file,
#'                          verbose=FALSE)
# file.remove(temp.pec.file)
#' @export
privateRF <- function(train.ds=NULL, holdout.ds=NULL, validation.ds=NULL,
                      label="phenos",
                      is.simulated=TRUE,
                      rf.importance.measure="MeanDecreaseGini",
                      rf.ntree=500,
                      rf.mtry=NULL,
                      pec.file=NULL,
                      update.freq=50,
                      threshold=4 / sqrt(nrow(holdout.ds)),
                      tolerance=1 / sqrt(nrow(holdout.ds)),
                      signal.names=NULL,
                      save.file=NULL,
                      verbose=FALSE) {
  if(is.null(train.ds) | is.null(holdout.ds)) {
    stop("At least train and holdout data sets must be provided")
  }
  n <- nrow(train.ds)
  d <- ncol(train.ds) - 1
  param.mtry <- rf.mtry
  if(is.null(rf.mtry)) {
    param.mtry <- floor(sqrt(d))
  }
  if(param.mtry < 1 | param.mtry > d) {
    stop(paste("mtry parameter", param.mtry, "out of range 1 < mtry < d"))
  }
  if(is.simulated & is.null(signal.names)) {
    stop("privateRF: No signal names provided")
  }
  if(is.null(pec.file)) {
    stop("privateRF: No previous privateEC results file in pec.file argument")
  }
  if(!file.exists(pec.file)) {
    stop("privateRF: Previous privateEC results file in pec.file argument:", pec.file)
  }

  ptm <- proc.time()

  if(verbose) cat("loading resuls from privacy EC", pec.file, "\n")
  load(pec.file)

  predictors.train <- as.matrix(train.ds[, 1:d])
  predictors.holdout <- as.matrix(holdout.ds[, 1:d])
  if(is.simulated) {
    predictors.validation <- as.matrix(validation.ds[, 1:d])
  }

  train.rf <- randomForest::randomForest(x=predictors.train,
                                         y=train.ds[, label],
                                         ntree=rf.ntree,
                                         mtry=param.mtry,
                                         importance=T)
  train.imp <- train.rf$importance[, rf.importance.measure]
  holdout.rf <- randomForest::randomForest(x=predictors.holdout,
                                           y=holdout.ds[, label],
                                           ntree=rf.ntree,
                                           mtry=param.mtry,
                                           importance=T)
  holdout.imp <- holdout.rf$importance[, rf.importance.measure]
  trainanswers <- train.imp
  holdoutanswers <- holdout.imp

  diffs <- abs(trainanswers - holdoutanswers)
  noise <- stats::rnorm(d, 0, tolerance)
  abovethr <- diffs > threshold + noise
  holdoutanswers[!abovethr] <- trainanswers[!abovethr]
  holdoutanswers[abovethr] <- (holdoutanswers + stats::rnorm(d, 0, tolerance))[abovethr]
  trainpos <- trainanswers > 1 / sqrt(n)
  holdoutpos <- holdoutanswers > 1 / sqrt(n)
  trainneg <- trainanswers < -1 / sqrt(n)
  holdoutneg <- holdoutanswers < -1 / sqrt(n)

  selected <- (trainpos & holdoutpos) | (trainneg & holdoutneg)
  trainanswers[!selected] <- 0
  sortanswers <- order(abs(trainanswers))
  alg.steps <- seq(d, 1, -update.freq)[-1]
  vars.remain <- alg.steps
  numks <- length(vars.remain)
  noisy_vals <- matrix(-1, numks, 3)
  correct.detect.rf <- vector(mode="numeric")
  var.names <- list()
  if(verbose) cat("loop for", numks, ": ")
  for (i in 1:numks) {
    if(verbose) cat(i, " ")
    k <- vars.remain[i]
    if (k == 0) {
      ftrain <- 0.5
      fholdout <- 0.5
      fvalidation <- 0.5
    } else {
      topk <- utils::tail(sortanswers, k)
      var.names[[i]] <- colnames(train.ds)[topk]
      # rf classifier:
      data.k <- train.ds[, topk, drop=F]
      if(param.mtry >= ncol(data.k)) {
        break
      }
      rf.model.k <- randomForest::randomForest(x=data.k,
                                               y=train.ds[, label],
                                               ntree=rf.ntree,
                                               mtry=param.mtry)
      # ytrain.pred <- stats::predict(rf.model.k, newdata=train.ds)
      yholdout.pred <- stats::predict(rf.model.k, newdata=holdout.ds)
      if(is.simulated) {
        yvalidation.pred <- stats::predict(rf.model.k, newdata=validation.ds)
      }
      ftrain <- 1 - mean(rf.model.k$confusion[, "class.error"])
      # ftrain <- mean(ytrain.pred == train.ds[, label])
      # ftrain <- 1- rf.model.k$prediction.error
      fholdout <- mean(yholdout.pred == holdout.ds[, label])

      if (abs(ftrain - fholdout) < threshold + stats::rnorm(1, 0, tolerance)) {
        fholdout <- ftrain
      } else {
        fholdout <- fholdout + stats::rnorm(1, 0, tolerance)
      }
      #   fvalidation <- sum((sign(predictors.validation %*% weights)) == pheno.validation)/n
      if(is.simulated) {
        fvalidation <- mean(yvalidation.pred == validation.ds[, label])
      } else {
        fvalidation <- 0
      }
    }
    noisy_vals[i,] <- c(ftrain, fholdout, fvalidation)
    if(is.simulated) {
      correct.detect.rf <- c(correct.detect.rf,
                           sum(var.names[[i]] %in% signal.names))
    }
  }
  if(verbose) cat(i, "\n")

  colnames(noisy_vals) <- c("train.acc", "holdout.acc", "validation.acc")
  rfplots <- data.frame(vars.remain, noisy_vals, alg=3)

  melted.rfs <- reshape2::melt(rfplots, id=c("vars.remain", "alg"))

  if(!is.null(save.file)) {
    if(verbose) cat("saving plot data", save.file, "\n")
    save(rfplots, melted.rfs, correct.detect.rf, file=save.file)
  }

  if(verbose) cat("privateRF elasped time:", (proc.time() - ptm)[3], "\n")

  list(algo.acc=rfplots,
       ggplot.data=melted.rfs,
       correct=correct.detect.rf,
       elasped=(proc.time() - ptm)[3])
}

#' Standard random forests algorithm serves as a baseline model
#'
#' @param train.ds A data frame with training data and class labels
#' @param holdout.ds A data frame with holdout data and class labels
#' @param validation.ds A data frame with validation data and class labels
#' @param label A character vector of the class variable column name
#' @param rf.ntree An integer the number of trees in the random forest
#' @param rf.mtry An integer the number of variables sampled at each random forest node split
#' @param is.simulated Is the data simulated (or real?)
#' @param signal.names A character vector of signal names in simulated data
#' @param save.file A character vector for results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list containing:
#' \describe{
#'   \item{algo.acc}{data frame of results, a row for each update}
#'   \item{ggplot.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @examples
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' sim.data <- createSimulation(num.vars=num.variables, n=num.samples,
#'                              sim.type="mainEffect", verbose=FALSE)
#' rra.results <- standardRF(train.ds=sim.data$train,
#'                           holdout.ds=sim.data$holdout,
#'                           validation.ds=sim.data$validation,
#'                           label=sim.data$class.label,
#'                           is.simulated=TRUE,
#'                           verbose=FALSE,
#'                           signal.names=sim.data$signal.names)
#' @family classification
#' @export
standardRF <- function(train.ds=NULL, holdout.ds=NULL, validation.ds=NULL,
                       label="phenos",
                       rf.ntree=500,
                       rf.mtry=NULL,
                       is.simulated=TRUE,
                       signal.names=NULL,
                       save.file=NULL,
                       verbose=FALSE) {
  if(is.null(train.ds) | is.null(holdout.ds)) {
    stop("At least train and holdout data sets must be provided")
  }
  n <- nrow(train.ds)
  d <- ncol(train.ds) - 1
  param.mtry <- rf.mtry
  if(is.null(rf.mtry)) {
    param.mtry <- floor(sqrt(d))
  }
  if(param.mtry < 1 | param.mtry > d) {
    stop(paste("mtry parameter", param.mtry, "out of range 1 < mtry < d"))
  }
  if(is.simulated & is.null(signal.names)) {
    stop("regularRF: No signal names provided")
  }
  ptm <- proc.time()
  bag.simul <- randomForest::randomForest(stats::as.formula(paste(label, "~ .", sep="")),
                                          data=rbind(train.ds, holdout.ds),
                                          ntree=rf.ntree,
                                          mtry=param.mtry,
                                          importance=T)
  # rf.holdo.accu <- 1- bag.simul$prediction.error
  rf.holdout.accu <- 1 - mean(bag.simul$confusion[, "class.error"])
  if(is.simulated) {
    rf.pred <- stats::predict(bag.simul, newdata=validation.ds)
    rf.validation.accu <- mean(rf.pred == validation.ds[, label])
  } else {
    rf.validation.accu  <- 0
  }
  if(verbose) cat("accuracies", rf.holdout.accu, rf.validation.accu, "\n")
  if(!is.null(save.file)) {
    save(rf.validation.accu, rf.holdout.accu, file=save.file)
  }

  rRaplots <- data.frame(vars.remain=ncol(train.ds),
                         train.acc=1,
                         holdout.acc=rf.holdout.accu,
                         validation.acc=rf.validation.accu,
                         alg=4)
  melted.rra <- reshape2::melt(rRaplots, id=c("vars.remain", "alg"))

  if(verbose) cat("regularRF elapsed time:", (proc.time() - ptm)[3], "\n")

  list(algo.acc=rRaplots,
       ggplot.data=melted.rra,
       correct=ifelse(is.simulated, ncol(train.ds), 0),
       elasped=(proc.time() - ptm)[3])
}

#' Call C++ inbix Evaporative Cooling Privacy
#'
#' Assumes the inbix executable is in the PATH.
#'
#' @param train.ds A data frame with training data and class labels
#' @param holdout.ds A data frame with holdout data and class labels
#' @param validation.ds A data frame with validation data and class labels
#' @param label A character vector of the class variable column name
#' @param is.simulated Is the data simulated (or real?)
#' @param bias A numeric for effect size in simulated signal variables
#' @param update.freq A integer for the number of steps before update
#' @param start.temp A numeric for EC starting temperature
#' @param final.temp A numeric for EC final temperature
#' @param tau.param A numeric for tau to control reduction schedule
#' @param rf.ntree An integer the number of trees in the random forest
#' @param rf.mtry An integer the number of variables sampled at each random forest node split
#' @param save.file A character vector for results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @note inbix must be instaled in the path
#' @return A list containing:
#' \describe{
#'   \item{algo.acc}{data frame of results, a row for each update}
#'   \item{ggplot.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @examples
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' sim.data <- createSimulation(num.vars=num.variables, n=num.samples,
#'                              sim.type="mainEffect", verbose=FALSE)
#' pec.results <- privateECinbix(train.ds=sim.data$train,
#'                               holdout.ds=sim.data$holdout,
#'                               validation.ds=sim.data$validation,
#'                               verbose=FALSE)
#' @family classification
#' @export
privateECinbix <- function(train.ds=NULL, holdout.ds=NULL, validation.ds=NULL,
                           label="phenos",
                           is.simulated=TRUE,
                           bias=0.4,
                           update.freq=50,
                           start.temp=0.1,
                           final.temp=10 ^ (-5),
                           tau.param=100,
                           rf.ntree=500,
                           rf.mtry=NULL,
                           save.file=NULL,
                           verbose=FALSE) {
  # check for inbix in the PATH or stop with an error
  if(Sys.which("inbix") == "") {
    stop("inbix is not installed or not in the PATH")
  }
  # check input parameters
  if(is.null(train.ds) | is.null(holdout.ds)) {
    stop("At least train and holdout data sets must be provided")
  }
  d <- ncol(train.ds) - 1
  param.mtry <- rf.mtry
  if(is.null(rf.mtry)) {
    param.mtry <- floor(sqrt(d))
  }
  if(param.mtry < 1 | param.mtry > d) {
    stop(paste("mtry parameter", param.mtry, "out of range 1 < mtry < d"))
  }
  ptm <- proc.time()
  unique.sim.prefix <- paste(save.file, bias, sep="_")
  correct.detect.inbix <- vector(mode="numeric")
  # write simple tab-separated vales (tsv) files
  if(verbose) cat("Writing ", unique.sim.prefix, ".sim.(train|holdout|validation).tab files\n")
  utils::write.table(train.ds, paste(unique.sim.prefix, ".sim.train.tab", sep=""), sep="\t",
                     row.names=FALSE, col.names=TRUE, quote=FALSE)
  utils::write.table(holdout.ds, paste(unique.sim.prefix, ".sim.holdout.tab", sep=""), sep="\t",
                     row.names=FALSE, col.names=TRUE, quote=FALSE)
  utils::write.table(validation.ds, paste(unique.sim.prefix, ".sim.test.tab", sep=""), sep="\t",
                     row.names=FALSE, col.names=TRUE, quote=FALSE)
  if(verbose) cat("Running inbix privacy EC on saved simulation data sets\n")
  system.command <- paste("OMP_NUM_THREADS=4 inbix ",
                          "--ec-privacy ",
                          "--depvarname ", label, " ",
                          "--ec-privacy-train-file ", unique.sim.prefix, ".sim.train.tab ",
                          "--ec-privacy-holdout-file ", unique.sim.prefix, ".sim.holdout.tab ",
                          "--ec-privacy-test-file ", unique.sim.prefix, ".sim.test.tab ",
                          "--ec-privacy-start-temp ", start.temp, " ",
                          "--ec-privacy-final-temp ", final.temp, " ",
                          "--ec-privacy-tau ", tau.param, " ",
                          "--ec-privacy-update-frequency ", update.freq, " ",
                          "--ntree ", rf.ntree, " ",
                          "--mtry ", param.mtry, " ",
                          "--out ", unique.sim.prefix, ".privateec", sep="")
  if(verbose) cat(system.command, "\n")
  system(system.command, intern=FALSE, ignore.stdout=TRUE, ignore.stderr=TRUE)
  # if(verbose) cat("Reading inbix privacy EC attributes used results\n")
  # attrs.file <- paste(unique.sim.prefix, ".privateec.selected.attributes.tab", sep="")
  # if(verbose) cat(attrs.file, "\n")
  # attrs.table <- utils::read.table(attrs.file, sep="\t", header=FALSE, stringsAsFactors=FALSE)

  if(verbose) cat("Reading inbix privacy EC algorithm run details\n")
  iters.file <- paste(unique.sim.prefix, ".privateec.privacy.iterations.tab", sep="")
  if(verbose) cat(iters.file, "\n")
  iters.table <- utils::read.table(iters.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  if(is.simulated) {
    correct.detect.inbix <- as.integer(iters.table$Correct)
  }
  fxplots <- data.frame(vars.remain=as.integer(iters.table$Keep),
                        ftrain=iters.table$TrainAcc,
                        fholdout=iters.table$HoldoutAcc,
                        fvalidation=iters.table$TestAcc,
                        alg=5)
  melted.fx <- reshape2::melt(fxplots, id=c("vars.remain", "alg"))
  if(verbose) cat("Cleaning up inbix private EC results files\n")
  inbix.temp.files <- c(paste(unique.sim.prefix, ".sim.train.tab", sep=""),
                        paste(unique.sim.prefix, ".sim.holdout.tab", sep=""),
                        paste(unique.sim.prefix, ".sim.test.tab", sep=""),
                        paste(unique.sim.prefix, ".privateec.privacy.iterations.tab", sep=""),
                        paste(unique.sim.prefix, ".privateec.log", sep=""))
  #file.remove(c(inbix.conv.files, rf.files, attrs.file, iters.file))
  file.remove(inbix.temp.files)
  if(!is.null(save.file)) {
    if(verbose) {
      cat("saving results to ", save.file, "\n")
    }
    save(fxplots, melted.fx, correct.detect.inbix, bias,
         train.ds, holdout.ds, validation.ds, file=save.file)
  }
  if(verbose) cat("privacyECinbix elapsed time:", (proc.time() - ptm)[3], "\n")

  list(algo.acc=fxplots,
       ggplot.data=melted.fx,
       correct=correct.detect.inbix,
       elasped=(proc.time() - ptm)[3])
}
