# -----------------------------------------------------------------------------
# classification.R - Trang Le and Bill White - Fall 2016/Spring 2017
#
# Privacy preserving evaporative cooling and comparison classification
# algorithms used in the paper XXXXX.
# -----------------------------------------------------------------------------

#' Compute and return importance scores (Relief-F scores)
#'
#' @param train.set A training data frame with last column as class
#' @param holdo.set A holdout data frame with last column as class
#' @param imp.estimator A character vector CORElearn attribute importance estimator
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list with two data frames representing the importance scores
#' (Relief-F scores) for the train and holdout data sets.
getImportanceScores <- function(train.set=NULL,
                                holdo.set=NULL,
                                imp.estimator="ReliefFbestK",
                                verbose=FALSE) {
  if(is.null(train.set)) {
    stop("getImportanceScores: No training data set provided")
  }
  if(is.null(holdo.set)) {
    stop("getImportanceScores: No holdout data set provided")
  }
  if(verbose) cat("\ttrain\n")
  train.relief <- CORElearn::attrEval("pheno", data=train.set,
                                      estimator=imp.estimator)
  if(verbose) cat("\tholdout\n")
  holdo.relief <- CORElearn::attrEval("pheno", data=holdo.set,
                                      estimator=imp.estimator)
  list(data.frame(train.relief), data.frame(holdo.relief))
}

#' Private Evaporative Cooling feature selection and classification
#'
#' @param data.sets A list of train, holdout and test data frames
#' @param is.simulated Is the data simulated (or real?)
#' @param n An integer for the number of samples
#' @param d An integer for the number of variables
#' @param shortname A character vector of a parameters separated by '_'
#' @param bias A numeric for bias in signal variables simulation
#' @param type A character vector of the type of simulation: sva|er|inte
#' @param myrun A character vector of a unique run identifier
#' @param update.freq An integer the number of steps before update
#' @param corelearn.estimator CORElearn Relief-F estimator
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
#'   \item{plots.data}{data frame of results, a row for each update}
#'   \item{melted.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @note Within thresholdout, we choose a threshold of 4 / sqrt(n) and
#' tolerance of 1 / sqrt(n) as suggested in the thresholdout’s supplementary
#' material (Dwork, et al., 2015).
#' @export
privateEC <- function(data.sets=NULL,
                      is.simulated=TRUE,
                      n=100,
                      d=100,
                      shortname="paramstring",
                      bias=0.4,
                      type="sva",
                      myrun="001",
                      update.freq=50,
                      corelearn.estimator="ReliefFbestK",
                      ntree=100,
                      start.temp=0.1,
                      final.temp=10 ^ (-5),
                      tau.param=100,
                      threshold=4 / sqrt(n),
                      tolerance=1 / sqrt(n),
                      signal.names=NULL,
                      save.file=NULL,
                      verbose=FALSE) {
  if(is.null(data.sets)) {
    stop("privateEC: No data sets provided as first argument")
  }
  if(is.simulated & is.null(signal.names)) {
    stop("privateEC: No signal names provided")
  }
  ptm <- proc.time()
  if(verbose) cat("running Relief-F importance for training and holdout sets\n")
  important.scores <- getImportanceScores(train.set=data.sets$train,
                                          holdo.set=data.sets$holdout,
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
  ftests <- 0.5
  ftrains <- 0.5
  correct.detect.ec <- vector(mode="numeric")
  oldAccuracy <- 0.5
  num.att <- length(att.names)

  num.atts <- c(0, num.att)
  kept.atts <- att.names
  var.names <- list()

  if(verbose) cat("private EC optimization loop\n")
  num.updates <- 0
  while ((myT > Tmin) && (tail(num.atts,1) > 2) &&
         (tail(num.atts,2)[1] != tail(num.atts,2)[2])) {

    diff <- diff.scores * (diff.scores > 10^(-3)) + 10^(-3) * (diff.scores < 10^(-3))
    PAs <- exp(-q1.scores / (2 * diff * myT))
    PAs <- PAs[kept.atts,]
    sumPAs <- sum(PAs)
    scaled.PAs <- PAs / sum(PAs)
    cum.scaled.PAs <- cumsum(scaled.PAs)
    num.remv <- 1 # only remove 1 attribute
    prob.rands <- sort(runif(num.remv, min=0, max=1))
    remv.atts <- kept.atts[prob.rands < cum.scaled.PAs][1]
    if (num.remv >= length(kept.atts)){
      kept.atts <- NULL
      break
    } else {
      kept.atts <- setdiff(kept.atts, remv.atts)
    }

    # Compute train and holdo accuracy for new S_A attributes:
    att.names <- kept.atts
    if ((i %% update.freq) == 1) {
      num.updates <- num.updates + 1
      if(verbose) cat("step", i, "update", num.updates, "myT > Tmin",
                      myT, Tmin, "?\n")
      if(verbose) cat("\trecomputing scores with evaporated attributes removed\n")
      new.X_train <- data.sets$train[, c(kept.atts, "pheno"), drop=F]
      new.X_holdo <- data.sets$holdout[, c(kept.atts, "pheno"), drop=F]
      new.X_test <- data.sets$test[, c(kept.atts, "pheno"), drop=F]
      new.scores <- getImportanceScores(new.X_train, new.X_holdo)
      q1.scores <- new.scores[[1]]
      q2.scores <- new.scores[[2]]
      diff.scores <- abs(q1.scores - q2.scores)
      delta.q <- max(diff.scores)
      if(verbose) cat("\trunning randomForest\n")
      result.rf <- randomForest::randomForest(pheno ~ ., data=new.X_train)
      ftrain <- 1 - mean(result.rf$confusion[,"class.error"])
      if(verbose) cat("\tpredict\n")
      holdo.pred <- predict(result.rf, newdata=new.X_holdo)
      fholdo <- mean(holdo.pred == data.sets$holdout$pheno)
      if(is.simulated) {
        test.pred <- predict(result.rf, newdata=new.X_test)
        ftest <- mean((test.pred == data.sets$test$pheno))
      } else {
        ftest <- 0
      }
      if(abs(ftrain - fholdo) < (threshold + rnorm(1, 0, tolerance))) {
        fholdo <- ftrain
      } else {
        if(verbose) cat("\tadjust holdout with rnorm\n")
        fholdo <- fholdo + rnorm(1, 0, tolerance)
      }

      fholds <- c(fholds, fholdo)
      ftests <- c(ftests, ftest)
      ftrains <- c(ftrains, ftrain)

      if(verbose) cat("\tadjusting temperature\n")
      myT <- myT * exp(-1 / tau) # dropping T
      if(verbose) cat("\t", i - 1, ': ', myT, '\n')

      if(verbose) cat("\tcollecting results\n")
      num.att <- length(att.names)
      num.atts <- c(num.atts, num.att)
      var.names[[i]] <- kept.atts
      if(is.simulated) {
        correct.detect.ec <- c(correct.detect.ec,
                               sum(var.names[[i]] %in% signal.names))
      }
    }
    i <- i + 1
  }
  if(verbose) cat("private EC optimization loop elapsed time:",
                  (proc.time() - ptm)[3], "\n")

  num.atts <- num.atts[-1] # remove the first value 0
  fplots <- data.frame(num.atts,
                       ftrain=ftrains,
                       fholdo=fholds,
                       ftest=ftests,
                       alg=1)
  fplots <- fplots[-1, ] # remove the first row
  melted.fs <- reshape2::melt(fplots, id=c("num.atts", "alg"))
  num.atts <- num.atts[-1]
  if(!is.null(save.file)) {
    if(verbose) {
      cat("saving results to ", save.file, "\n")
    }
    save(fplots, melted.fs, correct.detect.ec, n, d, signal.names,
         bias, shortname, threshold, tolerance, bias, file=save.file)
  }
  cat("privateEC elapsed time:", (proc.time() - ptm)[3], "\n")

  list(plots.data=fplots,
       melted.data=melted.fs,
       correct=correct.detect.ec,
       elasped=(proc.time() - ptm)[3])
}

#' Original privacy algorithm
#'
#' Original Thresholdout with Dwork’s linear classifier
#' (Dwork, et al., 2015)
#'
#' @param data.sets A list of train, holdout and test data frames
#' @param is.simulated Is the data simulated (or real?)
#' @param n An integer for the number of samples
#' @param d An integer for the number of variables
#' @param shortname A character vector of a parameters separated by '_'
#' @param type A character vector of the type of simulation: sva|er|inte
#' @param myrun A character vector of a unique run identifier
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
#'   \item{plots.data}{data frame of results, a row for each update}
#'   \item{melted.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @export
originalPrivacy <- function(data.sets=NULL,
                            is.simulated=TRUE,
                            n=100,
                            d=100,
                            shortname="paramstring",
                            type="sva",
                            myrun="001",
                            update.freq=50,
                            pec.file=NULL,
                            threshold=4 / sqrt(n),
                            tolerance=1 / sqrt(n),
                            signal.names=NULL,
                            save.file=NULL,
                            verbose=FALSE) {
  if(is.null(data.sets)) {
    stop("originalPrivacy: No data sets provided as data.sets argument")
  }
  if(is.simulated & is.null(signal.names)) {
    stop("originalPrivacy: No signal names provided")
  }
  if(is.null(pec.file)) {
    stop("originalPrivacy: No previous privateEC results file in pec.file argument")
  }
  if(!file.exists(pec.file)) {
    stop("originalPrivacy: privateEC results file expected in pec.file argument:", pec.file)
  }

  ptm <- proc.time()

  # compare with original privacy::
  if(verbose) cat("loading resuls from privacy EC", pec.file, "\n")
  load(pec.file)

  predictors.train <- as.matrix(data.sets$train[, 1:(d - 1)])
  predictors.holdo <- as.matrix(data.sets$holdout[, 1:(d - 1)])
  if(is.simulated) {
    predictors.test <- as.matrix(data.sets$test[, 1:(d - 1)])
  }
  train.pheno <- data.sets$train[, d]
  holdo.pheno <- data.sets$holdout[, d]
  if(is.simulated) {
    test.pheno <- data.sets$test[, d]
    test.pheno <- as.numeric(levels(test.pheno))[test.pheno]
  }
  train.pheno <- as.numeric(levels(train.pheno))[train.pheno]
  holdo.pheno <- as.numeric(levels(holdo.pheno))[holdo.pheno]

  trainanswers <- (t(predictors.train) %*% train.pheno) / nrow(data.sets$train)
  holdoanswers <- (t(predictors.holdo) %*% holdo.pheno) / nrow(data.sets$holdout)
  diffs <- abs(trainanswers - holdoanswers)
  noise <- rnorm(d - 1, 0, tolerance)
  abovethr <- diffs > threshold + noise
  holdoanswers[!abovethr] <- trainanswers[!abovethr]
  holdoanswers[abovethr] <- (holdoanswers +
                               rnorm(d - 1, 0, tolerance))[abovethr]
  trainpos <- trainanswers > 1 / sqrt(n)
  holdopos <- holdoanswers > 1 / sqrt(n)
  trainneg <- trainanswers < -1 / sqrt(n)
  holdoneg <- holdoanswers < -1 / sqrt(n)

  selected <- (trainpos & holdopos) | (trainneg & holdoneg)
  trainanswers[!selected] <- 0
  sortanswers <- order(abs(trainanswers))
  # num.atts <- c(0,10,20,30,45,70,100,150,200,250,300,400,500)
  alg.steps <- seq(d + 1, 1, -update.freq)[-1]
  num.atts <- alg.steps
  numks <- length(num.atts)
  noisy_vals <- matrix(-1, numks, 3)
  correct.detect.ori <- vector(mode="numeric")
  var.names <- list()

  for (i in 1:numks){
    k <- num.atts[i]
    topk <- tail(sortanswers, k)
    weights <- matrix(0, d - 1, 1)
    weights[topk] <- sign(trainanswers[topk])
    ftrain <- mean((sign(predictors.train %*% weights)) == train.pheno)
    fholdo <- mean((sign(predictors.holdo %*% weights)) == holdo.pheno)
    if(abs(ftrain - fholdo) < threshold + rnorm(1, 0, tolerance)){
      fholdo <- ftrain
    } else {
      fholdo <- fholdo + rnorm(1, 0, tolerance)
    }
    if(is.simulated) {
      ftest <- mean((sign(predictors.test %*% weights)) == test.pheno)
    } else {
      ftest <- 0
    }
    if(k == 0) {
      ftrain <- 0.5
      fholdo <- 0.5
      ftest <- 0.5
    }
    noisy_vals[i,] <- c(ftrain, fholdo, ftest)
    var.names[[i]] <- colnames(data.sets$train)[topk]
    if(is.simulated) {
      correct.detect.ori <- c(correct.detect.ori,
                            sum(var.names[[i]] %in% signal.names))
    }
  }

  colnames(noisy_vals) <- c("ftrain", "fholdo", "ftest")
  pplots <- data.frame(num.atts, noisy_vals, alg=2)
  melted.ps <- reshape2::melt(pplots, id=c("num.atts", "alg"))

  if(!is.null(save.file)) {
    if(verbose) cat("saving to", save.file, "\n")
    save(pplots, melted.ps, correct.detect.ori, file=save.file)
  }
  cat("originalPrivacy elapsed time:", (proc.time() - ptm)[3], "\n")

  list(plots.data=pplots,
       melted.data=melted.ps,
       correct=correct.detect.ori,
       elasped=(proc.time() - ptm)[3])
}

# -----------------------------------------------------------------------------
#' Private random forests algorithm
#'
#' Random Forest Thresholdout, which is TO with the feature selection
#' and classifier replaced with Random Forest.
#'
#' @param data.sets A list of train, holdout and test data frames
#' @param is.simulated Is the data simulated (or real?)
#' @param n An integer for the number of samples
#' @param d An integer for the number of variables
#' @param shortname A character vector of a parameters separated by '_'
#' @param type A character vector of the type of simulation: sva|er|inte
#' @param rf.importance.measure A character vector for the random forest importance measure
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
#'   \item{plots.data}{data frame of results, a row for each update}
#'   \item{melted.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @export
privateRF <- function(data.sets=NULL,
                      is.simulated=TRUE,
                      n=100,
                      d=100,
                      shortname="paramstring",
                      type="sva",
                      rf.importance.measure="MeanDecreaseGini",
                      pec.file=NULL,
                      update.freq=50,
                      threshold=4 / sqrt(n),
                      tolerance=1 / sqrt(n),
                      signal.names=NULL,
                      save.file=NULL,
                      verbose=FALSE) {
  if(is.null(data.sets)) {
    stop("privateRF: No data sets provided as first argument")
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

  predictors.train <- as.matrix(data.sets$train[, 1:(d - 1)])
  predictors.holdo <- as.matrix(data.sets$holdout[, 1:(d - 1)])
  if(is.simulated) {
    predictors.test <- as.matrix(data.sets$test[, 1:(d - 1)])
  }

  train.rf <- randomForest::randomForest(x=predictors.train,
                                         y=data.sets$train$pheno,
                                         importance=T)
  train.imp <- train.rf$importance[, rf.importance.measure]
  holdo.rf <- randomForest::randomForest(x=predictors.holdo,
                                         y=data.sets$holdout$pheno,
                                         importance=T)
  holdo.imp <- holdo.rf$importance[, rf.importance.measure]
  trainanswers <- train.imp
  holdoanswers <- holdo.imp

  diffs <- abs(trainanswers - holdoanswers)
  noise <- rnorm(d - 1, 0, tolerance)
  abovethr <- diffs > threshold + noise
  holdoanswers[!abovethr] <- trainanswers[!abovethr]
  holdoanswers[abovethr] <- (holdoanswers + rnorm(d - 1, 0, tolerance))[abovethr]
  trainpos <- trainanswers > 1 / sqrt(n)
  holdopos <- holdoanswers > 1 / sqrt(n)
  trainneg <- trainanswers < -1 / sqrt(n)
  holdoneg <- holdoanswers < -1 / sqrt(n)

  selected <- (trainpos & holdopos) | (trainneg & holdoneg)
  trainanswers[!selected] <- 0
  sortanswers <- order(abs(trainanswers))
  alg.steps <- seq(d + 1, 1, -update.freq)[-1]
  num.atts <- alg.steps
  numks <- length(num.atts)
  noisy_vals <- matrix(-1, numks, 3)
  correct.detect.rf <- vector(mode="numeric")
  var.names <- list()
  if(verbose) cat("loop for", numks, ": ")
  for (i in 1:numks) {
    if(verbose) cat(i, " ")
    k <- num.atts[i]
    if (k == 0) {
      ftrain <- 0.5
      fholdo <- 0.5
      ftest <- 0.5
    } else {
      topk <- tail(sortanswers, k)
      var.names[[i]] <- colnames(data.sets$train)[topk]
      # rf classifier:
      data.k <- data.sets$train[, topk, drop=F]
      rf.model.k <- randomForest::randomForest(x=data.k,
                                               y=data.sets$train$pheno)
      # ytrain.pred <- predict(rf.model.k, newdata=data.sets$train)
      yholdo.pred <- predict(rf.model.k, newdata=data.sets$holdout)
      if(is.simulated) {
        ytest.pred <- predict(rf.model.k, newdata=data.sets$test)
      }
      ftrain <- 1 - mean(rf.model.k$confusion[, "class.error"])
      # ftrain <- mean(ytrain.pred == data.sets$train$pheno)
      # ftrain <- 1- rf.model.k$prediction.error
      fholdo <- mean(yholdo.pred == data.sets$holdout$pheno)

      if (abs(ftrain - fholdo) < threshold + rnorm(1, 0, tolerance)) {
        fholdo <- ftrain
      } else {
        fholdo <- fholdo + rnorm(1, 0, tolerance)
      }
      #   ftest <- sum((sign(predictors.test %*% weights)) == pheno.test)/n
      if(is.simulated) {
        ftest <- mean(ytest.pred == data.sets$test$pheno)
      } else {
        ftest <- 0
      }
    }
    noisy_vals[i,] <- c(ftrain, fholdo, ftest)
    if(is.simulated) {
      correct.detect.rf <- c(correct.detect.rf,
                           sum(var.names[[i]] %in% signal.names))
    }
  }
  if(verbose) cat(i, "\n")

  colnames(noisy_vals) <- c("ftrain", "fholdo", "ftest")
  rfplots <- data.frame(num.atts, noisy_vals, alg=3)

  melted.rfs <- reshape2::melt(rfplots, id=c("num.atts", "alg"))

  if(!is.null(save.file)) {
    if(verbose) cat("saving plot data", save.file, "\n")
    save(rfplots, melted.rfs, correct.detect.rf, file=save.file)
  }

  cat("privateRF elasped time:", (proc.time() - ptm)[3], "\n")

  list(plots.data=rfplots,
       melted.data=melted.rfs,
       correct=correct.detect.rf,
       elasped=(proc.time() - ptm)[3])
}

# -----------------------------------------------------------------------------
#' Standard random forests algorithm serves as a baseline model
#'
#' @param data.sets A list of train, holdout and test data frames
#' @param is.simulated Is the data simulated (or real?)
#' @param shortname A character vector of a parameters separated by '_'
#' @param type A character vector of the type of simulation: sva|er|inte
#' @param signal.names A character vector of signal names in simulated data
#' @param save.file A character vector for results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list containing:
#' \describe{
#'   \item{plots.data}{data frame of results, a row for each update}
#'   \item{melted.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @export
standardRF <- function(data.sets=NULL,
                       is.simulated=TRUE,
                       shortname="paramstring",
                       type="sva",
                       signal.names=NULL,
                       save.file=NULL,
                       verbose=FALSE) {
  if(is.null(data.sets)) {
    stop("regularRF: No data sets provided as first argument")
  }
  if(is.simulated & is.null(signal.names)) {
    stop("regularRF: No signal names provided")
  }
  ptm <- proc.time()
  bag.simul <- randomForest::randomForest(pheno ~.,
                                          data=rbind(data.sets$train,
                                                     data.sets$holdout),
                                          importance=T)
  # rf.holdo.accu <- 1- bag.simul$prediction.error
  rf.holdo.accu <- 1 - mean(bag.simul$confusion[, "class.error"])
  if(is.simulated) {
    rf.pred <- predict(bag.simul, newdata=data.sets$test)
    rf.test.accu <- mean(rf.pred == data.sets$test$pheno)
  } else {
    rf.test.accu  <- 0
  }
  if(verbose) cat("accuracies", rf.holdo.accu, rf.test.accu, "\n")
  if(!is.null(save.file)) {
    save(rf.test.accu, rf.holdo.accu, file=save.file)
  }

  rRaplots <- data.frame(num.atts=ncol(data.sets$train),
                         ftrain=1,
                         fholdo=rf.holdo.accu,
                         ftest=rf.test.accu,
                         alg=4)
  melted.rra <- reshape2::melt(rRaplots, id=c("num.atts", "alg"))

  cat("regularRF elapsed time:", (proc.time() - ptm)[3], "\n")

  list(plots.data=rRaplots,
       melted.data=melted.rra,
       correct=ifelse(is.simulated, ncol(data.sets$train), 0),
       elasped=(proc.time() - ptm)[3])
}

#' Call C++ inbix Evaporative Cooling Privacy
#'
#' Assumes the inbix executable is in the PATH.
#'
#' @param data.sets A list of train, holdout and test data frames
#' @param is.simulated Is the data simulated (or real?)
#' @param n An integer for the number of samples
#' @param shortname A character vector of a parameters separated by '_'
#' @param bias A numeric for ?
#' @param type A character vector of the type of simulation: sva|er|inte
#' @param myrun A character vector of a unique run identifier
#' @param update.freq A integer for the number of steps before update
#' @param start.temp A numeric for EC starting temperature
#' @param final.temp A numeric for EC final temperature
#' @param tau.param A numeric for tau to control reduction schedule
#' @param save.file A character vector for results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @note inbix must be instaled in the path
#' @return A list containing:
#' \describe{
#'   \item{plots.data}{data frame of results, a row for each update}
#'   \item{melted.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @export
privateECinbix <- function(data.sets=NULL,
                           is.simulated=TRUE,
                           n=100,
                           shortname="paramstring",
                           bias=0.4,
                           type="sva",
                           myrun="001",
                           update.freq=50,
                           start.temp=0.1,
                           final.temp=10 ^ (-5),
                           tau.param=100,
                           save.file=NULL,
                           verbose=FALSE) {
  # ---------------------------------------------------------------------------
  if(is.null(data.sets)) {
    stop("privateECinbix: No data sets provided as first argument")
  }
  # check for inbix in the PATH or stop with an error
  if(Sys.which("inbix") == "") {
    stop("inbix is not installed or not in the PATH")
  }
  # ---------------------------------------------------------------------------
  ptm <- proc.time()
  base.sim.prefix <- paste(type, "_", shortname, sep="")
  unique.sim.prefix <- paste(type, "_", shortname, sep="")
  correct.detect.inbix <- vector(mode="numeric")
  # write simple tab-separated vales (tsv) files
  if(verbose) cat("Writing ", unique.sim.prefix, ".sim.(train|holdout|test).tab files\n")
  write.table(data.sets$train, paste(unique.sim.prefix, ".sim.train.tab", sep=""), sep="\t",
              row.names=FALSE, col.names=TRUE, quote=FALSE)
  write.table(data.sets$holdout, paste(unique.sim.prefix, ".sim.holdo.tab", sep=""), sep="\t",
              row.names=FALSE, col.names=TRUE, quote=FALSE)
  write.table(data.sets$test, paste(unique.sim.prefix, ".sim.test.tab", sep=""), sep="\t",
              row.names=FALSE, col.names=TRUE, quote=FALSE)
  # ---------------------------------------------------------------------------
  if(verbose) cat("Running inbix privacy EC on saved simulation data sets\n")
  system.command <- paste("OMP_NUM_THREADS=4 inbix ",
                          "--ec-privacy ",
                          "--depvarname pheno ",
                          "--ec-privacy-train-file ", unique.sim.prefix, ".sim.train.tab ",
                          "--ec-privacy-holdout-file ", unique.sim.prefix, ".sim.holdo.tab ",
                          "--ec-privacy-test-file ", unique.sim.prefix, ".sim.test.tab ",
                          "--ec-privacy-start-temp ", start.temp, " ",
                          "--ec-privacy-final-temp ", final.temp, " ",
                          "--ec-privacy-tau ", tau.param, " ",
                          "--ec-privacy-update-frequency ", update.freq, " ",
                          "--out ", unique.sim.prefix, ".privateec", sep="")
  if(verbose) cat(system.command, "\n")
  system(system.command, intern=FALSE, ignore.stdout=TRUE, ignore.stderr=TRUE)
  # ---------------------------------------------------------------------------
  # if(verbose) cat("Reading inbix privacy EC attributes used results\n")
  # attrs.file <- paste(unique.sim.prefix, ".privateec.selected.attributes.tab", sep="")
  # if(verbose) cat(attrs.file, "\n")
  # attrs.table <- read.table(attrs.file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  # ---------------------------------------------------------------------------
  if(verbose) cat("Reading inbix privacy EC algorithm run details\n")
  iters.file <- paste(unique.sim.prefix, ".privateec.privacy.iterations.tab", sep="")
  if(verbose) cat(iters.file, "\n")
  iters.table <- read.table(iters.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  #correct.detect.inbix <- iters.table[, 9]
  if(is.simulated) {
    correct.detect.inbix <- as.integer(iters.table$Correct)
  }
  fxplots <- data.frame(num.atts=as.integer(iters.table$Keep),
                        ftrain=iters.table$TrainAcc,
                        fholdo=iters.table$HoldoutAcc,
                        ftest=iters.table$TestAcc,
                        alg=5)
  melted.fx <- reshape2::melt(fxplots, id=c("num.atts", "alg"))
  # ---------------------------------------------------------------------------
  if(verbose) cat("Cleaning up inbix private EC results files\n")
  inbix.temp.files <- c(paste(unique.sim.prefix, ".sim.train.tab", sep=""),
                        paste(unique.sim.prefix, ".sim.holdo.tab", sep=""),
                        paste(unique.sim.prefix, ".sim.test.tab", sep=""),
                        paste(unique.sim.prefix, ".privateec.privacy.iterations.tab", sep=""),
                        paste(unique.sim.prefix, ".privateec.log", sep=""))
  #file.remove(c(inbix.conv.files, rf.files, attrs.file, iters.file))
  file.remove(inbix.temp.files)
  # ---------------------------------------------------------------------------
  if(!is.null(save.file)) {
    if(verbose) {
      cat("saving results to ", save.file, "\n")
    }
    save(fxplots, melted.fx, correct.detect.inbix, n, bias, shortname,
         type, data.sets$train, data.sets$holdout, data.sets$test,
         file=save.file)
  }
  # ---------------------------------------------------------------------------
  cat("privacyECinbix elapsed time:", (proc.time() - ptm)[3], "\n")

  list(plots.data=fxplots,
       melted.data=melted.fx,
       correct=correct.detect.inbix,
       elasped=(proc.time() - ptm)[3])
}
