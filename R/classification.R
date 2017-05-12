# classification.R - Trang Le and Bill White - Fall 2016/Spring 2017
# Privacy preserving evaporative cooling and comparison classification
# algorithms used in the Bioinformatics paper.

#' Compute and return importance scores (Relief-F scores)
#'
#' @param train.set A training data frame with last column as class
#' @param holdo.set A holdout data frame with last column as class
#' @param label A character vector of the class variable column name
#' @param imp.estimator A character vector CORElearn attribute importance estimator
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list with two data frames representing the importance scores
#' (Relief-F scores) for the train and holdout data sets.
getImportanceScores <- function(train.set=NULL,
                                holdo.set=NULL,
                                label="phenos",
                                imp.estimator="ReliefFbestK",
                                verbose=FALSE) {
  if(is.null(train.set)) {
    stop("getImportanceScores: No training data set provided")
  }
  if(is.null(holdo.set)) {
    stop("getImportanceScores: No holdout data set provided")
  }
  if(verbose) cat("\ttrain\n")
  train.relief <- CORElearn::attrEval(label, data=train.set, estimator=imp.estimator)
  if(verbose) cat("\tholdout\n")
  holdo.relief <- CORElearn::attrEval(label, data=holdo.set, estimator=imp.estimator)

  list(data.frame(train.relief), data.frame(holdo.relief))
}

#' Private Evaporative Cooling feature selection and classification
#'
#'
#' @param data.sets A list of train, holdout and validation data frames
#' @param label A character vector of the class variable column name
#' @param is.simulated Is the data simulated (or real?)
#' @param bias A numeric for effect size in simulated signal variables
#' @param sim.type A character vector of the type of simulation:
#' mainEffect/interactionErdos/interactionScalefree
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
#' @param shortname A character vector of a parameters separated by '_'
#' @param save.file A character vector for results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list with:
#' \describe{
#'   \item{plots.data}{data frame of results, a row for each update}
#'   \item{melted.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @examples
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' nbias <- pct.signals * num.variables
#' signals <- sprintf("gene%04d", 1:nbias)
#' sim.data <- createSimulation(d=num.variables, n=num.samples,
#'                              sim.type="mainEffect", verbose=FALSE)
#' ec.results <- privateEC(data.sets=sim.data, is.simulated=TRUE,
#'                         signal.names=signals, verbose=FALSE)
#' @note Within thresholdout, we choose a threshold of 4 / sqrt(n) and
#' tolerance of 1 / sqrt(n) as suggested in the thresholdout’s supplementary
#' material (Dwork, et al., 2015).
#' @references
#' Trang Le, W. K. Simmons, M. Misaki, B.C. White, J. Savitz, J. Bodurka,
#' and B. A. McKinney. “Privacy preserving evaporative cooling feature
#' selection and classification with Relief-F and Random Forests,”
#' Bioinformatics. Accepted. https://doi.org/10.1093/bioinformatics/btx298. 2017
#' @family classification
#' @export
privateEC <- function(data.sets=NULL,
                      label="phenos",
                      is.simulated=TRUE,
                      bias=0.4,
                      sim.type="mainEffect",
                      myrun="001",
                      update.freq=50,
                      corelearn.estimator="ReliefFbestK",
                      start.temp=0.1,
                      final.temp=10 ^ (-5),
                      tau.param=100,
                      threshold=4 / sqrt(nrow(data.sets$holdout)),
                      tolerance=1 / sqrt(nrow(data.sets$holdout)),
                      signal.names=NULL,
                      shortname="paramstring",
                      save.file=NULL,
                      verbose=FALSE) {
  if(is.null(data.sets) | (length(data.sets) < 2)) {
    stop("At least two data sets provided as first argument")
  }
  if(is.simulated & is.null(signal.names)) {
    warning("No signal names provided")
  }
  ptm <- proc.time()
  if(verbose) cat("running Relief-F importance for training and holdout sets\n")
  important.scores <- getImportanceScores(train.set=data.sets$train,
                                          holdo.set=data.sets$holdout,
                                          label=label,
                                          imp.estimator=corelearn.estimator,
                                          verbose=verbose)
  if(verbose) cat("private EC importance:", corelearn.estimator,
                  "elapsed time:", (proc.time() - ptm)[3], "\n")

  n <- nrow(data.sets[[1]])
  d <- ncol(data.sets[[1]]) - 1

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
      new.X_train <- data.sets$train[, c(kept.atts, label), drop=F]
      new.X_holdo <- data.sets$holdout[, c(kept.atts, label), drop=F]
      new.X_validation <- data.sets$validation[, c(kept.atts, label), drop=F]
      new.scores <- getImportanceScores(new.X_train, new.X_holdo, label=label)
      q1.scores <- new.scores[[1]]
      q2.scores <- new.scores[[2]]
      diff.scores <- abs(q1.scores - q2.scores)
      delta.q <- max(diff.scores)
      if(verbose) cat("\trunning randomForest\n")
      result.rf <- randomForest::randomForest(label ~ ., data=new.X_train)
      ftrain <- 1 - mean(result.rf$confusion[,"class.error"])
      if(verbose) cat("\tpredict\n")
      holdo.pred <- predict(result.rf, newdata=new.X_holdo)
      fholdo <- mean(holdo.pred == data.sets$holdout[, label])
      if(is.simulated) {
        validation.pred <- predict(result.rf, newdata=new.X_validation)
        fvalidation <- mean(validation.pred == data.sets$validation[, label])
      } else {
        fvalidation <- 0
      }
      if(abs(ftrain - fholdo) < (threshold + rnorm(1, 0, tolerance))) {
        fholdo <- ftrain
      } else {
        if(verbose) cat("\tadjust holdout with rnorm\n")
        fholdo <- fholdo + rnorm(1, 0, tolerance)
      }

      ftrains <- c(ftrains, ftrain)
      fholds <- c(fholds, fholdo)
      fvalidations <- c(fvalidations, fvalidation)

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
                       fvalidation=fvalidations,
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

  if(verbose) cat("privateEC elapsed time:", (proc.time() - ptm)[3], "\n")

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
#' @param data.sets A list of train, holdout and validation data frames
#' @param label A character vector of the class variable column name
#' @param is.simulated Is the data simulated (or real?)
#' @param sim.type A character vector of the type of simulation:
#' mainEffect/interactionErdos/interactionScalefree
#' @param myrun A character vector of a unique run identifier
#' @param update.freq A integer for the number of steps before update
#' @param pec.file A character vector filename of privateEC results
#' @param threshold A numeric, default 4 / sqrt(n) suggested in the
#'  thresholdout’s supplementary material (Dwork, et al.,2015)
#' @param tolerance A numeric, default 1 / sqrt(n) suggested in the
#'  thresholdout’s supplementary material (Dwork, et al.,2015)
#' @param signal.names A character vector of signal names in simulated data
#' @param shortname A character vector of a parameters separated by '_'
#' @param save.file A character vector for results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list containing:
#' \describe{
#'   \item{plots.data}{data frame of results, a row for each update}
#'   \item{melted.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @examples
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' nbias <- pct.signals * num.variables
#' signals <- sprintf("gene%04d", 1:nbias)
#' temp.pec.file <- tempfile(pattern="pEc_temp", tmpdir=tempdir())
#'
#' data.sets <- createSimulation(d=num.variables,
#'                               n=num.samples,
#'                               sim.type="mainEffect",
#'                               verbose=FALSE)
#' pec.results <- privateEC(data.sets=data.sets,
#'                          is.simulated=TRUE,
#'                          n=num.samples,
#'                          d=num.variables,
#'                          signal.names=signals,
#'                          save.file=temp.pec.file,
#'                          verbose=FALSE)
#' por.results <- originalThresholout(data.sets=data.sets,
#'                                is.simulated=TRUE,
#'                                n=num.samples,
#'                                d=num.variables,
#'                                signal.names=signals,
#'                                pec.file=temp.pec.file,
#'                                verbose=FALSE)
#' file.remove(temp.pec.file)
#' @family classification
#' @export
originalThresholout <- function(data.sets=NULL,
                                label="phenos",
                                is.simulated=TRUE,
                                sim.type="mainEffect",
                                myrun="001",
                                update.freq=50,
                                pec.file=NULL,
                                threshold=4 / sqrt(nrow(data.sets$holdout)),
                                tolerance=1 / sqrt(nrow(data.sets$holdout)),
                                signal.names=NULL,
                                shortname="paramstring",
                                save.file=NULL,
                                verbose=FALSE) {
  if(is.null(data.sets)) {
    stop("No data sets provided as data.sets argument")
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

  n <- nrow(data.sets[[1]])
  d <- ncol(data.sets[[1]]) - 1

  # compare with original privacy::
  if(verbose) cat("loading resuls from privacy EC", pec.file, "\n")
  load(pec.file)

  predictors.train <- as.matrix(data.sets$train[, 1:d])
  predictors.holdo <- as.matrix(data.sets$holdout[, 1:d])
  if(is.simulated) {
    predictors.validation <- as.matrix(data.sets$validation[, 1:d])
  }
  train.pheno <- data.sets$train[, label]
  holdo.pheno <- data.sets$holdout[, label]
  if(is.simulated) {
    validation.pheno <- data.sets$validation[, label]
    validation.pheno <- as.numeric(levels(validation.pheno))[validation.pheno]
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
      fvalidation <- mean((sign(predictors.validation %*% weights)) == validation.pheno)
    } else {
      fvalidation <- 0
    }
    if(k == 0) {
      ftrain <- 0.5
      fholdo <- 0.5
      fvalidation <- 0.5
    }
    noisy_vals[i,] <- c(ftrain, fholdo, fvalidation)
    var.names[[i]] <- colnames(data.sets$train)[topk]
    if(is.simulated) {
      correct.detect.ori <- c(correct.detect.ori,
                            sum(var.names[[i]] %in% signal.names))
    }
  }

  colnames(noisy_vals) <- c("ftrain", "fholdo", "fvalidation")
  pplots <- data.frame(num.atts, noisy_vals, alg=2)
  melted.ps <- reshape2::melt(pplots, id=c("num.atts", "alg"))

  if(!is.null(save.file)) {
    if(verbose) cat("saving to", save.file, "\n")
    save(pplots, melted.ps, correct.detect.ori, file=save.file)
  }

  if(verbose) cat("originalThresholout elapsed time:", (proc.time() - ptm)[3], "\n")

  list(plots.data=pplots,
       melted.data=melted.ps,
       correct=correct.detect.ori,
       elasped=(proc.time() - ptm)[3])
}

#' Private random forests algorithm
#'
#' Random Forest Thresholdout, which is TO with the feature selection
#' and classifier replaced with Random Forest.
#'
#' @param data.sets A list of train, holdout and validation data frames
#' @param label A character vector of the class variable column name
#' @param is.simulated Is the data simulated (or real?)
#' @param sim.type A character vector of the type of simulation:
#' mainEffect/interactionErdos/interactionScalefree
#' @param rf.importance.measure A character vector for the random forest importance measure
#' @param pec.file A character vector filename of privateEC results
#' @param update.freq A integer for the number of steps before update
#' @param threshold A numeric, default 4 / sqrt(n) suggested in the
#'  thresholdout’s supplementary material (Dwork, et al.,2015)
#' @param tolerance A numeric, default 1 / sqrt(n) suggested in the
#'  thresholdout’s supplementary material (Dwork, et al.,2015)
#' @param signal.names A character vector of signal names in simulated data
#' @param shortname A character vector of a parameters separated by '_'
#' @param save.file A character vector for results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list containing:
#' \describe{
#'   \item{plots.data}{data frame of results, a row for each update}
#'   \item{melted.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @examples
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' nbias <- pct.signals * num.variables
#' signals <- sprintf("gene%04d", 1:nbias)
#' temp.pec.file <- tempfile(pattern="pEc_temp", tmpdir=tempdir())
#'
#' data.sets <- createSimulation(d=num.variables,
#'                               n=num.samples,
#'                               sim.type="mainEffect",
#'                               verbose=FALSE)
#' pec.results <- privateEC(data.sets=data.sets,
#'                          is.simulated=TRUE,
#'                          n=num.samples,
#'                          d=num.variables,
#'                          signal.names=signals,
#'                          save.file=temp.pec.file,
#'                          verbose=FALSE)
#' prf.results <- privateRF(data.sets=data.sets,
#'                          is.simulated=TRUE,
#'                          n=num.samples,
#'                          d=num.variables,
#'                          signal.names=signals,
#'                          pec.file=temp.pec.file,
#'                          verbose=FALSE)
#' file.remove(temp.pec.file)
#' @export
privateRF <- function(data.sets=NULL,
                      label="phenos",
                      is.simulated=TRUE,
                      sim.type="mainEffect",
                      rf.importance.measure="MeanDecreaseGini",
                      pec.file=NULL,
                      update.freq=50,
                      threshold=4 / sqrt(nrow(data.sets$holdout)),
                      tolerance=1 / sqrt(nrow(data.sets$holdout)),
                      signal.names=NULL,
                      shortname="paramstring",
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

  n <- nrow(data.sets[[1]])
  d <- ncol(data.sets[[1]]) - 1

  if(verbose) cat("loading resuls from privacy EC", pec.file, "\n")
  load(pec.file)

  predictors.train <- as.matrix(data.sets$train[, 1:d])
  predictors.holdo <- as.matrix(data.sets$holdout[, 1:d])
  if(is.simulated) {
    predictors.validation <- as.matrix(data.sets$validation[, 1:d])
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
      fvalidation <- 0.5
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
        yvalidation.pred <- predict(rf.model.k, newdata=data.sets$validation)
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
      #   fvalidation <- sum((sign(predictors.validation %*% weights)) == pheno.validation)/n
      if(is.simulated) {
        fvalidation <- mean(yvalidation.pred == data.sets$validation$pheno)
      } else {
        fvalidation <- 0
      }
    }
    noisy_vals[i,] <- c(ftrain, fholdo, fvalidation)
    if(is.simulated) {
      correct.detect.rf <- c(correct.detect.rf,
                           sum(var.names[[i]] %in% signal.names))
    }
  }
  if(verbose) cat(i, "\n")

  colnames(noisy_vals) <- c("train.acc", "holdo.acc", "validation.acc")
  rfplots <- data.frame(num.atts, noisy_vals, alg=3)

  melted.rfs <- reshape2::melt(rfplots, id=c("num.atts", "alg"))

  if(!is.null(save.file)) {
    if(verbose) cat("saving plot data", save.file, "\n")
    save(rfplots, melted.rfs, correct.detect.rf, file=save.file)
  }

  if(verbose) cat("privateRF elasped time:", (proc.time() - ptm)[3], "\n")

  list(plots.data=rfplots,
       melted.data=melted.rfs,
       correct=correct.detect.rf,
       elasped=(proc.time() - ptm)[3])
}

#' Standard random forests algorithm serves as a baseline model
#'
#' @param data.sets A list of train, holdout and validation data frames
#' @param label A character vector of the class variable column name
#' @param is.simulated Is the data simulated (or real?)
#' @param sim.type A character vector of the type of simulation:
#' mainEffect/interactionErdos/interactionScalefree
#' @param signal.names A character vector of signal names in simulated data
#' @param shortname A character vector of a parameters separated by '_'
#' @param save.file A character vector for results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list containing:
#' \describe{
#'   \item{plots.data}{data frame of results, a row for each update}
#'   \item{melted.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @examples
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' nbias <- pct.signals * num.variables
#' signals <- sprintf("gene%04d", 1:nbias)
#' data.sets <- createSimulation(d=num.variables, n=num.samples,
#'                               sim.type="mainEffect", verbose=FALSE)
#' rra.results <- standardRF(data.sets=data.sets,
#'                           is.simulated=TRUE,
#'                           type=type,
#'                           verbose=FALSE,
#'                           signal.names=signals)
#' @family classification
#' @export
standardRF <- function(data.sets=NULL,
                       label="phenos",
                       is.simulated=TRUE,
                       sim.type="mainEffect",
                       signal.names=NULL,
                       shortname="paramstring",
                       save.file=NULL,
                       verbose=FALSE) {
  if(is.null(data.sets)) {
    stop("regularRF: No data sets provided as first argument")
  }
  if(is.simulated & is.null(signal.names)) {
    stop("regularRF: No signal names provided")
  }
  ptm <- proc.time()
  bag.simul <- randomForest::randomForest(label ~.,
                                          data=rbind(data.sets$train,
                                                     data.sets$holdout),
                                          importance=T)
  # rf.holdo.accu <- 1- bag.simul$prediction.error
  rf.holdo.accu <- 1 - mean(bag.simul$confusion[, "class.error"])
  if(is.simulated) {
    rf.pred <- predict(bag.simul, newdata=data.sets$validation)
    rf.validation.accu <- mean(rf.pred == data.sets$validation[, label])
  } else {
    rf.validation.accu  <- 0
  }
  if(verbose) cat("accuracies", rf.holdo.accu, rf.validation.accu, "\n")
  if(!is.null(save.file)) {
    save(rf.validation.accu, rf.holdo.accu, file=save.file)
  }

  rRaplots <- data.frame(num.atts=ncol(data.sets$train),
                         ftrain=1,
                         fholdo=rf.holdo.accu,
                         fvalidation=rf.validation.accu,
                         alg=4)
  melted.rra <- reshape2::melt(rRaplots, id=c("num.atts", "alg"))

  if(verbose) cat("regularRF elapsed time:", (proc.time() - ptm)[3], "\n")

  list(plots.data=rRaplots,
       melted.data=melted.rra,
       correct=ifelse(is.simulated, ncol(data.sets$train), 0),
       elasped=(proc.time() - ptm)[3])
}

#' Call C++ inbix Evaporative Cooling Privacy
#'
#' Assumes the inbix executable is in the PATH.
#'
#' @param data.sets A list of train, holdout and validation data frames
#' @param label A character vector of the class variable column name
#' @param is.simulated Is the data simulated (or real?)
#' @param bias A numeric for effect size in simulated signal variables
#' @param sim.type A character vector of the type of simulation:
#' mainEffect/interactionErdos/interactionScalefree
#' @param myrun A character vector of a unique run identifier
#' @param update.freq A integer for the number of steps before update
#' @param start.temp A numeric for EC starting temperature
#' @param final.temp A numeric for EC final temperature
#' @param tau.param A numeric for tau to control reduction schedule
#' @param shortname A character vector of a parameters separated by '_'
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
#' @examples
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' nbias <- pct.signals * num.variables
#' signals <- sprintf("gene%04d", 1:nbias)
#' sim.data <- createSimulation(d=num.variables, n=num.samples,
#'                              sim.type="mainEffect", verbose=FALSE)
#' pec.results <- privateECinbix(data.sets=sim.data,
#'                               n=num.samples,
#'                               verbose=FALSE)
#' @family classification
#' @export
privateECinbix <- function(data.sets=NULL,
                           label="phenos",
                           is.simulated=TRUE,
                           bias=0.4,
                           sim.type="mainEffect",
                           myrun="001",
                           update.freq=50,
                           start.temp=0.1,
                           final.temp=10 ^ (-5),
                           tau.param=100,
                           shortname="paramstring",
                           save.file=NULL,
                           verbose=FALSE) {
  if(is.null(data.sets)) {
    stop("privateECinbix: No data sets provided as first argument")
  }
  # check for inbix in the PATH or stop with an error
  if(Sys.which("inbix") == "") {
    stop("inbix is not installed or not in the PATH")
  }
  ptm <- proc.time()
  base.sim.prefix <- paste(sim.type, "_", shortname, sep="")
  unique.sim.prefix <- paste(sim.type, "_", shortname, sep="")
  correct.detect.inbix <- vector(mode="numeric")
  # write simple tab-separated vales (tsv) files
  if(verbose) cat("Writing ", unique.sim.prefix, ".sim.(train|holdout|validation).tab files\n")
  write.table(data.sets$train, paste(unique.sim.prefix, ".sim.train.tab", sep=""), sep="\t",
              row.names=FALSE, col.names=TRUE, quote=FALSE)
  write.table(data.sets$holdout, paste(unique.sim.prefix, ".sim.holdo.tab", sep=""), sep="\t",
              row.names=FALSE, col.names=TRUE, quote=FALSE)
  write.table(data.sets$validation, paste(unique.sim.prefix, ".sim.test.tab", sep=""), sep="\t",
              row.names=FALSE, col.names=TRUE, quote=FALSE)
  if(verbose) cat("Running inbix privacy EC on saved simulation data sets\n")
  system.command <- paste("OMP_NUM_THREADS=4 inbix ",
                          "--ec-privacy ",
                          "--depvarname ", label, " ",
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
  # if(verbose) cat("Reading inbix privacy EC attributes used results\n")
  # attrs.file <- paste(unique.sim.prefix, ".privateec.selected.attributes.tab", sep="")
  # if(verbose) cat(attrs.file, "\n")
  # attrs.table <- read.table(attrs.file, sep="\t", header=FALSE, stringsAsFactors=FALSE)

  if(verbose) cat("Reading inbix privacy EC algorithm run details\n")
  iters.file <- paste(unique.sim.prefix, ".privateec.privacy.iterations.tab", sep="")
  if(verbose) cat(iters.file, "\n")
  iters.table <- read.table(iters.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  if(is.simulated) {
    correct.detect.inbix <- as.integer(iters.table$Correct)
  }
  fxplots <- data.frame(num.atts=as.integer(iters.table$Keep),
                        ftrain=iters.table$TrainAcc,
                        fholdo=iters.table$HoldoutAcc,
                        fvalidation=iters.table$TestAcc,
                        alg=5)
  melted.fx <- reshape2::melt(fxplots, id=c("num.atts", "alg"))
  if(verbose) cat("Cleaning up inbix private EC results files\n")
  inbix.temp.files <- c(paste(unique.sim.prefix, ".sim.train.tab", sep=""),
                        paste(unique.sim.prefix, ".sim.holdo.tab", sep=""),
                        paste(unique.sim.prefix, ".sim.test.tab", sep=""),
                        paste(unique.sim.prefix, ".privateec.privacy.iterations.tab", sep=""),
                        paste(unique.sim.prefix, ".privateec.log", sep=""))
  #file.remove(c(inbix.conv.files, rf.files, attrs.file, iters.file))
  file.remove(inbix.temp.files)
  if(!is.null(save.file)) {
    if(verbose) {
      cat("saving results to ", save.file, "\n")
    }
    save(fxplots, melted.fx, correct.detect.inbix, bias, shortname,
         type, data.sets$train, data.sets$holdout, data.sets$validation,
         file=save.file)
  }
  if(verbose) cat("privacyECinbix elapsed time:", (proc.time() - ptm)[3], "\n")

  list(plots.data=fxplots,
       melted.data=melted.fx,
       correct=correct.detect.inbix,
       elasped=(proc.time() - ptm)[3])
}
