# classification.R - Trang Le and Bill White - Fall 2016/Spring 2017
#
# Classification algorithms used in the Bioinformatics paper:
# Differential privacy-based Evaporative Cooling feature selection and
# classification with Relief-F and Random Forests
# https://doi.org/10.1093/bioinformatics/btx298

#' Compute and return epistasis rank scores
#'
#' Modified SNPrank function
#'
#' @param G Am adjacency matrix
#' @param Gamma_vec A vector of prior knowledge
#' @family classification
#' @export
epistasisRank <- function(G, Gamma_vec) {
  n <- nrow(G)
  geneNames <- colnames(G)
  Gdiag <- diag(G)
  Gtrace <- sum(Gdiag)
  #colsum <- colSums(G)
  diag(G) <- 0
  Gtrace <- Gtrace * n
  colsumG <- colSums(G)
  #rowSumG <- rowSums(G)
  rowsum_denom <- matrix(0, n, 1)
  for (i in 1:n) {
    localSum <- 0
    for (j in 1:n) {
      factor <- ifelse(G[i, j] != 0, 1, 0)
      localSum <- localSum + (factor * colsumG[j])
    }
    rowsum_denom[i] <- localSum
  }
  #   gamma_vec <- rep(gamma, n)
  gamma_vec <- Gamma_vec
  gamma_matrix <- matrix(nrow = n, ncol = n, data = rep(gamma_vec, n))
  if (Gtrace) {
    b <- ((1 - gamma_vec)/n) + (Gdiag/Gtrace)
  }
  else {
    b <- ((1 - gamma_vec)/n)
  }
  D <- matrix(nrow = n, ncol = n, data = c(0))
  diag(D) <- 1/colsumG
  I <- diag(n)
  temp <- I - gamma_matrix * G %*% D
  r <- solve(temp, b)
  snpranks <- r/sum(r)
  saveTable <- data.frame(gene = geneNames, snprank = snpranks)
  sortedTable <- saveTable[order(saveTable$snprank, decreasing = TRUE), ]
  sortedTable
}

#' Compute and return importance scores (Relief-F scores)
#'
#' @param train.set A training data frame with last column as class
#' @param holdout.set A holdout data frame with last column as class
#' @param label A character vector of the class variable column name
#' @param importance.options A list importance operation parameters
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list with two data frames representing the importance scores
#' (Relief-F scores) for the train and holdout data sets.
#' @family classification
#' @export
getImportanceScores <- function(train.set=NULL,
                                holdout.set=NULL,
                                label="phenos",
                                importance.options=list(name = "relieff",
                                                        corelearn.estimator = "ReliefFbestK",
                                                        xgb.train.model = NULL,
                                                        xgb.holdout.model = NULL,
                                                        feature.names = NULL,
                                                        train.regain = NULL,
                                                        holdout.regain = NULL,
                                                        priorknowledge = NULL),
                                verbose=FALSE) {
  if (is.null(train.set)) {
    stop("getImportanceScores: No training data set provided")
  }
  if (is.null(holdout.set)) {
    stop("getImportanceScores: No holdout data set provided")
  }
  if (is.null(importance.options$feature.names)) {
    stop("getImportanceScores: No feature names provided")
  }
  if (verbose) cat("\tComputing importance scores\n")
  good.results <- FALSE
  if (importance.options$name == "relieff") {
    if (verbose) cat("\tRelief-F train\n")
    train.importance <- CORElearn::attrEval(label,
                                            data = train.set,
                                            estimator = importance.options$corelearn.estimator)
    if (verbose) cat("\tRelief-F holdout\n")
    holdout.importance <- CORElearn::attrEval(label,
                                              data = holdout.set,
                                              estimator = importance.options$corelearn.estimator)
    good.results <- TRUE
  }
  if (importance.options$name == "xgboost") {
    if (verbose) cat("\txgboost train\n")
    if (is.null(importance.options$xgb.train.model)) {
      stop("xgboost feature ranker [ ", importance.options$name, " ] requires a training model")
    }
    xgbResults <- xgboostRF(train.ds = train.set,
                            holdout.ds = holdout.set,
                            validation.ds = NULL,
                            label = label)
    train.importance <- xgboost::xgb.importance(feature_names = importance.options$feature.names,
                                                model = importance.options$xgb.train.model)
    if (verbose) cat("\txgboost holdout\n")
    if (is.null(importance.options$xgb.holdout.model)) {
      stop("xgboost feature ranker [ ", importance.options$name, " ] requires a holdout model")
    }
    holdout.importance <- xgboost::xgb.importance(feature_names = importance.options$feature.names,
                                                  model = xgbResults)
    good.results <- TRUE
  }
  if (importance.options$name == "epistasisrank") {
    if (verbose) cat("\tEpistasisRank train\n")
    train.importance <- epistasisRank(importance.options$train.regain, importance.options$priorknowledge)
    if (verbose) cat("\tEpistasisRank holdout\n")
    holdout.importance <- epistasisRank(importance.options$holdout.regain, importance.options$priorknowledge)
  }
  if (!good.results) {
    stop("Invalid feature importance ranker name [ ", importance.options$name, " ]")
  }

  list(data.frame(train.importance), data.frame(holdout.importance))
}

#' Private Evaporative Cooling feature selection and classification
#'
#' @param train.ds A data frame with training data and class labels
#' @param holdout.ds A data frame with holdout data and class labels
#' @param validation.ds A data frame with validation data and class labels
#' @param label A character vector of the class variable column name
#' @param learner.options A list of learner options
#' @param importance.options A list of importance options
#' @param is.simulated Is the data simulated (or real?)
#' @param bias A numeric for effect size in simulated signal variables
#' @param update.freq An integer the number of steps before update
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
#' sim.data <- createSimulation(num.samples = num.samples,
#'                              num.variables = num.variables,
#'                              pct.signals = pct.signals,
#'                              pct.train = 1 / 3,
#'                              pct.holdout = 1 / 3,
#'                              pct.validation = 1 /3,
#'                              sim.type = "mainEffect",
#'                              verbose = FALSE)
#' pec.results <- privateEC(train.ds = sim.data$train,
#'                          holdout.ds = sim.data$holdout,
#'                          validation.ds = sim.data$validation,
#'                          label = sim.data$class.label,
#'                          learner.options = list(name = "randomforest",
#'                                                 rf.ntree = 500,
#'                                                 rf.mtry = NULL),
#'                          importance.options = list(name = "relieff",
#'                                                    feature.names = colnames(sim.data$train),
#'                                                    corelearn.estimator = "ReliefFbestK"),
#'                          is.simulated = TRUE,
#'                          signal.names = sim.data$signal.names,
#'                          verbose = FALSE)
#' pec.results <- privateEC(train.ds = sim.data$train,
#'                          holdout.ds = sim.data$holdout,
#'                          validation.ds = sim.data$validation,
#'                          label = sim.data$class.label,
#'                          learner.options = list(name = "xgboost",
#'                                                 rf.ntree = 500,
#'                                                 rf.mtry = NULL),
#'                          importance.options = list(name = "relieff",
#'                                                    feature.names = colnames(sim.data$train),
#'                                                    corelearn.estimator = "ReliefFbestK"),
#'                          is.simulated = TRUE,
#'                          signal.names = sim.data$signal.names,
#'                          verbose = FALSE)
#' @note
#' Within thresholdout, we choose a threshold of 4 / sqrt(n) and
#' tolerance of 1 / sqrt(n) as suggested in the thresholdout’s supplementary
#' material (Dwork, et al., 2015).
#' @references
#' Trang Le, W. K. Simmons, M. Misaki, B.C. White, J. Savitz, J. Bodurka,
#' and B. A. McKinney. “Differential privacy-based Evaporative Cooling feature selection
#' and classification with Relief-F and Random Forests,”
#' Bioinformatics. Accepted. \url{https://doi.org/10.1093/bioinformatics/btx298}. 2017
#'
#' For more information see:
#' \href{http://insilico.utulsa.edu/index.php/privateec/}{Insilico Lab privateEC Page}
#' @family classification
#' @export
privateEC <- function(train.ds=NULL,
                      holdout.ds=NULL,
                      validation.ds=NULL,
                      label="phenos",
                      is.simulated=TRUE,
                      bias=0.4,
                      update.freq=50,
                      learner.options=list(name = "randomforest",
                                           rf.ntree = 500,
                                           rf.mtry = NULL),
                      importance.options=list(name = "relieff",
                                              feature.names = colnames(train.ds),
                                              corelearn.estimator = "ReliefFbestK"),
                      start.temp=0.1,
                      final.temp=10 ^ (-5),
                      tau.param=100,
                      threshold=4 / sqrt(nrow(train.ds)),
                      tolerance=1 / sqrt(nrow(train.ds)),
                      signal.names=NULL,
                      save.file=NULL,
                      verbose=FALSE) {
  if (is.null(train.ds) | is.null(holdout.ds)) {
    stop("At least train and holdout data sets must be provided")
  }
  if (is.simulated & is.null(signal.names)) {
    warning("For correct detection of signals, 'signal.names' parameter is required")
  }
  n <- nrow(train.ds)
  d <- ncol(train.ds) - 1
  param.mtry <- floor(sqrt(d))
  if (!is.null(learner.options$rf.mtry)) {
    param.mtry <- learner.options$rf.mtry
  }
  if ((param.mtry < 1) | (param.mtry > d)) {
      stop(paste("mtry parameter", param.mtry, "out of range 1 < mtry < d"))
  }
  ptm <- proc.time()
  if (verbose) cat("running Relief-F importance for training and holdout sets\n")
  if (importance.options$name == "relieff" |
      importance.options$name == "xgboost" |
      importance.options$name == "epistasisrank") {
    important.scores <- getImportanceScores(train.set = train.ds,
                                            holdout.set = holdout.ds,
                                            label = label,
                                            importance.options = importance.options,
                                            verbose = verbose)
  } else {
    stop("Importance method [", importance.options$name, "] is not recognized")
  }
  if (verbose) cat("private EC importance:", importance.options$name,
                   "elapsed time:", (proc.time() - ptm)[3], "\n")

  q1.scores <- important.scores[[1]]
  q2.scores <- important.scores[[2]]
  att.names <- rownames(q1.scores)
  diff.scores <- abs(q1.scores - q2.scores)
  delta.q <- max(diff.scores)
  # NOTE: bcw, below not used?
  # q1.scores.plot <- q1.scores
  fholds <- 0.5
  fvalidations <- 0.5
  ftrains <- 0.5
  correct.detect.ec <- vector(mode = "numeric")
  # NOTE: bcw, below not used?
  # oldAccuracy <- 0.5
  cur.vars.remain <- length(att.names)
  vars.remain <- c(0, cur.vars.remain)
  kept.atts <- att.names
  var.names <- list()

  if (verbose) cat("private EC optimization loop\n")
  T0 <- start.temp
  Tmin <- final.temp
  tau <- tau.param
  i <- 1
  myT <- T0
  num.updates <- 0
  while ((myT > Tmin) && (utils::tail(vars.remain, 1) > 2) &&
         (utils::tail(vars.remain, 2)[1] != utils::tail(vars.remain, 2)[2])) {
    diff <- diff.scores * (diff.scores > 10^(-3)) + 10^(-3) * (diff.scores < 10^(-3))
    PAs <- exp(-q1.scores / (2 * diff * myT))
    PAs <- PAs[kept.atts, ]
    # NOTE: bcw, sumPAs not used?
    sumPAs <- sum(PAs)
    scaled.PAs <- PAs / sumPAs
    cum.scaled.PAs <- cumsum(scaled.PAs)
    num.remv <- 1 # only remove 1 attribute
    prob.rands <- sort(stats::runif(num.remv, min = 0, max = 1))
    remv.atts <- kept.atts[prob.rands < cum.scaled.PAs][1]
    if (num.remv >= length(kept.atts)) {
      kept.atts <- NULL
      break
    } else {
      kept.atts <- setdiff(kept.atts, remv.atts)
    }
    # Compute train and holdout accuracy for new S_A attributes:
    att.names <- kept.atts
    if ((i %% update.freq) == 1) {
      num.updates <- num.updates + 1
      if (verbose) cat("step", i, "update", num.updates, "myT > Tmin",
                      myT, Tmin, "?\n")
      if (verbose) cat("\trecomputing scores with evaporated attributes removed\n")
      new.X_train <- train.ds[, c(kept.atts, label), drop = F]
      new.X_holdout <- holdout.ds[, c(kept.atts, label), drop = F]
      new.X_validation <- validation.ds[, c(kept.atts, label), drop = F]
      # new.scores <- getImportanceScores(new.X_train, new.X_holdout, label = label)
      if (importance.options$name == "relieff" |
          importance.options$name == "xgboost" |
          importance.options$name == "epistasisrank") {
        importance.options$feature.names <- colnames(new.X_train)
        new.scores <- getImportanceScores(train.set = new.X_train,
                                          holdout.set = new.X_holdout,
                                          label = label,
                                          importance.options = importance.options,
                                          verbose = verbose)
      } else {
        stop("Importance method [", importance.options$name, "] is not recognized")
      }

      q1.scores <- new.scores[[1]]
      q2.scores <- new.scores[[2]]
      diff.scores <- abs(q1.scores - q2.scores)
      delta.q <- max(diff.scores)
      if (param.mtry >= ncol(new.X_train) - 1) {
        kept.atts <- NULL
        break
      }
      method.valid <- FALSE
      if (learner.options$name == "randomforest") {
        if (verbose) cat("\trunning randomForest\n")
        model.formula <- stats::as.formula(paste(label, "~.", sep = ""))
        result.rf <- randomForest::randomForest(formula = model.formula,
                                                data = new.X_train,
                                                ntree = learner.options$rf.ntree,
                                                mtry = param.mtry)
        ftrain <- 1 - mean(result.rf$confusion[,"class.error"])
        if (verbose) cat("\tpredict\n")
        holdout.pred <- stats::predict(result.rf, newdata = new.X_holdout)
        fholdout <- mean(holdout.pred == holdout.ds[, label])
        if (is.simulated) {
          validation.pred <- stats::predict(result.rf, newdata = new.X_validation)
          fvalidation <- mean(validation.pred == validation.ds[, label])
        } else {
          fvalidation <- 0
        }
        method.valid <- TRUE
      }
      # bcw - 1/31/18 - xgboost package
      if (learner.options$name == "xgboost") {
        if (verbose) cat("\trunning xgboost\n")
        xgbResults <- xgboostRF(train.ds = new.X_train,
                                holdout.ds = new.X_holdout,
                                validation.ds = new.X_validation,
                                rf.ntree = learner.options$rf.ntree,
                                label = label)
        ftrain <- xgbResults$algo.acc$train.acc
        fholdout <- xgbResults$algo.acc$holdout.acc
        fvalidation <- xgbResults$algo.acc$validation.acc
        method.valid <- TRUE
      }
      if (!method.valid) {
        stop("Learner method [ ", learner.options$name, " ] is not recognized")
      }

      if (verbose) {
        cat("\tthreshold: ", threshold, "\n")
        cat("\ttolerance:", threshold, "\n")
        cat("\tftrain:    ", ftrain, "\n")
        cat("\tfholdout:  ", fholdout, "\n")
      }

      if (abs(ftrain - fholdout) < (threshold + stats::rnorm(1, 0, tolerance))) {
        if (verbose) {
          cat("\tClassification error difference is less than threshold + random tolerance\n")
          cat("\tAdjusting holdout error to training error\n")
        }
        fholdout <- ftrain
      } else {
        if (verbose) cat("\tadjust holdout by adding normal random value 0 to tolerance\n")
        fholdout <- fholdout + stats::rnorm(1, 0, tolerance)
      }
      ftrains <- c(ftrains, ftrain)
      fholds <- c(fholds, fholdout)
      fvalidations <- c(fvalidations, fvalidation)

      if (verbose) cat("\tadjusting temperature\n")
      myT <- myT * exp(-1 / tau) # dropping T
      if (verbose) cat("\t", i - 1, ': ', myT, '\n')

      if (verbose) cat("\tcollecting results\n")
      cur.vars.remain <- length(att.names)
      vars.remain <- c(vars.remain, cur.vars.remain)
      var.names[[num.updates]] <- kept.atts
      if (is.simulated) {
        correct.detect.ec <- c(correct.detect.ec,
                               sum(var.names[[num.updates]] %in% signal.names))
      }
    }
    i <- i + 1
  }
  elapsed.time <- (proc.time() - ptm)[3]
  if (verbose) cat("private EC optimization loop performed", num.updates,
                   "updates in", elapsed.time, " seconds\n")

  # prepare results for returning to caller
  vars.remain <- vars.remain[-1] # remove the first value 0
  fplots <- data.frame(vars.remain,
                       train.acc = ftrains,
                       holdout.acc = fholds,
                       validation.acc = fvalidations,
                       alg = 1)
  fplots <- fplots[-1, ] # remove the first row
  melted.fs <- reshape2::melt(fplots, id = c("vars.remain", "alg"))

  # save the results to an Rdata file if requested
  if (!is.null(save.file)) {
    if (verbose) {
      cat("saving results to ", save.file, "\n")
    }
    save(fplots, melted.fs, correct.detect.ec, elapsed.time, n, d,
         signal.names, threshold, tolerance, bias, file = save.file)
  }

  if (verbose) cat("privateEC elapsed time:", elapsed.time, "\n")

  list(algo.acc = fplots,
       ggplot.data = melted.fs,
       correct = correct.detect.ec,
       elasped = elapsed.time)
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
#' temp.pec.file <- tempfile(pattern = "pEc_temp", tmpdir = tempdir())
#' sim.data <- createSimulation(num.variables = num.variables,
#'                              num.samples = num.samples,
#'                              sim.type = "mainEffect",
#'                              pct.train = 1 / 3,
#'                              pct.holdout = 1 / 3,
#'                              pct.validation = 1 / 3,
#'                              verbose = FALSE)
#' pec.results <- privateEC(train.ds = sim.data$train,
#'                          holdout.ds = sim.data$holdout,
#'                          validation.ds = sim.data$validation,
#'                          label = sim.data$class.label,
#'                          is.simulated = TRUE,
#'                          signal.names = sim.data$signal.names,
#'                          save.file = temp.pec.file,
#'                          verbose = FALSE)
#' por.results <- originalThresholdout(train.ds = sim.data$train,
#'                                     holdout.ds = sim.data$holdout,
#'                                     validation.ds = sim.data$validation,
#'                                     label = sim.data$class.label,
#'                                     is.simulated = TRUE,
#'                                     signal.names = sim.data$signal.names,
#'                                     pec.file = temp.pec.file,
#'                                     verbose = FALSE)
#' file.remove(temp.pec.file)
#' @family classification
#' @export
originalThresholdout <- function(train.ds=NULL,
                                 holdout.ds=NULL,
                                 validation.ds=NULL,
                                 label="phenos",
                                 is.simulated=TRUE,
                                 update.freq=50,
                                 pec.file=NULL,
                                 threshold=4 / sqrt(nrow(holdout.ds)),
                                 tolerance=1 / sqrt(nrow(holdout.ds)),
                                 signal.names=NULL,
                                 save.file=NULL,
                                 verbose=FALSE) {
  if (is.null(train.ds) | is.null(holdout.ds)) {
    stop("At least train and holdout data sets must be provided")
  }
  if (is.simulated & is.null(signal.names)) {
    stop("No signal names provided")
  }
  if (is.null(pec.file)) {
    stop("No previous privateEC results file in pec.file argument")
  }
  if (!file.exists(pec.file)) {
    stop("privateEC results file expected in pec.file argument:", pec.file)
  }

  ptm <- proc.time()

  n <- nrow(train.ds)
  d <- ncol(train.ds) - 1

  # compare with original privacy::
  if (verbose) cat("loading resuls from privacy EC", pec.file, "\n")
  load(pec.file)

  predictors.train <- as.matrix(train.ds[, 1:d])
  predictors.holdout <- as.matrix(holdout.ds[, 1:d])
  if (is.simulated) {
    predictors.validation <- as.matrix(validation.ds[, 1:d])
  }
  train.pheno <- train.ds[, label]
  holdout.pheno <- holdout.ds[, label]
  if (is.simulated) {
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
  correct.detect.ori <- vector(mode = "numeric")
  var.names <- list()

  for (i in 1:numks) {
    k <- vars.remain[i]
    topk <- utils::tail(sortanswers, k)
    weights <- matrix(0, d, 1)
    weights[topk] <- sign(trainanswers[topk])
    ftrain <- mean((sign(predictors.train %*% weights)) == train.pheno)
    fholdout <- mean((sign(predictors.holdout %*% weights)) == holdout.pheno)
    if (abs(ftrain - fholdout) < threshold + stats::rnorm(1, 0, tolerance)) {
      fholdout <- ftrain
    } else {
      fholdout <- fholdout + stats::rnorm(1, 0, tolerance)
    }
    if (is.simulated) {
      fvalidation <- mean((sign(predictors.validation %*% weights)) == validation.pheno)
    } else {
      fvalidation <- 0
    }
    if (k == 0) {
      ftrain <- 0.5
      fholdout <- 0.5
      fvalidation <- 0.5
    }
    noisy_vals[i,] <- c(ftrain, fholdout, fvalidation)
    var.names[[i]] <- colnames(train.ds)[topk]
    if (is.simulated) {
      correct.detect.ori <- c(correct.detect.ori,
                              sum(var.names[[i]] %in% signal.names))
    }
  }

  colnames(noisy_vals) <- c("train.acc", "holdout.acc", "validation.acc")
  pplots <- data.frame(vars.remain, noisy_vals, alg = 2)
  melted.ps <- reshape2::melt(pplots, id = c("vars.remain", "alg"))

  if (!is.null(save.file)) {
    if (verbose) cat("saving to", save.file, "\n")
    save(pplots, melted.ps, correct.detect.ori, file = save.file)
  }

  if (verbose) cat("originalThresholout elapsed time:", (proc.time() - ptm)[3], "\n")

  list(algo.acc = pplots,
       ggplot.data = melted.ps,
       correct = correct.detect.ori,
       elasped = (proc.time() - ptm)[3])
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
#' temp.pec.file <- tempfile(pattern = "pEc_temp", tmpdir = tempdir())
#' sim.data <- createSimulation(num.variables = num.variables,
#'                              num.samples = num.samples,
#'                              sim.type = "mainEffect",
#'                              pct.train = 1 / 3,
#'                              pct.holdout = 1 / 3,
#'                              pct.validation = 1 / 3,
#'                              verbose = FALSE)
#' pec.results <- privateEC(train.ds = sim.data$train,
#'                          holdout.ds = sim.data$holdout,
#'                          validation.ds = sim.data$validation,
#'                          label = sim.data$class.label,
#'                          is.simulated = TRUE,
#'                          signal.names = sim.data$signal.names,
#'                          save.file = temp.pec.file,
#'                          verbose = FALSE)
#' prf.results <- privateRF(train.ds = sim.data$train,
#'                          holdout.ds = sim.data$holdout,
#'                          validation.ds = sim.data$validation,
#'                          label = sim.data$class.label,
#'                          is.simulated = TRUE,
#'                          signal.names = sim.data$signal.names,
#'                          pec.file = temp.pec.file,
#'                          verbose = FALSE)
# file.remove(temp.pec.file)
#' @family classification
#' @export
privateRF <- function(train.ds=NULL,
                      holdout.ds=NULL,
                      validation.ds=NULL,
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
  if (is.null(train.ds) | is.null(holdout.ds)) {
    stop("At least train and holdout data sets must be provided")
  }
  n <- nrow(train.ds)
  d <- ncol(train.ds) - 1
  param.mtry <- rf.mtry
  if (is.null(rf.mtry)) {
    param.mtry <- floor(sqrt(d))
  }
  if (param.mtry < 1 | param.mtry > d) {
    stop(paste("mtry parameter", param.mtry, "out of range 1 < mtry < d"))
  }
  if (is.simulated & is.null(signal.names)) {
    stop("privateRF: No signal names provided")
  }
  if (is.null(pec.file)) {
    stop("privateRF: No previous privateEC results file in pec.file argument")
  }
  if (!file.exists(pec.file)) {
    stop("privateRF: Previous privateEC results file in pec.file argument:", pec.file)
  }

  ptm <- proc.time()

  if (verbose) cat("loading resuls from privacy EC", pec.file, "\n")
  load(pec.file)

  predictors.train <- as.matrix(train.ds[, 1:d])
  predictors.holdout <- as.matrix(holdout.ds[, 1:d])
  # if (is.simulated) {
  #   if (!is.null(validation.ds)) {
  #     # NOTE: bcw, predictors.validation not used?
  #     predictors.validation <- as.matrix(validation.ds[, 1:d])
  #   }
  # }
  train.rf <- randomForest::randomForest(x = predictors.train,
                                         y = train.ds[, label],
                                         ntree = rf.ntree,
                                         mtry = param.mtry,
                                         importance = T)
  train.imp <- train.rf$importance[, rf.importance.measure]
  holdout.rf <- randomForest::randomForest(x = predictors.holdout,
                                           y = holdout.ds[, label],
                                           ntree = rf.ntree,
                                           mtry = param.mtry,
                                           importance = T)
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
  correct.detect.rf <- vector(mode = "numeric")
  var.names <- list()
  if (verbose) cat("loop for", numks, ": ")
  for (i in 1:numks) {
    if (verbose) cat(i, " ")
    k <- vars.remain[i]
    if (k == 0) {
      ftrain <- 0.5
      fholdout <- 0.5
      fvalidation <- 0.5
    } else {
      topk <- utils::tail(sortanswers, k)
      var.names[[i]] <- colnames(train.ds)[topk]
      # rf classifier:
      data.k <- train.ds[, topk, drop = F]
      if (param.mtry >= ncol(data.k)) {
        break
      }
      rf.model.k <- randomForest::randomForest(x = data.k,
                                               y = train.ds[, label],
                                               ntree = rf.ntree,
                                               mtry = param.mtry)
      # ytrain.pred <- stats::predict(rf.model.k, newdata = train.ds)
      yholdout.pred <- stats::predict(rf.model.k, newdata = holdout.ds)
      if (is.simulated) {
        if (!is.null(validation.ds)) {
          yvalidation.pred <- stats::predict(rf.model.k, newdata = validation.ds)
        }
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
      if (is.simulated) {
        if (!is.null(validation.ds)) {
          fvalidation <- mean(yvalidation.pred == validation.ds[, label])
        } else {
          fvalidation <- 0
        }
      } else {
        fvalidation <- 0
      }
    }
    noisy_vals[i,] <- c(ftrain, fholdout, fvalidation)
    if (is.simulated) {
      correct.detect.rf <- c(correct.detect.rf,
                           sum(var.names[[i]] %in% signal.names))
    }
  }
  if (verbose) cat(i, "\n")

  colnames(noisy_vals) <- c("train.acc", "holdout.acc", "validation.acc")
  rfplots <- data.frame(vars.remain, noisy_vals, alg = 3)

  melted.rfs <- reshape2::melt(rfplots, id = c("vars.remain", "alg"))

  if (!is.null(save.file)) {
    if (verbose) cat("saving plot data", save.file, "\n")
    save(rfplots, melted.rfs, correct.detect.rf, file = save.file)
  }

  if (verbose) cat("privateRF elasped time:", (proc.time() - ptm)[3], "\n")

  list(algo.acc = rfplots,
       ggplot.data = melted.rfs,
       correct = correct.detect.rf,
       elasped = (proc.time() - ptm)[3])
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
#' sim.data <- createSimulation(num.variables = num.variables,
#'                              num.samples = num.samples,
#'                              sim.type = "mainEffect",
#'                              pct.train = 1 / 3,
#'                              pct.holdout = 1 / 3,
#'                              pct.validation = 1 / 3,
#'                              verbose = FALSE)
#' rra.results <- standardRF(train.ds=sim.data$train,
#'                           holdout.ds=sim.data$holdout,
#'                           validation.ds=sim.data$validation,
#'                           label=sim.data$class.label,
#'                           is.simulated=TRUE,
#'                           verbose=FALSE,
#'                           signal.names=sim.data$signal.names)
#' @family classification
#' @export
standardRF <- function(train.ds=NULL,
                       holdout.ds=NULL,
                       validation.ds=NULL,
                       label="phenos",
                       rf.ntree=500,
                       rf.mtry=NULL,
                       is.simulated=TRUE,
                       signal.names=NULL,
                       save.file=NULL,
                       verbose=FALSE) {
  if (is.null(train.ds) | is.null(holdout.ds)) {
    stop("At least train and holdout data sets must be provided")
  }
  # NOTE: bcw, 'n' not used?
  # n <- nrow(train.ds)
  d <- ncol(train.ds) - 1
  param.mtry <- rf.mtry
  if (is.null(rf.mtry)) {
    param.mtry <- floor(sqrt(d))
  }
  if (param.mtry < 1 | param.mtry > d) {
    stop(paste("mtry parameter", param.mtry, "out of range 1 < mtry < d"))
  }
  if (is.simulated & is.null(signal.names)) {
    stop("regularRF: No signal names provided")
  }
  ptm <- proc.time()
  bag.simul <- randomForest::randomForest(stats::as.formula(paste(label, "~ .", sep = "")),
                                          data = rbind(train.ds, holdout.ds),
                                          ntree = rf.ntree,
                                          mtry = param.mtry,
                                          importance = T)
  # rf.holdo.accu <- 1- bag.simul$prediction.error
  rf.holdout.accu <- 1 - mean(bag.simul$confusion[, "class.error"])
  if (is.simulated) {
    rf.pred <- stats::predict(bag.simul, newdata = validation.ds)
    rf.validation.accu <- mean(rf.pred == validation.ds[, label])
  } else {
    rf.validation.accu  <- 0
  }
  if (verbose) cat("accuracies", rf.holdout.accu, rf.validation.accu, "\n")
  if (!is.null(save.file)) {
    save(rf.validation.accu, rf.holdout.accu, file = save.file)
  }

  rRaplots <- data.frame(vars.remain = ncol(train.ds),
                         train.acc = 1,
                         holdout.acc = rf.holdout.accu,
                         validation.acc = rf.validation.accu,
                         alg = 4)
  melted.rra <- reshape2::melt(rRaplots, id = c("vars.remain", "alg"))

  elapsed.time <- (proc.time() - ptm)[3]
  if (verbose) cat("regularRF elapsed time:", elapsed.time, "\n")

  list(algo.acc = rRaplots,
       ggplot.data = melted.rra,
       correct = ifelse(is.simulated, ncol(train.ds), 0),
       elasped = elapsed.time)
}

#' xgboost random forests algorithm
#'
#' Scalable and Flexible Gradient Boosting
#' XGBoost is short for “Extreme Gradient Boosting”, where the term “Gradient Boosting” is proposed in the paper
#' Greedy Function Approximation: A Gradient Boosting Machine, by Friedman. XGBoost is based on this original model.
#' This is a function using gradient boosted trees for privacyEC.
#' .
#' @param train.ds A data frame with training data and class labels
#' @param holdout.ds A data frame with holdout data and class labels
#' @param validation.ds A data frame with validation data and class labels
#' @param label A character vector of the class variable column name
#' @param rf.ntree An integer the number of trees in the random forest
#' @param num.iter An integer number of xgboost iterations
#' @param num.threads An integer for OpenMP number of cores
#' @param max.depth An integer aximum tree depth
#' @param learn.rate A numeric gradient learning rate 0-1
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
#' sim.data <- createSimulation(num.variables = num.variables,
#'                              num.samples = num.samples,
#'                              sim.type = "mainEffect",
#'                              pct.train = 1 / 3,
#'                              pct.holdout = 1 / 3,
#'                              pct.validation = 1 / 3,
#'                              verbose = FALSE)
#' rra.results <- xgboostRF(train.ds=sim.data$train,
#'                          holdout.ds=sim.data$holdout,
#'                          validation.ds=sim.data$validation,
#'                          label=sim.data$class.label,
#'                          num.iter = 3,
#'                          rf.ntree = 500,
#'                          max.depth = 10,
#'                          is.simulated=TRUE,
#'                          verbose=FALSE,
#'                          signal.names=sim.data$signal.names)
#' @family classification
#' @export
xgboostRF <- function(train.ds=NULL,
                      holdout.ds=NULL,
                      validation.ds=NULL,
                      label="phenos",
                      rf.ntree=500,
                      num.iter=1,
                      num.threads=1,
                      max.depth=4,
                      learn.rate=1,
                      save.file=NULL,
                      verbose=FALSE) {
  if (is.null(train.ds) | is.null(holdout.ds)) {
    stop("At least train and holdout data sets must be provided")
  }
  pheno.col <- which(colnames(train.ds) == label)
  if (length(pheno.col) > 1) {
    stop("More than one phenotype column with label [ ", label, " ]")
  }
  if ((pheno.col < 1) | (pheno.col > ncol(train.ds))) {
    stop("Could not find phenotype column with label [ ", label, " ]")
  }
  ptm <- proc.time()
  if (verbose) cat("Running xgboost\n")
  var.names <- colnames(train.ds[, -pheno.col])
  train.data <- as.matrix(train.ds[, -pheno.col])
  colnames(train.data) <- var.names
  holdout.data <- as.matrix(holdout.ds[, -pheno.col])
  colnames(holdout.data) <- var.names
  train.pheno <- as.integer(train.ds[, label]) - 1
  if (verbose) print(table(train.pheno))
  holdo.pheno <- as.integer(holdout.ds[, label]) - 1
  if (verbose) print(table(holdo.pheno))
  dtrain <- xgboost::xgb.DMatrix(data = train.data, label = train.pheno)
  dholdo <- xgboost::xgb.DMatrix(data = holdout.data, label = holdo.pheno)
  bst <- xgboost::xgboost(data = dtrain,
                          eta = learn.rate,
                          nthread = num.threads,
                          nrounds = num.iter,
                          max.depth = max.depth,
                          num_parallel_tree = rf.ntree,
                          objective = "binary:logistic")
  pred.prob <- predict(bst, dtrain)
  pred.bin <- as.numeric(pred.prob > 0.5)
  rf.train.accu <- 1 - mean(pred.bin != train.pheno)
  cat("training-accuracy:", rf.train.accu, "\n")
  if (verbose) cat("Predict using the best tree\n")
  pred.prob <- predict(bst, dholdo)
  # The most important thing to remember is that to do a classification, you just do a
  # regression to the label and then apply a threshold.
  pred.bin <- as.numeric(pred.prob > 0.5)
  if (verbose) print(table(pred.bin))
  # compute classification accuracy
  rf.holdo.accu <- 1 - mean(pred.bin != holdo.pheno)
  cat("holdout-accuracy:", rf.holdo.accu, "\n")
  if (verbose) cat("Computing variable importance\n")
  importance_matrix <- xgboost::xgb.importance(model = bst)
  if (verbose) print(importance_matrix)
  if (!is.null(validation.ds)) {
    validation.pheno <- as.integer(validation.ds[, label]) - 1
    validation.data <- as.matrix(validation.ds[, -pheno.col])
    colnames(validation.data) <- var.names
    dvalidation <- xgboost::xgb.DMatrix(data = validation.data, label = validation.pheno)
    rf.prob <- stats::predict(bst, newdata = dvalidation)
    rf.bin <- as.numeric(rf.prob > 0.5)
    if (verbose) print(table(rf.bin))
    rf.validation.accu <- 1 - mean(rf.bin != validation.pheno)
    cat("validation-accuracy:", rf.validation.accu, "\n")
  } else {
    rf.validation.accu  <- 0
  }
  if (verbose) cat("accuracies", rf.train.accu, rf.holdo.accu, rf.validation.accu, "\n")
  if (!is.null(save.file)) {
    save(rf.validation.accu, rf.holdout.accu, file = save.file)
  }
  xgb.plots <- data.frame(vars.remain = ncol(train.ds),
                          train.acc = rf.train.accu,
                          holdout.acc = rf.holdo.accu,
                          validation.acc = rf.validation.accu,
                          alg = 5)
  melted.rra <- reshape2::melt(xgb.plots, id = c("vars.remain", "alg"))
  elapsed.time <- (proc.time() - ptm)[3]
  if (verbose) cat("xgboostRF elapsed time:", elapsed.time, "\n")

  list(algo.acc = xgb.plots,
       ggplot.data = melted.rra,
       elasped = elapsed.time)
}
