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
  for (current.iteration in 1:n) {
    localSum <- 0
    for (j in 1:n) {
      factor <- ifelse(G[current.iteration, j] != 0, 1, 0)
      localSum <- localSum + (factor * colsumG[j])
    }
    rowsum_denom[current.iteration] <- localSum
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
  final.ranks <- r/sum(r)
  saveTable <- data.frame(gene = geneNames, rank = final.ranks)
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
  if (dim(train.set)[2] != dim(holdout.set)[2]) {
    print(dim(train.set))
    print(dim(holdout.set))
    stop("Training and holdout data sets have different sized columns")
  }
  if (verbose) cat("\tComputing importance scores\n")
  good.results <- FALSE
  if (importance.options$name == "relieff") {
    if (verbose) cat("\tRelief-F train\n")
    train.importance <- CORElearn::attrEval(label, data = train.set,
                                            estimator = importance.options$corelearn.estimator)
    if (verbose) cat("\tRelief-F holdout\n")
    holdout.importance <- CORElearn::attrEval(label, data = holdout.set,
                                              estimator = importance.options$corelearn.estimator)
    good.results <- TRUE
  }
  # ----------------------------------------------------------------
  if (importance.options$name == "epistasisrank") {
    if (verbose) cat("\tEpistasisRank train\n")
    train.imp <- epistasisRank(importance.options$train.regain, importance.options$priorknowledge)
    train.importance <- as.numeric(train.imp$rank)
    names(train.importance) <- importance.options$feature.names
    if (verbose) cat("\tEpistasisRank holdout\n")
    holdout.imp <- epistasisRank(importance.options$holdout.regain, importance.options$priorknowledge)
    holdout.importance <- as.numeric(holdout.imp)
    names(holdout.importance) <- importance.options$feature.names
    good.results <- TRUE
  }
  # ----------------------------------------------------------------
  if (!good.results) {
    stop("Variable importance ranker name [ ", importance.options$name, " ] not found or failed")
  }
  list(train.scores = train.importance, holdout.scores = holdout.importance)
}

#' Private Evaporative Cooling feature selection and classification
#'
#' @param train.ds A data frame with training data and class labels
#' @param holdout.ds A data frame with holdout data and class labels
#' @param validation.ds A data frame with validation data and class labels
#' @param label A character vector of the class variable column name
#' @param learner.options A list of learner parameters
#' @param importance.options A list of importance parameters
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
#' ntree <- 500
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
#'                          importance.options = list(name = "relieff",
#'                                                    feature.names = colnames(sim.data$train),
#'                                                    corelearn.estimator = "ReliefFbestK"),
#'                          learner.options = list(name = "randomforest",
#'                                                 rf.ntree = ntree),
#'                          is.simulated = TRUE,
#'                          signal.names = sim.data$signal.names,
#'                          verbose = FALSE)
#' pec.results <- privateEC(train.ds = sim.data$train,
#'                          holdout.ds = sim.data$holdout,
#'                          validation.ds = sim.data$validation,
#'                          label = sim.data$class.label,
#'                          importance.options = list(name = "relieff",
#'                                                    feature.names = colnames(sim.data$train),
#'                                                    corelearn.estimator = "ReliefFbestK"),
#'                          learner.options = list(name = "xgboost", rf.ntree = ntree),
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
                      update.freq=5,
                      importance.options=list(name = "relieff",
                                              feature.names = colnames(train.ds),
                                              corelearn.estimator = "ReliefFbestK"),
                      learner.options=list(name = "randomforest",
                                           rf.ntree = 500,
                                           rf.mtry = 5,
                                           num.rounds = 1,
                                           num.threads = 2,
                                           max.depth = 6,
                                           learn.rate = 1,
                                           obj.func = "binary:logistic"), # max_depth maximum depth of a tree. Default: 6
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
  if (ncol(train.ds) != ncol(holdout.ds)) {
    stop("Training and holdout data sets must have the same number of variables (columns)")
  }
  if (sum(colnames(train.ds) != colnames(holdout.ds)) > 0) {
    stop("Training and holdout data sets must have the same variable (column) names")
  }
  if (!is.null(validation.ds)) {
    if (ncol(train.ds) != ncol(validation.ds)) {
      stop("Training and validation data sets must have the same number of variables (columns)")
    }
    if (sum(colnames(train.ds) != colnames(validation.ds)) > 0) {
      stop("Training and validation data sets must have the same variable (column) names")
    }
  }
  if (is.simulated & is.null(signal.names)) {
    warning("For correct detection of signals, 'signal.names' parameter is required")
  }
  num.data.rows <- nrow(train.ds)
  num.data.cols <- ncol(train.ds) - 1
  orig.var.names <- colnames(train.ds)
  if (is.null(learner.options$rf.mtry)) {
    param.mtry <- floor(sqrt(num.data.cols))
    if ((param.mtry < 1) | (param.mtry > num.data.cols)) {
      stop(paste("mtry parameter", param.mtry, "out of range 1 < mtry < num.data.cols"))
    }
    learner.options$rf.mtry <- param.mtry
  }
  if (is.null(importance.options$feature.names)) {
    cat("WARNING: No features names in importance options list, using column names of training data set\n")
    importance.options$feature.names <- colnames(train.ds)
  }
  # ---------------------------------------------------------------------------
  ptm <- proc.time()
  cat("\nRunning the privateEC algorithm...\n")
  if (verbose) cat("running [", importance.options$name, "] importance for training and holdout sets\n")
  if (importance.options$name == "relieff" |
      importance.options$name == "epistasisrank") {
    # pass-through the importance function options list
    initial.scores <- getImportanceScores(train.set = train.ds,
                                          holdout.set = holdout.ds,
                                          label = label,
                                          importance.options = importance.options,
                                          verbose = verbose)
  } else {
    stop("Importance method [", importance.options$name, "] is not recognized")
  }
  if (length(initial.scores) != 2) {
    stop("Could not get importance scores")
  }
  if (verbose) cat("private EC importance:", importance.options$name,
                   "elapsed time:", (proc.time() - ptm)[3], "\n")
  q1.scores <- initial.scores[[1]]
  q2.scores <- initial.scores[[2]]
  if (length(q1.scores) != length(q2.scores)) {
    stop("Initial imoprtance scores for train and holdout data sets are not the same size")
  }
  orig.var.names <- names(q1.scores)
  diff.scores <- abs(q1.scores - q2.scores)
  delta.q <- max(diff.scores)
  all.train.acc <- vector(mode = "numeric")
  all.holdout.acc <- vector(mode = "numeric")
  all.validation.acc <- vector(mode = "numeric")
  vars.remain.per.update <- vector(mode = "numeric")
  correct.detect.ec <- vector(mode = "numeric")
  current.var.names <- orig.var.names
  prev.var.names <- current.var.names
  if (verbose) cat("private EC optimization loop\n")
  T0 <- start.temp
  Tmin <- final.temp
  tau <- tau.param
  current.temp <- T0
  current.iteration <- 1
  num.updates <- 0
  keep.looping <- ifelse(length(current.var.names) > 1, TRUE, FALSE)
  while (keep.looping & (current.temp > Tmin)) {
    diff <- diff.scores * (diff.scores > 10^(-3)) + 10^(-3) * (diff.scores < 10^(-3))
    PAs <- exp(-q1.scores / (2 * diff * current.temp))
    PAs <- PAs[current.var.names]
    sumPAs <- sum(PAs)
    scaled.PAs <- PAs / sumPAs
    cum.scaled.PAs <- cumsum(scaled.PAs)
    num.remv <- 1 # only remove 1 attribute
    prob.rands <- sort(stats::runif(num.remv, min = 0, max = 1))
    var.names.to.remove <- current.var.names[prob.rands < cum.scaled.PAs][1]
    if (length(var.names.to.remove) != num.remv) {
      cat("\tremove by probability distribution: cannot remove any variables")
      keep.looping <- FALSE
      next
    } else {
      current.var.names <- setdiff(current.var.names, var.names.to.remove)
      if (length(current.var.names) < 2) {
        cat("\tsetdiff remove: last variable removed\n")
        keep.looping <- FALSE
        next
      }
      if (verbose) cat("\t[", length(current.var.names), "] remaining variables\n")
    }
    if (length(prev.var.names) == length(current.var.names)) {
      cat("\tCONVERGENCE to a set of variables that are unchanging\n")
      keep.looping <- FALSE
      break
    }
    # run on the current set of variables - "update"
    if ((current.iteration == 1) | (current.iteration %% update.freq) == 0) {
      num.updates <- num.updates + 1
      if (verbose) cat("\n\niteration [", current.iteration, " ] UPDATE [", num.updates,
                       "] ", length(current.var.names), "current.temp > Tmin? [", current.temp, "] > [", Tmin, "]?\n")
      # re-compute imnportance
      if (verbose) cat("\trecomputing scores with evaporated attributes removed\n")
      new.train.ds <- train.ds[, c(current.var.names, label)]
      new.holdout.ds <- holdout.ds[, c(current.var.names, label)]
      new.validation.ds <- validation.ds[, c(current.var.names, label)]
      if (importance.options$name == "relieff" |
          importance.options$name == "epistasisrank") {
        importance.options$feature.names <- current.var.names
        new.scores <- getImportanceScores(train.set = new.train.ds,
                                          holdout.set = new.holdout.ds,
                                          label = label,
                                          importance.options = importance.options,
                                          verbose = verbose)
      } else {
        stop("Importance method [", importance.options$name, "] is not recognized")
      }
      q1.scores <- new.scores[[1]]
      q2.scores <- new.scores[[2]]
      if (length(q1.scores) != length(q2.scores)) {
        cat("WARNING: q1 and q2 scores of different lengths [", length(q1.scores), "] vs. [",
            length(q2.scores), "], exiting EC loop\n")
        keep.looping <- FALSE
        next
      }
      diff.scores <- abs(q1.scores - q2.scores)
      delta.q <- max(diff.scores)
      # run the learner
      method.valid <- FALSE
      if (learner.options$name == "randomforest") {
        if (verbose) cat("\trunning randomForest\n")
        rf.train.ds <- train.ds[, c(current.var.names, label)]
        rf.train.pheno <- factor(rf.train.ds[, label])
        if (any(is.na(rf.train.pheno))) {
          cat("Phenos:\n")
          print(rf.train.pheno)
          stop("NAs detected in training phenotype")
        }
        result.rf <- randomForest::randomForest(as.formula(paste(label, "~ .")),
                                                data = rf.train.ds,
                                                ntree = learner.options$rf.ntree,
                                                mtry = learner.options$rf.mtry)
        if (verbose) print(result.rf)
        current.train.acc <- 1 - mean(result.rf$confusion[, "class.error"])

        if (verbose) cat("\tpredict holdout\n")
        #new.holdout.ds[, label] <- factor(new.holdout.ds[, label])
        holdout.pred <- predict(result.rf, newdata = new.holdout.ds)
        if (verbose) print(holdout.pred)
        if (verbose) print(new.holdout.ds[, label])
        current.holdout.acc <- mean(holdout.pred == new.holdout.ds[, label])

        if (!is.null(new.validation.ds)) {
          if (verbose) cat("\tpredict validation\n")
          if (verbose) print(new.validation.ds[, label])
          validation.pred <- predict(result.rf, newdata = new.validation.ds)
          current.validation.acc <- mean(validation.pred == new.validation.ds[, label])
        } else {
          current.validation.acc <- 0
        }
        method.valid <- TRUE
      }
      if (learner.options$name == "xgboost") {
        if (verbose) cat("\trunning xgboost\n")
        xgb.results <- xgboostRF(train.ds = new.train.ds,
                                 holdout.ds = new.holdout.ds,
                                 validation.ds = new.validation.ds,
                                 rf.ntree = learner.options$rf.ntree,
                                 num.rounds = learner.options$num.rounds,
                                 num.threads = learner.options$num.threads,
                                 max.depth =  learner.options$max.depth,
                                 learn.rate = learner.options$learn.rate,
                                 label = label,
                                 verbose = verbose)
        current.train.acc <- xgb.results$algo.acc$train.acc
        current.holdout.acc <- xgb.results$algo.acc$holdout.acc
        current.validation.acc <- xgb.results$algo.acc$validation.acc
        method.valid <- TRUE
      }
      if (!method.valid) {
        stop("Learner method [ ", learner.options$name, " ] is not recognized")
      }
      if (verbose) {
        cat("\tthreshold:          ", threshold, "\n")
        cat("\ttolerance:          ", tolerance, "\n")
        cat("\tlearner function:   ", learner.options$name, "\n")
        cat("\timportance function:", learner.options$name, "\n")
        cat("\tftrain:             ", current.train.acc, "\n")
        cat("\tfholdout:           ", current.holdout.acc, "\n")
      }
      # make adjustments
      lil.bit <- stats::rnorm(1, 0, tolerance)
      if (abs(current.train.acc - current.holdout.acc) < (threshold + lil.bit)) {
        if (verbose) {
          cat("\tClassification error difference is less than threshold + random tolerance\n")
          cat("\tAdjusting holdout error to training error\n")
        }
        current.holdout.acc <- current.train.acc
      } else {
        if (verbose) cat("\tadjust holdout by adding normal random value 0 to tolerance\n")
        current.holdout.acc <- current.holdout.acc + lil.bit
      }
      # save results
      all.train.acc <- c(all.train.acc, current.train.acc)
      all.holdout.acc <- c(all.holdout.acc, current.holdout.acc)
      all.validation.acc <- c(all.validation.acc, current.validation.acc)
      vars.remain.per.update <- c(vars.remain.per.update, length(current.var.names))
      # adjust the temperature
      if (verbose) cat("\tadjusting temperature\n")
      current.temp <- current.temp * exp(-1 / tau) # dropping T
      if (verbose) cat("\t", current.iteration - 1, ': ', current.temp, '\n')
      # check for simulated variables
      if (verbose) cat("\tcollecting results\n")
      if (is.simulated) {
        correct.detect.ec <- c(correct.detect.ec, sum(current.var.names %in% signal.names))
      }
      if (length(current.var.names) < 2) {
        cat("\tvariables list exhausted: no more variables to remove\n")
        keep.looping <- FALSE
      } else {

      }
    }
    current.iteration <- current.iteration + 1
    prev.var.names <- current.var.names
  }
  elapsed.time <- (proc.time() - ptm)[3]
  if (verbose) cat("private EC optimization loop performed", num.updates,
                   "updates in", elapsed.time, " seconds\n")

  # prepare results for returning to caller
  plot.data <- data.frame(vars.remain = vars.remain.per.update,
                          train.acc = all.train.acc,
                          holdout.acc = all.holdout.acc,
                          validation.acc = all.validation.acc,
                          alg = 1)
  if (verbose) print(plot.data)
  plot.data.narrow <- reshape2::melt(plot.data, id = c("vars.remain", "alg"))

  # save the results to an Rdata file if requested
  if (!is.null(save.file)) {
    if (verbose) {
      cat("saving results to [", save.file, "] \n")
    }
    save(plot.data, plot.data.narrow, correct.detect.ec, num.data.cols, num.data.rows,
         signal.names, threshold, tolerance, bias, elapsed.time, file = save.file)
  }

  elapsed.time <- (proc.time() - ptm)[3]
  if (verbose) cat("privateEC elapsed time:", elapsed.time, "\n")

  list(algo.acc = plot.data,
       ggplot.data = plot.data.narrow,
       correct = correct.detect.ec,
       elasped = elapsed.time,
       atts.remain = current.var.names)
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
                      rf.mtry=5,
                      num.rounds=1,
                      num.threads=1,
                      max.depth=6,
                      learn.rate=1,
                      obj.func="binary:logistic",
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
  # ----------------------------------------------------------------
  if (verbose) cat("Running xgboostRF\n")
  ptm <- proc.time()
  var.names <- colnames(train.ds[, -pheno.col])
  train.data <- as.matrix(train.ds[, -pheno.col])
  colnames(train.data) <- var.names
  holdout.data <- as.matrix(holdout.ds[, -pheno.col])
  colnames(holdout.data) <- var.names
  train.pheno <- as.integer(train.ds[, pheno.col]) - 1
  if (verbose) print(table(train.pheno))
  holdout.pheno <- as.integer(holdout.ds[, pheno.col]) - 1
  if (verbose) print(table(holdout.pheno))
  # ----------------------------------------------------------------
  xgb.train.ds <- xgboost::xgb.DMatrix(data = train.data, label = train.pheno)
  xgb.holdo.ds <- xgboost::xgb.DMatrix(data = holdout.data, label = holdout.pheno)
  cat("-----------------------------------\n",
      "eta =",  learn.rate, "\n",
      "nthread =",  num.threads, "\n",
      "nrounds =", num.rounds, "\n",
      "max.depth =", max.depth, "\n",
      "objective =", obj.func, "\n",
      "------------------------------------\n")
  train.model <- xgboost::xgboost(data = xgb.train.ds,
                                  eta = learn.rate,
                                  nthread = num.threads,
                                  nrounds = num.rounds,
                                  max.depth = max.depth,
                                  objective = obj.func)
  if (verbose) cat("Computing variable importance\n")
  importance.xgb <- xgboost::xgb.importance(model = train.model,
                                            feature_names = colnames(xgb.train.ds))
  if (verbose) print(importance.xgb)
  importance.values <- as.numeric(importance.xgb$Gain)
  names(importance.values) <- as.character(importance.xgb$Feature)
  pred.prob <- predict(train.model, xgb.train.ds)
  pred.bin <- as.numeric(pred.prob > 0.5)
  rf.train.accu <- 1 - mean(pred.bin != train.pheno)
  cat("training-accuracy:", rf.train.accu, "\n")
  # ----------------------------------------------------------------
  if (verbose) cat("Predict using the best tree\n")
  pred.prob <- predict(train.model, xgb.holdo.ds)
  pred.bin <- as.numeric(pred.prob > 0.5)
  if (verbose) print(table(pred.bin))
  rf.holdo.accu <- 1 - mean(pred.bin != holdout.pheno)
  cat("holdout-accuracy:", rf.holdo.accu, "\n")
  # ----------------------------------------------------------------
  if (!is.null(validation.ds)) {
    cat("Preparing dating for prediction\n")
    # print(validation.ds[, pheno.col])
    validation.pheno <- as.integer(validation.ds[, pheno.col]) - 1
    # print(validation.pheno)
    validation.data <- as.matrix(validation.ds[, -pheno.col])
    print(dim(validation.data))
    colnames(validation.data) <- var.names
    xgb.validation.ds <- xgboost::xgb.DMatrix(data = validation.data, label = validation.pheno)
    rf.prob <- predict(train.model, xgb.validation.ds)
    rf.bin <- as.numeric(rf.prob > 0.5)
    if (verbose) print(table(rf.bin))
    rf.validation.accu <- 1 - mean(rf.bin != validation.pheno)
    cat("validation-accuracy:", rf.validation.accu, "\n")
  } else {
    rf.validation.accu <- 0
  }
  # ----------------------------------------------------------------
  if (verbose) cat("accuracies", rf.train.accu, rf.holdo.accu, rf.validation.accu, "\n")
  if (!is.null(save.file)) {
    save(rf.validation.accu, rf.holdout.accu, file = save.file)
  }
  xgb.plots <- data.frame(vars.remain = length(var.names),
                          train.acc = rf.train.accu,
                          holdout.acc = rf.holdo.accu,
                          validation.acc = rf.validation.accu,
                          alg = 5)
  melted.xgb <- reshape2::melt(xgb.plots, id = c("vars.remain", "alg"))
  # ----------------------------------------------------------------
  if (file.exists("xgboost.model")) {
    file.remove("xgboost.model")
  }
  elapsed.time <- (proc.time() - ptm)[3]
  if (verbose) cat("xgboostRF elapsed time:", elapsed.time, "\n")
  list(algo.acc = xgb.plots,
       ggplot.data = melted.xgb,
       trn.model = train.model,
       trn.importance = importance.values,
       elasped = elapsed.time)
}
