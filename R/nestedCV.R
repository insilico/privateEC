# 
# nested cross-validation for feature-selection/parameter-tuning
# Saeid Parvandeh, April 2018
#
#############################################
##### Consensus nested cross validation #####
#--------------------------------------------

#' Consensus nested cross validation for feature selection and parameter tuning
#' @param train.ds A training data frame with last column as outcome
#' @param validation.ds A validation data frame with last column as outcome
#' @param label A character vector of the outcome variable column name. class/qtrait for classification/regression
#' @param is.simulated A TRUE or FALSE character for data type
#' @param ncv_folds A numeric vector to indicate nested cv folds: c(k_outer, k_inner)
#' @param param.tune A TRUE or FALSE character for tuning parameters
#' @param learning_method Name of the method: glmnet/xgbTree/rf
#' @param importance.algorithm A character vestor containing a specific importance algorithm subtype
#' @param relief.k.method A character of numeric to indicate number of nearest neighbors for relief algorithm.
#' Possible characters are: k_half_sigma (floor((num.samp-1)*0.154)), m6 (floor(num.samp/6)), 
#' myopic (floor((num.samp-1)/2)), and m4 (floor(num.samp/4))
#' @param num_tree Number of trees in random forest and xgboost methods
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list with:
#' \describe{
#'   \item{Train}{Training data accuracy}
#'   \item{Validation}{Validation data accuracy}
#'   \item{Features}{number of variables detected correctly in nested cross validation}
#'   \item{Elapsed}{total elapsed time}
#' } 
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' label <- "class"
#' sim.data <- createSimulation(num.samples = num.samples,
#'                              num.variables = num.variables,
#'                              pct.signals = pct.signals,
#'                              sim.type = "mainEffect",
#'                              label = label,
#'                              verbose = FALSE)                           
#' cnCV.results <- consensus_nestedCV(train.ds = sim.data$train, 
#'                                    validation.ds = sim.data$holdout, 
#'                                    label = label,
#'                                    is.simulated = TRUE,
#'                                    ncv_folds = c(10, 10),
#'                                    param.tune = FALSE,
#'                                    learning_method = "rf",
#'                                    importance.algorithm = "ReliefFbestK",
#'                                    num_tree = 500, 
#'                                    verbose = FALSE)
#' @family nestedCV
#' @export
consensus_nestedCV <- function(train.ds = NULL, 
                               validation.ds = NULL, 
                               label = "class",
                               is.simulated = TRUE,
                               ncv_folds = c(10, 10),
                               param.tune = FALSE,
                               learning_method = "rf",
                               importance.algorithm = "ReliefFequalK",
                               relief.k.method = "k_half_sigma",
                               num_tree = 500, 
                               verbose = FALSE){
  if (is.numeric(relief.k.method)) {
    if (relief.k.method > floor((dim(train.ds)[1]-1)/2)){
      warning("ReliefF k too large. Using maximum.")
      k <- floor((dim(train.ds)[1]-1)/2) 
    } else {
      k <- relief.k.method
    }
    # if someone specifies a numeric value (integer hopefully), use this value for k.
    # However, make sure it is not larger than floor((num.samp.min-1)/2), where
    # num.samp.min  is the min of the train, holdout.. sample sizes.
    # Or you could test the inequality when you encounter each data split.
    # If someone does exceed the threshold, set k to floor((num.samp.min-1)/2) 
    # and writing warning that says 
    # "ReliefF k too large. Using maximum." 
  } else if (relief.k.method ==  "myopic"){
    k <- floor((dim(train.ds)[1]-1)/2)
    # where k may change based on num.samp for train, holdout...
  } else if (relief.k.method ==  "m6") { # default "m6" method
    k <- floor(dim(train.ds)[1]/6)
    # where k may change based on num.samp for train, holdout...
  } else if (relief.k.method == "m4") {
    k <- floor(dim(train.ds)[1]/4)
  } else {
    k <-  floor((dim(train.ds)[1]-1)*0.154)
  }
  if (sum(ncv_folds)>(dim(train.ds)[1])/3){
    stop("There are less than three observations in each fold")
  }
  tune_params <- NULL; accu_vec <- NULL
  Train_accu <- NULL; Test_accu <- NULL
  relief_atts <- list()
  ptm <- proc.time()
  if(verbose){cat("\nRunning consensus nested cross-validation...\n")}
  if(verbose){cat("\n Create [",ncv_folds[1],"] outer folds\n")}
  outer_folds <- caret::createFolds(train.ds[, label], ncv_folds[1], list = FALSE)
  for (i in 1:ncv_folds[1]){
    atts <- list()
    if(verbose){cat("\n Create [",ncv_folds[2],"] inner folds of outer fold[",i,"]\n")}
    inner_folds <- caret::createFolds(train.ds[, label][outer_folds!=i], ncv_folds[2], list = TRUE)
    if(verbose){cat("\n Feature Selection...\n")} 
    for (j in 1:length(inner_folds)){
      inner_idx <- which(outer_folds!=i)[-inner_folds[[j]]]
      rf_scores <- CORElearn::attrEval(label, train.ds[inner_idx, ], 
                                       estimator = importance.algorithm,
                                       costMatrix = NULL, 
                                       outputNumericSplits=FALSE,
                                       kNearestEqual = k)
      atts[[j]] <- names(which(rf_scores>0))
    }
    relief_atts[[i]] <- Reduce(intersect, atts)
    if(param.tune){
      outer_idx <- which(outer_folds!=i)
      trn.data <- as.matrix(train.ds[outer_idx, relief_atts[[i]]])
      tst.data <- as.matrix(train.ds[-outer_idx, relief_atts[[i]]])
      if(label == "class"){
        trn.pheno <- as.factor(train.ds[, label][outer_idx])
      } else {
        trn.pheno <- train.ds[, label][outer_idx]}
      if(label == "class"){
        tst.pheno <- as.factor(train.ds[, label][-outer_idx])
      } else {
        tst.pheno <- train.ds[, label][-outer_idx]}
      if(verbose){cat("\n Parameter tuning...\n")}
      train_model <- caret::train(x = trn.data,
                                  y = trn.pheno,
                                  method = learning_method,
                                  metric = ifelse(is.factor(trn.pheno), "Accuracy", "RMSE"),
                                  trControl = caret::trainControl(method = ifelse(learning_method == "glmnet", "cv", "adaptive_cv"),
                                                           number = 10))
      pred.class <- predict(train_model, tst.data)
      accu <- ifelse(label == "class",1 - mean(pred.class != tst.pheno), cor(tst.pheno, pred.class)^2)
      tune_params <- rbind(tune_params, data.frame(train_model$bestTune, accu))
    }
  }
  nCV_atts <- Reduce(intersect, relief_atts)
  tuneParam <- tune_params[which.max(tune_params$accu), ]
  if(verbose){cat("\n Validating...\n")}
  train.data <- as.matrix(train.ds[, nCV_atts])
  test.data  <- as.matrix(validation.ds[, nCV_atts])
  if(is.simulated && label == "class"){
    train.pheno <- as.integer(train.ds[, label])-1
    test.pheno <- as.integer(validation.ds[, label])-1
  } else if (!is.simulated && label == "class"){
    train.pheno <- as.integer(train.ds[, label])
    test.pheno <- as.integer(validation.ds[, label])
  } else if (label == "qtrait"){
    train.pheno <- train.ds[, label]
    test.pheno <- validation.ds[, label]
  }
  if(verbose){cat("\n Perform [",learning_method,"]\n")}
  if(learning_method == "glmnet"){
    # glmnet - train
    if(param.tune){alpha = tuneParam$alpha; lambda = tuneParam$lambda}else{alpha = 0.5; lambda = min(glmnet.model$lambda)}
    glmnet.model <- glmnet::glmnet(train.data, train.pheno, family = ifelse(label == "class", "binomial", "gaussian"), 
                                   alpha = alpha, lambda = lambda)
    train.pred <- predict(glmnet.model, train.data, s = lambda, type = ifelse(label == "class","class", "response"))
    Train_accu <- ifelse(label == "class" , 1 - mean(as.integer(train.pred) != train.pheno), cor(as.numeric(train.pred), train.pheno)^2)
    # test
    test.pred <- predict(glmnet.model, test.data, s = lambda, type = "class")
    Test_accu <- ifelse(label == "class", 1 - mean(as.integer(test.pred) != test.pheno), cor(test.pred, test.pheno)^2)
  } else if(learning_method == "xgbTree"){
    # xgboost - train
    dtrain <- xgboost::xgb.DMatrix(data = train.data, label = train.pheno)
    dtest <- xgboost::xgb.DMatrix(data = test.data, label = test.pheno)
    if(param.tune){
      shrinkage = tuneParam$eta; max_depth = tuneParam$max_depth; gamma = tuneParam$gamma 
      subsample = tuneParam$subsample; colsample_bytree = tuneParam$colsample_bytree
      min_child_weight = tuneParam$min_child_weight
    } else {
      shrinkage = 0.1; max_depth = 4; gamma = 0; subsample = 1; colsample_bytree = 1;  min_child_weight = 1
    }
    xgb.model <- xgboost::xgboost(data = dtrain,
                                  eta = shrinkage,
                                  nrounds = 4,
                                  max_depth = max_depth,
                                  gamma = gamma,
                                  subsample = subsample,
                                  colsample_bytree = colsample_bytree,
                                  min_child_weight = min_child_weight,
                                  objective = ifelse(label == "class", "binary:logistic", "reg:linear"))
    train.pred.prob <- predict(xgb.model, dtrain)
    train.pred.bin <- as.numeric(train.pred.prob > 0.5)
    Train_accu <- ifelse(label == "class", 1 - mean(train.pred.bin != train.pheno), cor(train.pred.prob, train.pheno)^2)
    # test
    test.pred.prob <- predict(xgb.model, dtest)
    test.pred.bin <- as.numeric(test.pred.prob > 0.5)
    Test_accu <- ifelse(label == "class", 1 - mean(test.pred.bin != test.pheno), cor(test.pred.prob, test.pheno)^2)
  } else if(learning_method == "rf") {
    # random forest - train
    rf.model <- randomForest::randomForest(train.data, 
                                           y = if(label == "class"){as.factor(train.pheno)}else{train.pheno},
                                           mtry = if(param.tune && tuneParam$mtry > 1 && tuneParam$mtry < ncol(train.data))
                                           {tuneParam$mtry} else if (label == "class"){
                                             max(floor(ncol(train.data)/3), 1)
                                           } else {
                                             floor(sqrt(ncol(train.data)))
                                           },
                                           ntree = num_tree)
    Train_accu <- ifelse(label == "class", 1 - mean(rf.model$confusion[, "class.error"]), 
                         cor(as.numeric(as.vector(rf.model$predicted)), train.pheno)^2)
    
    # test
    test.pred <- predict(rf.model, newdata = test.data)
    Test_accu <- ifelse(label == "class", 1 - mean(test.pred != test.pheno),
                        cor(as.numeric(as.vector(test.pred)), test.pheno)^2)
  }
  elapsed.time <- (proc.time() - ptm)[3]
  if(verbose){cat("nestedCV elapsed time", elapsed.time, "\n")}
  list(cv.acc = Train_accu, Validation = Test_accu, Features = nCV_atts, Elapsed = elapsed.time)
}

###########################################
##### Regular nested cross validation #####
#------------------------------------------

#' Regular nested cross validation for feature selection and parameter tuning
#' @param train.ds A training data frame with last column as outcome
#' @param validation.ds A validation data frame with last column as outcome
#' @param label A character vector of the outcome variable column name. class/qtrait for classification/regression
#' @param is.simulated A TRUE or FALSE character for data type
#' @param ncv_folds A numeric vector to indicate nested cv folds: c(k_outer, k_inner)
#' @param param.tune A TRUE or FALSE character for tuning parameters
#' @param learning_method Name of the method: glmnet/xgbTree/rf
#' @param xgb.obj Name of xgboost algorithm
#' @param importance.algorithm A character vestor containing a specific importance algorithm subtype
#' @param relief.k.method A character of numeric to indicate number of nearest neighbors for relief algorithm.
#' Possible characters are: k_half_sigma (floor((num.samp-1)*0.154)), m6 (floor(num.samp/6)), 
#' myopic (floor((num.samp-1)/2)), and m4 (floor(num.samp/4))
#' @param num_tree Number of trees in random forest and xgboost methods
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list with:
#' \describe{
#'   \item{Train}{Training data accuracy}
#'   \item{Validation}{Validation data accuracy}
#'   \item{Features}{number of variables detected correctly in nested cross validation}
#'   \item{Elapsed}{total elapsed time}
#' } 
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' label <- "class"
#' sim.data <- createSimulation(num.samples = num.samples,
#'                              num.variables = num.variables,
#'                              pct.signals = pct.signals,
#'                              sim.type = "mainEffect",
#'                              label = label,
#'                              verbose = FALSE)                         
#' rnCV.results <- regular_nestedCV(train.ds = sim.data$train, 
#'                                  validation.ds = sim.data$holdout, 
#'                                  label = label,
#'                                  is.simulated = TRUE,
#'                                  ncv_folds = c(10, 10),
#'                                  param.tune = FALSE,
#'                                  learning_method = "rf",
#'                                  importance.algorithm = "ReliefFbestK",
#'                                  num_tree = 500, 
#'                                  verbose = FALSE)
#' @family nestedCV
#' @export
regular_nestedCV <- function(train.ds = NULL, 
                             validation.ds = NULL, 
                             label = "class",
                             is.simulated = TRUE,
                             ncv_folds = c(10, 10),
                             param.tune = FALSE,
                             learning_method = "rf",
                             xgb.obj = "binary:logistic",
                             importance.algorithm = "ReliefFequalK",
                             relief.k.method = "k_half_sigma",
                             num_tree = 500, 
                             verbose = FALSE){
  if (is.numeric(relief.k.method)) {
    if (relief.k.method > floor((dim(train.ds)[1]-1)/2)){
      warning("ReliefF k too large. Using maximum.")
      k <- floor((dim(train.ds)[1]-1)/2) 
    } else {
      k <- relief.k.method
    }
    # if someone specifies a numeric value (integer hopefully), use this value for k.
    # However, make sure it is not larger than floor((num.samp.min-1)/2), where
    # num.samp.min  is the min of the train, holdout.. sample sizes.
    # Or you could test the inequality when you encounter each data split.
    # If someone does exceed the threshold, set k to floor((num.samp.min-1)/2) 
    # and writing warning that says 
    # "ReliefF k too large. Using maximum." 
  } else if (relief.k.method ==  "myopic"){
    k <- floor((dim(train.ds)[1]-1)/2)
    # where k may change based on num.samp for train, holdout...
  } else if (relief.k.method ==  "m6") { # default "m6" method
    k <- floor(dim(train.ds)[1]/6)
    # where k may change based on num.samp for train, holdout...
  } else if (relief.k.method == "m4") {
    k <- floor(dim(train.ds)[1]/4)
  } else {
    k <- floor((dim(train.ds)[1]-1)*0.154)
  }
  if (sum(ncv_folds)>(dim(train.ds)[1])/3){
    stop("There are less than three observations in each fold")
  }
  relief_atts <- list(); tune_params <- NULL; 
  Train_accu <- NULL; Test_accu <- NULL
  ptm <- proc.time()
  if(verbose){cat("\nRunning regular nested cross-validation...\n")}
  outer_folds <- caret::createFolds(train.ds[, label], ncv_folds[1], list = FALSE)
  for (i in 1:ncv_folds[1]){
    outer_idx <- which(outer_folds!=i)
    rf_scores <- CORElearn::attrEval(label, train.ds[outer_idx, ], 
                                     estimator = importance.algorithm,
                                     costMatrix = NULL, 
                                     outputNumericSplits=FALSE,
                                     kNearestEqual = k)
    relief_atts[[i]] <- names(which(rf_scores>0))
    trn.data <- as.matrix(train.ds[outer_idx, relief_atts[[i]]])
    tst.data <- as.matrix(train.ds[-outer_idx, relief_atts[[i]]])
    if(label == "class"){
      trn.pheno <- as.factor(train.ds[, label][outer_idx])
    } else {
      trn.pheno <- train.ds[, label][outer_idx]}
    if(label == "class"){
      tst.pheno <- as.factor(train.ds[, label][-outer_idx])
    } else {
      tst.pheno <- train.ds[, label][-outer_idx]}
    train_model <- caret::train(x = trn.data,
                                y = trn.pheno,
                                method = learning_method,
                                metric = ifelse(is.factor(trn.pheno), "Accuracy", "RMSE"),
                                trControl = caret::trainControl(method = ifelse(learning_method == "glmnet", "cv", "adaptive_cv"),
                                                         number = ncv_folds[2]))
    pred.class <- predict(train_model, tst.data)
    accu <- 1 - mean(pred.class != tst.pheno)
    tune_params <- rbind(tune_params, data.frame(train_model$bestTune, accu))
  }
  nCV_atts <- relief_atts[[which.max(tune_params$accu)]]
  tuneParam <- tune_params[which.max(tune_params$accu), ]
  if(verbose){cat("\n Validating...\n")}
  train.data <- as.matrix(train.ds[, nCV_atts])
  test.data  <- as.matrix(validation.ds[, nCV_atts])
  if(is.simulated && label == "class"){
    train.pheno <- as.integer(train.ds[, label])-1
    test.pheno <- as.integer(validation.ds[, label])-1
  } else if (!is.simulated && label == "class"){
    train.pheno <- as.integer(train.ds[, label])
    test.pheno <- as.integer(validation.ds[, label])
  } else if (label == "qtrait"){
    train.pheno <- train.ds[, label]
    test.pheno <- validation.ds[, label]
  }
  if(verbose){cat("\n Perform [",learning_method,"]\n")}
  if(learning_method == "glmnet"){
    # glmnet - train
    if(param.tune){alpha = tuneParam$alpha; lambda = tuneParam$lambda}else{alpha = 0.5; lambda = min(glmnet.model$lambda)}
    glmnet.model <- glmnet::glmnet(train.data, train.pheno, family = ifelse(label == "class", "binomial", "gaussian"), 
                                   alpha = alpha, lambda = lambda)
    train.pred <- predict(glmnet.model, train.data, s = lambda, type = ifelse(label == "class","class", "response"))
    Train_accu <- ifelse(label == "class" , 1 - mean(as.integer(train.pred) != train.pheno), cor(train.pred, train.pheno)^2)
    # test
    test.pred <- predict(glmnet.model, test.data, s = lambda, type = "class")
    Test_accu <- ifelse(label == "class", 1 - mean(as.integer(test.pred) != test.pheno), cor(test.pred, test.pheno)^2)
  } else if(learning_method == "xgbTree"){
    # xgboost - train
    dtrain <- xgboost::xgb.DMatrix(data = train.data, label = train.pheno)
    dtest <- xgboost::xgb.DMatrix(data = test.data, label = test.pheno)
    if(param.tune){
      shrinkage = tuneParam$eta; max_depth = tuneParam$max_depth; gamma = tuneParam$gamma 
      subsample = tuneParam$subsample; colsample_bytree = tuneParam$colsample_bytree
      min_child_weight = tuneParam$min_child_weight
    } else {
      shrinkage = 0.1; max_depth = 4; gamma = 0; subsample = 1; colsample_bytree = 1;  min_child_weight = 1
    }
    xgb.model <- xgboost::xgboost(data = dtrain,
                                  eta = shrinkage,
                                  nrounds = 4,
                                  max_depth = max_depth,
                                  gamma = gamma,
                                  subsample = subsample,
                                  colsample_bytree = colsample_bytree,
                                  min_child_weight = min_child_weight,
                                  objective = xgb.obj)
    train.pred.prob <- predict(xgb.model, dtrain)
    train.pred.bin <- as.numeric(train.pred.prob > 0.5)
    Train_accu <- ifelse(label == "class", 1 - mean(train.pred.bin != train.pheno), cor(train.pred.prob, train.pheno)^2)
    # test
    test.pred.prob <- predict(xgb.model, dtest)
    test.pred.bin <- as.numeric(test.pred.prob > 0.5)
    Test_accu <- ifelse(label == "class", 1 - mean(test.pred.bin != test.pheno), cor(test.pred.prob, test.pheno)^2)
  } else if(learning_method == "rf") {
    # random forest - train
    rf.model <- randomForest::randomForest(train.data, 
                                           y = if(label == "class"){as.factor(train.pheno)}else{train.pheno}, 
                                           mtry = if(param.tune && tuneParam$mtry > 1 && tuneParam$mtry < ncol(train.data))
                                           {tuneParam$mtry} else if (label == "class"){
                                             max(floor(ncol(train.data)/3), 1)
                                           } else {
                                             floor(sqrt(ncol(train.data)))
                                           },
                                           ntree = num_tree)
    Train_accu <- ifelse(label == "class", 1 - mean(rf.model$confusion[, "class.error"]), 
                         cor(as.numeric(as.vector(rf.model$predicted)), train.pheno)^2)
    # test
    test.pred <- predict(rf.model, newdata = test.data)
    Test_accu <- ifelse(label == "class", 1 - mean(test.pred != test.pheno),
                        cor(as.numeric(as.vector(test.pred)), test.pheno)^2)
  }
  elapsed.time <- (proc.time() - ptm)[3]
  if(verbose){cat("nestedCV elapsed time", elapsed.time, "\n")}
  list(cv.acc = Train_accu, Validation = Test_accu, Features = nCV_atts, Elapsed = elapsed.time)
}


