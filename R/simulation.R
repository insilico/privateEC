# simulation.R - Trang Le and Bill White - Fall 2016/Spring 2017
#
# Simulated data for comparison of classification algorithms in:
# Classification algorithms used in the Bioinformatics paper:
# Differential privacy-based Evaporative Cooling feature selection and
# classification with Relief-F and Random Forests
# https://doi.org/10.1093/bioinformatics/btx298

#' Split a data set for machine learning classification
#'
#' Return data.sets as a list of training set, holdout set and validation set
#' according to the predefined percentage of each partition
#' default is a 50-50 split into training and holdout, no testing set
#' code class/label/phenotypes as 1 and -1.
#'
#' @param all.data A data frame of n rows by d colums of data plus a label column
#' @param pct.train A numeric percentage of samples to use for traning
#' @param pct.holdout A numeric percentage of samples to use for holdout
#' @param pct.validation A numeric percentage of samples to use for testing
#' @param class.label A character vector of the data column name for the class label
#' @return A list containing:
#' \describe{
#'   \item{train}{traing data set}
#'   \item{holdout}{holdout data set}
#'   \item{validation}{validation data set}
#' }
#' @examples
#' data("rsfMRIcorrMDD")
#' data.sets <- splitDataset(rsfMRIcorrMDD)
#' @family simulation
#' @export
splitDataset <- function(all.data=NULL,
                         pct.train=0.5,
                         pct.holdout=0.5,
                         pct.validation=0,
                         class.label="phenos") {
  if (is.null(all.data)) {
    # stop or warning and return list of length 0?
    stop("No data passed")
  }
  if (1.0 - (pct.train + pct.holdout + pct.validation) > 0.001 ) {
    stop("Proportions of training, holdout and testing must to sum to 1")
  }
  if (!(class.label %in% colnames(all.data))) {
    stop("Class label is not in the column names or used more than once in data set column names")
  }
  if (!is.factor(all.data[, class.label])) {
    all.data[, class.label] <- factor(all.data[, class.label])
  }
  if (nlevels(all.data[, class.label]) != 2) {
    stop("Cannot split data set with more than or less than 2 factor levels in the class label")
  }

  class.levels <- levels(all.data[, class.label])
  ind.case <- rownames(all.data)[all.data[, class.label] == class.levels[1]]
  ind.ctrl <- rownames(all.data)[all.data[, class.label] == class.levels[2]]

  n.case <- length(ind.case)
  n.ctrl <- length(ind.ctrl)

  n.validation.case <- floor(pct.validation * n.case)
  n.holdo.case <- floor(pct.holdout * n.case)
  n.train.case <- n.case - n.validation.case - n.holdo.case
  partition.case <- sample(c(rep(3, n.validation.case), rep(2, n.holdo.case),
                             rep(1, n.train.case)), n.case)

  n.validation.ctrl <- floor(pct.validation * n.ctrl)
  n.holdo.ctrl <- floor(pct.holdout * n.ctrl)
  n.train.ctrl <- n.ctrl - n.validation.ctrl - n.holdo.ctrl
  partition.ctrl <- sample(c(rep(3, n.validation.ctrl),
                             rep(2, n.holdo.ctrl),
                             rep(1, n.train.ctrl)), n.ctrl)

  all.data <- data.frame(all.data)
  all.data[, class.label] <- factor(all.data[, class.label])
  levels(all.data[, class.label]) <- c(-1, 1)
  data.case <- all.data[ind.case, ]
  data.ctrl <- all.data[ind.ctrl, ]
  X_train <- rbind(data.case[partition.case == 1, ], data.ctrl[partition.ctrl == 1, ])
  X_holdo <- rbind(data.case[partition.case == 2, ], data.ctrl[partition.ctrl == 2, ])
  X_validation <- rbind(data.case[partition.case == 3, ], data.ctrl[partition.ctrl == 3, ])

  # if(nrow(X_validation) == 0) {
  #   data.sets <- list(train=X_train, holdout=X_holdo)
  # } else {
  data.sets <- list(train = X_train, holdout = X_holdo, validation = X_validation)
  # }
  #
  data.sets
}

#' Create a differentially coexpressed data set with interactions
#'
#' @param M An integer for the number of samples (rows)
#' @param N An integer for the number of variables (columns)
#' @param class.label A character vector for the name of the class column
#' @param meanExpression A numeric for the mean expression levels
#' @param A A matrix representing a weighted, undirected network (adjacency)
#' @param randSdNoise Random noise for the background expression levels
#' @param sdNoise A numeric for the noise in the differential expression
#' @param sampleIndicesInteraction A vector of integers of significant variables
#' @param verbose A flag for sending verbose output to stdout
#' @family simulation
#' @return A matrix representing the new new data set.
createInteractions <- function(M=100,
                               N=100,
                               class.label="class",
                               meanExpression=7,
                               A=NULL,
                               randSdNoise=1,
                               sdNoise=0.4,
                               sampleIndicesInteraction=NULL,
                               verbose=FALSE) {
  if (is.null(A)) {
    stop("privacyEC: No adjacency matrix provided")
  }
  if (is.null(sampleIndicesInteraction)) {
    stop("privacyEC: No sample signal indices provided")
  }
  # create a random data matrix
  D <- matrix(nrow = M, ncol = N, data = stats::rnorm(M * N,
                                                      mean = meanExpression,
                                                      sd = randSdNoise))

  # add co-expression
  already_modified <- rep(0, M)
  already_modified[1] <- 1
  for (i in 1:(M - 1)) {
    for (j in (i + 1):M) {
      if (verbose) cat("Condidering A: row", i, "column", j, "\n")
      if ((A[i, j] == 1) && (!already_modified[j])) {
        if (verbose) cat("Making row", j, "from row", i, "\n")
        D[j, ] <- D[i, ] + stats::rnorm(N, mean = 0, sd = as.numeric(sdNoise))
        already_modified[j] <- 1
      } else {
        if (already_modified[j] == 1 && !already_modified[i]) {
          # if j is already modified, we want to modify i,
          # unless i is already modified then do nothing
          D[i, ] <- D[j, ] + stats::rnorm(N, mean = 0, sd = sdNoise)
        }
      }
    }
  }

  # perturb to get differential coexpression
  n1 <- N / 2
  mGenesToPerturb <- length(sampleIndicesInteraction)
  for (i in 1:mGenesToPerturb) {
    geneIdxInteraction <- sampleIndicesInteraction[i]

    # NOTE: bcw, these are never used?
    # g0 <- D[sampleIndicesInteraction, (n1 + 1):N]
    # g1 <- D[sampleIndicesInteraction, 1:n1]

    # get the group 2 gene expression and randomly order for differential coexpression
    x <- D[geneIdxInteraction, 1:n1]
    x <- x[order(stats::runif(length(x)))]
    D[geneIdxInteraction, 1:n1] <- x
  }

  # return a regression ready data frame
  dimN <- ncol(D)
  n1 <- dimN / 2
  n2 <- dimN / 2
  subIds <- c(paste("case", 1:n1, sep = ""), paste("ctrl", 1:n2, sep = ""))
  phenos <- c(rep(1, n1), rep(0, n2))
  newD <- cbind(t(D), phenos)
  colnames(newD) <- c(paste("var", sprintf("%04d", 1:M), sep = ""), class.label)
  rownames(newD) <- subIds

  data.frame(newD)
}

#' Create a simulated data set with main effects
#'
#' \eqn{X = BS + \Gamma G + U}
#'
#' S = Biological group                                                                                                                   m
#' G = Batch
#' U = random error
#'
#' NOTE:  If you use conf=TRUE, then you must have exactly two surrogate
#' variables in the database this function only allows for confounding in
#' the database, not confounding in the new samples
#'
#' @param n.e number of variables
#' @param n.db sample size in database
#' @param n.ns sample size in newsample
#' @param sv.db batches
#' @param sv.ns batches
#' @param sd.b a numeric for sd.b?
#' @param sd.gam a numeric for sd.gam?
#' @param sd.u a numeric for sd.u?
#' @param conf a flag for conf?
#' @param distr.db a numeric for distr.db?
#' @param p.b a numeric for p.b?
#' @param p.gam a numeric for p.gam?
#' @param p.ov a numeric for p.ov?
#' @return A list with:
#' \describe{
#'   \item{db}{database creation variables}
#'   \item{vars}{variables used in simulation}
#' }
#' @references
#' Leek,J.T. and Storey,J.D. (2007) Capturing heterogeneity in gene expression
#' studies by surrogate variable analysis. PLoS Genet., 3, 1724–1735
#' @family simulation
#' @export
createMainEffects <- function(n.e=1000,
                              n.db=70,
                              n.ns=30,
                              sv.db=c("A", "B"),
                              sv.ns=c("A", "B"),
                              sd.b=1,
                              sd.gam=1,
                              sd.u=1,
                              conf=FALSE,
                              distr.db=NA,
                              p.b=0.3,
                              p.gam=0.3,
                              p.ov=p.b / 2) {
  n <- n.db + n.ns
  # Create random error
  U <- matrix(nrow = n.e, ncol = n, stats::rnorm(n.e * n, sd = sd.u))

  # Create index for database vs. new sample #
  ind <- as.factor(c(rep("db", n.db), rep("ns", n.ns)))

  # Create outcome, surrogate variables #
  # Use distr option to show % overlap of outcome, surrogate variables.
  # Note that .5 means no confounding between outcome, surrogate variables.

  # biological variable (fixed at 50% for each outcome)
  S.db <- c(rep(0, round(.5 * n.db)), rep(1, n.db - round(.5 * n.db)))
  S.ns <- c(rep(0, round(.5 * n.ns)), rep(1, n.ns - round(.5 * n.ns)))
  S <- c(S.db, S.ns)

  len0 <- sum(S.db == 0)
  len1 <- sum(S.db == 1)

  if (conf == FALSE) {
    # surrogate variable (no confounding in this function)
    n.sv.db <- length(sv.db)
    prop.db <- 1 / n.sv.db
    # create surrogate variables for outcome 0 in database
    x1 <- c()
    for (i in 1:n.sv.db) {
      x1 <- c(x1, rep(sv.db[i], floor(prop.db * len0)))
    }
    # If the rounding has caused a problem, randomly assign to fill out vector
    while (length(x1) != len0) {
      x1 <- c(x1, sample(sv.db, 1))
    }
    # surrogate variables for outcome 1 will be the same
    # this helps control for the randomly assignment - makes sure there is no
    # added confounding #
    x2 <- x1
  }

  if (conf == TRUE) {
    x1 <- c(rep("A", round(distr.db * len0)),
            rep("B", len0 - round(distr.db * len0)))
    x2 <- c(rep("A", round((1 - distr.db) * len1)),
            rep("B", len1 - round((1 - distr.db) * len1)))
  }

  # create surrogate variables for outcome 0 in new samples
  n.sv.ns <- length(sv.ns)
  prop.ns <- 1 / n.sv.ns

  len0 <- sum(S.ns == 0)
  len1 <- sum(S.ns == 1)

  x3 <- c()
  for (i in 1:n.sv.ns) {
    x3 <- c(x3, rep(sv.ns[i], floor(prop.ns * len0)))
  }
  # If the rounding has caused a problem, randomly assign to fill out vector
  while (length(x3) != len0) {
    x3 <- c(x3, sample(sv.ns, 1))
  }

  # surrogate variables for outcome 1 will be the same
  # this helps control for the randomly assignment - makes sure there is no
  # added confounding
  x4 <- x3
  G <- c(x1, x2, x3, x4)
  G <- t(stats::model.matrix(~ as.factor(G)))[-1, ]
  if (is.null(dim(G))) {
    G <- matrix(G, nrow = 1, ncol = n)
  }

  # Determine which probes are affected by what:
  # 30% for biological, 30% for surrogate, 10% overlap
  # First 30% of probes will be affected by biological signal
  ind.B <- rep(0, n.e)
  ind.B[1:round(p.b * n.e)] <- 1
  # Probes 20% thru 50% will be affected by surrogate variable
  ind.Gam <- rep(0, n.e)
  ind.Gam[round((p.b - p.ov) * n.e):round((p.b - p.ov + p.gam) * n.e)] <- 1

  # figure out dimensions for Gamma

  # create parameters for signal, noise
  B <- matrix(nrow = n.e, ncol = 1, stats::rnorm(n.e, mean = 0, sd = sd.b) * ind.B)
  Gam <- matrix(nrow = n.e, ncol = dim(G)[1],
                stats::rnorm(n.e * dim(G)[1], mean = 0, sd = sd.gam) * ind.Gam)

  # simulate the data
  sim.dat <- B %*% S + Gam %*% G + U
  sim.dat <- sim.dat + abs(min(sim.dat)) + 0.0001

  # simulate data without batch effects
  sim.dat.nobatch <- B %*% S + U
  sim.dat.nobatch <- sim.dat.nobatch + abs(min(sim.dat)) + 0.0001

  # divide parts into database, new samples
  db <- list()
  db$dat <- sim.dat[, ind == "db"]
  db$datnobatch <- sim.dat.nobatch[, ind == "db"]
  db$U <- U[, ind == "db"]
  db$B <- B
  db$S <- S[ind == "db"]
  db$Gam <- Gam
  db$G <- G[ind == "db"]

  vars <- list(n.e  =  n.e, n.db = n.db, n.ns = n.ns, sv.db = sv.db,
               sv.ns = sv.ns, sd.b = sd.b, sd.gam = sd.gam, sd.u = sd.u,
               conf = conf, distr.db = distr.db, p.b = p.b, p.gam = p.gam,
               p.ov = p.ov)

  list(db = db, vars = vars)
}

#' Create a data simulation and return train/holdout/validation data sets.
#'
#' @param num.variables An integer for the number of variables
#' @param num.samples An integer for the number of samples
#' @param pct.signals A numeric for proportion of simulated signal variables
#' @param bias A numeric for effect size in simulated signal variables
#' @param class.label A character vector for the name of the class column
#' @param sim.type A character vector of the type of simulation:
#' mainEffect/interactionErdos/interactionScalefree
#' @param pct.train A numeric percentage of samples to use for traning
#' @param pct.holdout A numeric percentage of samples to use for holdout
#' @param pct.validation A numeric percentage of samples to use for testing
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @param save.file A filename or NULL indicating whether to save the simulations to file
#' @return A list with:
#' \describe{
#'   \item{train}{traing data set}
#'   \item{holdout}{holdout data set}
#'   \item{validation}{validation data set}
#'   \item{class.label}{the class label/column name}
#'   \item{signal.names}{the variable names with simulated signals}
#'   \item{elapsed}{total elapsed time}
#' }
#' @examples
#' num.variables <- 100
#' num.samples <- 100
#' pct.signals <- 0.1
#' bias <- 0.4
#' sim.type <- "mainEffect"
#' sim.data <- createSimulation(num.samples=num.samples,
#'                              num.variables=num.variables,
#'                              pct.signals=pct.signals,
#'                              bias=bias,
#'                              sim.type=sim.type,
#'                              verbose=FALSE)
#' @family simulation
#' @export
createSimulation <- function(num.samples=100,
                             num.variables=100,
                             pct.signals=0.1,
                             bias=0.4,
                             class.label="class",
                             sim.type="mainEffect",
                             pct.train=0.5,
                             pct.holdout=0.5,
                             pct.validation=0,
                             save.file=NULL,
                             verbose=FALSE) {
  ptm <- proc.time()
  nbias <- pct.signals * num.variables
  if (sim.type == "mainEffect") {
    # new simulation:
    # sd.b sort of determines how large the signals are
    # p.b=0.1 makes 10% of the variables signal, bias <- 0.5
    my.sim.data <- simulateData(n.e=num.variables,
                                n.db=num.samples,
                                sd.b=bias,
                                p.b=pct.signals)$db
    dataset <- cbind(t(my.sim.data$datnobatch), my.sim.data$S)
  } else if (sim.type == "interactionScalefree") {
    # interaction simulation: scale-free
    g <- igraph::barabasi.game(num.variables, directed = F)
    A <- igraph::get.adjacency(g)
    myA <- as.matrix(A)
    dataset <- createDiffCoexpMatrixNoME(M=num.variables,
                                         N=num.samples,
                                         meanExpression=7,
                                         A=myA,
                                         randSdNoise=1,
                                         sdNoise=bias,
                                         sampleIndicesInteraction=1:nbias)
  } else if(sim.type == "interactionErdos") {
    attach.prob <- 0.1
    g <- igraph::erdos.renyi.game(num.variables, attach.prob)
    #   foo <- printIGraphStats(g)
    A <- igraph::get.adjacency(g)
    # degrees <- rowSums(A)
    myA <- as.matrix(A)
    dataset <- createDiffCoexpMatrixNoME(M=num.variables,
                                         N=num.samples,
                                         meanExpression=7,
                                         A=myA,
                                         randSdNoise=1,
                                         sdNoise=bias,
                                         sampleIndicesInteraction=1:nbias)
  }
  # make numeric matrix into a data frame for splitting and subsequent ML algorithms
  dataset <- as.data.frame(dataset)
  # BEWARE - DO NOT BE TEMPTED TO USE COLNAMES WITH '.' IN THE NAME 'sim.var1',
  # Removed it altogether; put 'sim' together with 'var' for 'simvar1' instead
  # Not technically allowed but tolerated by R, see make.names() - bcw - 3/4/18
  # Exanple: vgboost does not allow for certain functions using these names
  signal.names <- paste("simvar", 1:nbias, sep = "")
  background.names <- paste("var", 1:(num.variables - nbias), sep = "")
  var.names <- c(signal.names, background.names, class.label)
  colnames(dataset) <- var.names
  split.data <- splitDataset(all.data = dataset,
                             pct.train = pct.train,
                             pct.holdout = pct.holdout,
                             pct.validation = pct.validation,
                             class.label = class.label)

  if (!is.null(save.file)) {
    if (verbose) cat("saving to data/", save.file, ".Rdata\n", sep = "")
    save(split.data, pct.signals, bias, sim.type, file = save.file)
  }

  elapsed <- (proc.time() - ptm)[3]
  if (verbose) cat("createSimulation elapsed time:", elapsed, "\n")

  list(train = split.data$train,
       holdout = split.data$holdout,
       validation = split.data$validation,
       class.label = class.label,
       signal.names = signal.names,
       elapsed = elapsed)
}
