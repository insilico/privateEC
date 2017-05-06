# simulation.R - Trang Le and Bill White - Fall 2016/Spring 2017
# Simulated data for comparison of classification algorithms.

#' Create a differentially coexpressed data set without main effects
#'
#' @param M An integer for the number of samples (rows)
#' @param N An integer for the number of variables (columns)
#' @param meanExpression A numeric for the mean expression levels
#' @param A A matrix representing a weighted, undirected network (adjacency)
#' @param randSdNoise Random noise for the background expression levels
#' @param sdNoise A numeric for the noise in the differential expression
#' @param sampleIndicesInteraction A vector of integers of significant variables
#' @param verbose A flag for sending verbose output to stdout
#' @family simulation
#' @return A matrix representing the new new data set.
createDiffCoexpMatrixNoME <- function(M=100,
                                      N=100,
                                      meanExpression=7,
                                      A=NULL,
                                      randSdNoise=1,
                                      sdNoise=0.4,
                                      sampleIndicesInteraction=NULL,
                                      verbose=FALSE) {
  if(is.null(A)) {
    stop("privacyEC: No adjacency matrix provided")
  }
  if(is.null(sampleIndicesInteraction)) {
    stop("privacyEC: No sample signal indices provided")
  }
  # create a random data matrix
  D <- matrix(nrow=M, ncol=N, data=rnorm(M*N,
                                         mean=meanExpression,
                                         sd=randSdNoise))

  # add co-expression
  already_modified <- rep(0, M)
  already_modified[1] <- 1
  for(i in 1:(M-1)) {
    for(j in (i+1):M) {
      if(verbose) cat("Condidering A: row", i, "column", j, "\n")
      if((A[i, j] == 1) && (!already_modified[j])) {
        if(verbose) cat("Making row", j, "from row", i, "\n")
        D[j, ] <- D[i, ] + rnorm(N, mean=0, sd=as.numeric(sdNoise))
        already_modified[j] <- 1
      } else {
        if(already_modified[j]==1 && !already_modified[i]) {
          # if j is already modified, we want to modify i,
          # unless i is already modified then do nothing
          D[i,] <- D[j,] + rnorm(N, mean=0, sd=sdNoise)
        }
      }
    }
  }

  # perturb to get differential coexpression
  n1 <- N / 2;
  mGenesToPerturb <- length(sampleIndicesInteraction)
  for(i in 1:mGenesToPerturb) {
    geneIdxInteraction <- sampleIndicesInteraction[i]

    g0 <- D[sampleIndicesInteraction, (n1 + 1):N]
    g1 <- D[sampleIndicesInteraction, 1:n1]

    # get the group 2 gene expression and randomly order for differential coexpression
    x <- D[geneIdxInteraction, 1:n1]
    x <- x[order(runif(length(x)))]
    D[geneIdxInteraction, 1:n1] <- x
  }

  # return a regression ready data frame
  dimN <- ncol(D)
  n1 <- dimN / 2
  n2 <- dimN / 2
  subIds <- c(paste("case", 1:n1, sep=""), paste("ctrl", 1:n2, sep=""))
  phenos <- c(rep(1, n1), rep(0, n2))
  newD <- cbind(t(D), phenos)
  colnames(newD) <- c(paste("gene", sprintf("%04d", 1:M), sep=""), "Class")
  rownames(newD) <- subIds
  newD
}

#' Create a simulated data set
#'
#' simulating \eqn{X = BS + \Gamma G + U}
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
#'   \item{db}{?}
#'   \item{new}{?}
#'   \item{varst}{?}
#' }
#' @family simulation
simulateData <- function(n.e=1000,
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
  U <- matrix(nrow=n.e, ncol=n, rnorm(n.e * n, sd=sd.u))

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

  if(conf == FALSE){
    # surrogate variable (no confounding in this function)
    n.sv.db <- length(sv.db)
    prop.db <- 1 / n.sv.db
    # create surrogate variables for outcome 0 in database
    x1 <- c()
    for(i in 1:n.sv.db) {
      x1 <- c(x1, rep(sv.db[i], floor(prop.db * len0)))
    }
    # If the rounding has caused a problem, randomly assign to fill out vector
    while(length(x1) != len0) {
      x1 <- c(x1, sample(sv.db, 1))
    }
    # surrogate variables for outcome 1 will be the same
    # this helps control for the randomly assignment - makes sure there is no
    # added confounding #
    x2 <- x1
  }

  if(conf == TRUE) {
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
  for(i in 1:n.sv.ns) {
    x3 <- c(x3, rep(sv.ns[i], floor(prop.ns * len0)))
  }
  # If the rounding has caused a problem, randomly assign to fill out vector
  while(length(x3) != len0) {
    x3 <- c(x3, sample(sv.ns, 1))
  }

  # surrogate variables for outcome 1 will be the same
  # this helps control for the randomly assignment - makes sure there is no
  # added confounding
  x4 <- x3
  G <- c(x1, x2, x3, x4)
  G <- t(model.matrix(~ as.factor(G)))[-1, ]
  if(is.null(dim(G))) {
    G <- matrix(G, nrow=1, ncol=n)
  }

  # Determine which probes are affected by what:
  # 30% for biological, 30% for surrogate, 10% overlap
  # First 30% of probes will be affected by biological signal
  ind.B <- rep(0, n.e)
  ind.B[1:round(p.b * n.e)] <- 1
  # Probes 20% thru 50% will be affected by surrogate variable
  ind.Gam <- rep(0,n.e)
  ind.Gam[round((p.b-p.ov) * n.e):round((p.b - p.ov + p.gam) * n.e)] <- 1

  # figure out dimensions for Gamma

  # create parameters for signal, noise
  B <- matrix(nrow=n.e, ncol=1, rnorm(n.e, mean=0, sd=sd.b) * ind.B)
  Gam <- matrix(nrow=n.e, ncol=dim(G)[1],
                rnorm(n.e * dim(G)[1], mean=0, sd=sd.gam) * ind.Gam)

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

  vars <- list(n.e=n.e, n.db=n.db, n.ns=n.ns, sv.db=sv.db, sv.ns=sv.ns,
               sd.b=sd.b, sd.gam=sd.gam, sd.u=sd.u, conf=conf,
               distr.db=distr.db, p.b=p.b, p.gam=p.gam, p.ov=p.ov)

  return(list(db=db, new=new, vars=vars))
}

#' Create a data simulation and return train/holdout/test data sets.
#'
#' @param n An integer for the number of samples
#' @param d An integer for the number of variables
#' @param pb A numeric for proportion of functional variables
#' @param bias A numeric for bias in data simulation
#' @param shortname A character vector of a parameters separated by '_'
#' @param type A character vector of the type of simulation: sva|er|inte
#' @param myrun A character vector of a unique run identifier
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @param save.file A flag indicating whther to save the results to file
#' @return A list with:
#' \describe{
#'   \item{plots.data}{data frame of results, a row for each update}
#'   \item{melted.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @examples
#' sim.type <- "sva"
#' num.samples <- 100
#' num.variables <- 100
#' pb <- 0.1
#' bias <- 0.4
#' nbias <- pb * num.variables
#' signals <- sprintf("gene%04d", 1:nbias)
#' sim.data <- createSimulation(d=num.variables, n=num.samples, pb=pb,
#'                              bias=bias, type=sim.type, verbose=FALSE)
#' @family simulation
#' @export
createSimulation <- function(n=100,
                             d=100,
                             pb=0.1,
                             bias=0.4,
                             shortname="paramstring",
                             type="sva",
                             myrun="001",
                             verbose=FALSE,
                             save.file=FALSE) {
  ptm <- proc.time()
  nbias <- pb * d
  if(type == "sva"){
    # new simulation:
    # sd.b sort of determines how large the signals are
    # p.b=0.1 makes 10% of the variables signal, bias <- 0.5
    my.sim.data <- simulateData(n.e=d - 1, n.db=3 * n, sd.b=bias, p.b=pb)$db
    data <- cbind(t(my.sim.data$datnobatch), my.sim.data$S)
  } else if(type == "pri"){
    # old simulation:
    data <- rnorm(n * d * 3, 0, 1)
    data <- matrix(data, n * 3, d)
    data[ , d] <- sign(data[ , d])
    data[, 1:nbias] <- data[, 1:nbias] + bias * data[, d]
  } else if(type == "inte"){
    # interaction simulation: scale-free
    g <- igraph::barabasi.game(d - 1, directed=F)
    A <- igraph::get.adjacency(g)
    myA <- as.matrix(A)
    data <- createDiffCoexpMatrixNoME(M=d - 1,  N=n * 3, meanExpression=7,
                                      A=myA, randSdNoise=1, sdNoise=bias,
                                      1:nbias)
  } else if(type == "er"){
    p <- 0.1
    g <- igraph::erdos.renyi.game(d - 1, p)
    #   foo <- printIGraphStats(g)
    A <- igraph::get.adjacency(g)
    # degrees <- rowSums(A)
    myA <- as.matrix(A)
    data <- createDiffCoexpMatrixNoME(M=d - 1,  N=n * 3, meanExpression=7,
                                      A=myA, randSdNoise=1, sdNoise=bias,
                                      1:nbias)
  }

  ind.case <- sample(3, n * 3, replace=T)[1:floor(n * 3 / 2)]
  ind.ctrl <- sample(ind.case, floor(n * 3 / 2))
  data <- data.frame(data)
  data[, d] <- factor(data[, d])
  levels(data[, d]) <- c(-1, 1)
  colnames(data) <- c(paste("gene", sprintf("%04d", 1:(d - 1)), sep=""), "pheno")
  data.case <- data[data[d] == 1, ]
  data.ctrl <- data[data[d] == -1, ]
  X_train <- rbind(data.case[ind.case == 1,], data.ctrl[ind.ctrl == 1,])
  X_holdo <- rbind(data.case[ind.case == 2,], data.ctrl[ind.ctrl == 2,])
  X_test  <- rbind(data.case[ind.case == 3,], data.ctrl[ind.ctrl == 3,])

  if(save.file) {
    myfile <- paste("data/", type, "_", shortname, "_data.Rdata", sep="")
    if(verbose) cat("saving to data/", myfile, ".Rdata\n", sep="")
    save(n, d, pb, X_train, X_holdo, X_test, bias, type,
         shortname, file=myfile)
  }

  elapsed <- (proc.time() - ptm)[3]
  if(verbose) cat("createSimulation elapsed time:", elapsed, "\n")

  list(train=X_train, holdout=X_holdo, test=X_test, elapsed=elapsed)
}

#' Write inbix numeric and phenotype files (PLINK format)
#'
#' @param data.sets A list of train, holdout and test data frames
#' @param base.sim.prefix A character vector for the input and saved file prefixes
#' @param verbose A flag for sending berbose output to stdout
#' @return NULL
saveSimAsInbixNative <- function(data.sets=NULL,
                                 base.sim.prefix,
                                 verbose=FALSE) {
  if(is.null(data.sets)) {
    stop("privateEC: No data sets provided as first argument")
  }
  X_train <- data.sets$train
  X_holdo <- data.sets$holdout
  X_test <- data.sets$test

  # train
  train.expr.matrix <- X_train[, 1:(ncol(X_train)-1)]
  var.names <- colnames(train.expr.matrix)
  train.num.subj <- nrow(train.expr.matrix)
  train.subj.names <- paste("subj", 1:train.num.subj, sep="")
  train.phenotype <- ifelse(X_train[, ncol(X_train)] == -1, 0, 1)
  train.inbix <- cbind(train.subj.names, train.subj.names, train.expr.matrix)
  colnames(train.inbix) <- c("FID", "IID", var.names)
  write.table(train.inbix, file=paste(base.sim.prefix, ".train.sim.num", sep=""),
              quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
  train.inbix.pheno <- cbind(train.subj.names, train.subj.names, train.phenotype)
  write.table(train.inbix.pheno, file=paste(base.sim.prefix, ".train.sim.pheno", sep=""),
              quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
  # holdout
  holdo.expr.matrix <- X_holdo[, 1:(ncol(X_holdo)-1)]
  var.names <- colnames(holdo.expr.matrix)
  holdo.num.subj <- nrow(holdo.expr.matrix)
  holdo.subj.names <- paste("subj", 1:holdo.num.subj, sep="")
  holdo.phenotype <- ifelse(X_holdo[, ncol(X_holdo)] == -1, 0, 1)
  holdo.inbix <- cbind(holdo.subj.names, holdo.subj.names, holdo.expr.matrix)
  colnames(holdo.inbix) <- c("FID", "IID", var.names)
  write.table(holdo.inbix, file=paste(base.sim.prefix, ".holdo.sim.num", sep=""),
              quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
  holdo.inbix.pheno <- cbind(holdo.subj.names, holdo.subj.names, holdo.phenotype)
  write.table(holdo.inbix.pheno, file=paste(base.sim.prefix, ".holdo.sim.pheno", sep=""),
              quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
  # test
  test.expr.matrix <- X_test[, 1:(ncol(X_test)-1)]
  var.names <- colnames(test.expr.matrix)
  test.num.subj <- nrow(test.expr.matrix)
  test.subj.names <- paste("subj", 1:test.num.subj, sep="")
  test.phenotype <- ifelse(X_test[, ncol(X_test)] == -1, 0, 1)
  test.inbix <- cbind(test.subj.names, test.subj.names, test.expr.matrix)
  colnames(test.inbix) <- c("FID", "IID", var.names)
  write.table(test.inbix, file=paste(base.sim.prefix, ".test.sim.num", sep=""),
              quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
  test.inbix.pheno <- cbind(test.subj.names, test.subj.names, test.phenotype)
  write.table(test.inbix.pheno, file=paste(base.sim.prefix, ".test.sim.pheno", sep=""),
              quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
}
