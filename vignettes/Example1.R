## ------------------------------------------------------------------------
library(privateEC)
n.samples <- 100
n.variables <- 100
bias <- 0.4
type <- "mainEffect"
pct.signals <- 0.1
update.freq <- 5
verbose <- FALSE

data.sets <- createSimulation(num.samples = n.samples,
                              num.variables = n.variables,
                              pct.signals = pct.signals,
                              bias = bias,
                              pct.train = 1 / 3,
                              pct.holdout = 1 / 3,
                              pct.validation = 1 / 3,
                              sim.type = type,
                              save.file = NULL,
                              verbose = verbose)

pec.result <- privateEC(train.ds = data.sets$train,
                        holdout.ds = data.sets$holdout,
                        validation.ds = data.sets$validation,
                        label = data.sets$class.label,
                        is.simulated = TRUE,
                        bias = bias,
                        update.freq = update.freq,
                        save.file = NULL,
                        signal.names = data.sets$signal.names,
                        verbose = verbose)

## ---- echo=FALSE---------------------------------------------------------
knitr::kable(pec.result$algo.acc, caption="Algorithm Iterations",
             row.names=FALSE, digits=3)

## ---- echo=FALSE, fig.width=7, fig.width=7, fig.align='center'-----------
# library(ggplot2)
# ggplot(pec.result$melted.data, aes(x=num.atts, y=value, colour=variable)) +
#   geom_point(size=1) + geom_line()
plot(pec.result$algo.acc$vars.remain, 
     pec.result$algo.acc$holdout.acc, 
     col="red", pch=16, type='b', cex=0.75, 
     main="One run of privateEC",
     ylim=c(0.05, 1.0), 
     xlab="Number of Attributes in Model",
     ylab="Accuracy")
points(pec.result$algo.acc$vars.remain, 
       pec.result$algo.acc$train.acc, 
       col="green", pch=1, type='b', cex=0.75)
points(pec.result$algo.acc$vars.remain, 
       pec.result$algo.acc$validation.acc, 
       col="blue", pch=4, type='b', cex=0.75)
legend("topright", c("Train", "Holdout", "Test"), 
       pch=c(16, 1, 4), col=c("red", "green", "blue"), cex=0.75)

