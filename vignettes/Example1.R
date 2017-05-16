## ------------------------------------------------------------------------
library(privateEC)
update.freq <- 10
num.samples <- 100
num.variables <- 100
prob.bias <- 0.1
nbias <- prob.bias * num.variables
signals <- sprintf("gene%04d", 1:nbias)
sim.data <- createSimulation(d=num.variables, 
                             n=num.samples,
                             type="inte", 
                             verbose=FALSE)
pec.results <- privateEC(data.sets=sim.data, 
                         is.simulated=TRUE, 
                         n=num.samples,
                         signal.names=signals, 
                         verbose=FALSE, 
                         update.freq=update.freq)

## ---- echo=FALSE---------------------------------------------------------
knitr::kable(pec.results$plots.data, caption="Algorithm Iterations",
             row.names=FALSE, digits=3)

## ---- echo=FALSE, fig.width=14, fig.width=7, fig.align='center'----------
# library(ggplot2)
# ggplot(pec.results$melted.data, aes(x=num.atts, y=value, colour=variable)) +
#   geom_point(size=1) + geom_line()
plot(pec.results$plots.data$num.atts, 
     pec.results$plots.data$fholdo, 
     col="red", pch=16, type='b', cex=0.75, 
     main="One run of privateEC",
     xlab="Number of Attributes in Model",
     ylab="Accuracy")
points(pec.results$plots.data$num.atts, 
       pec.results$plots.data$ftrain, 
       col="green", pch=1, type='b', cex=0.75)
points(pec.results$plots.data$num.atts, 
       pec.results$plots.data$ftest, 
       col="blue", pch=4, type='b', cex=0.75)
legend("topright", c("Train", "Holdout", "Test"), pch=c(16, 1, 4), col=c("red", "green", "blue"), cex=0.75)

