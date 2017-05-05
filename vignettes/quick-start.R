## ---- echo=FALSE---------------------------------------------------------
library(privateEC)

## ------------------------------------------------------------------------
num.samples <- 100
num.variables <- 100
pb <- 0.1
nbias <- pb * num.variables
signals <- sprintf("gene%04d", 1:nbias)
sim.data <- createSimulation(d=num.variables, n=num.samples,
                             type="sva", verbose=FALSE)
pec.results <- privateEC(data.sets=sim.data, is.simulated=TRUE, n=num.samples,
                         signal.names=signals, verbose=FALSE, update.freq=5)

## ---- echo=FALSE---------------------------------------------------------
knitr::kable(pec.results$plots.data, caption="Algorithm Iterations")

## ---- echo=FALSE, fig.width=6, fig.align='center'------------------------
library(ggplot2)
ggplot(pec.results$melted.data, aes(x=num.atts, y=value, colour=variable)) +
  geom_point(size=1) + geom_line()

