# This is the data from EDASeq vignette

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
load("../../data/bees_data_48.rda")
load("../../data/bees_metrics.rda")
data = bees_data
metrics = bees_metrics
rm(bees_data,bees_metrics)
treatments = unique(sapply(colnames(data[,-1]), function(x) unlist(strsplit(x,"[.]"))[1]))
colNames = sapply(colnames(data), function(x) unlist(strsplit(x,"[.]"))[1])
outDir = "deg-48"

data[,-1] = log(data[,-1]+1)

plotDEG(data, metrics, outDir=outDir, threshVar="FDR", threshVal=0.05)
plotDEG(data, metrics, outDir=outDir, pointSize=1, option = "volcano", threshVar="FDR", threshVal=0.05)
plotDEG(data, metrics, outDir=outDir, option="scatterOrthogonal", threshVar="FDR", threshOrth = 1, threshVal=0.05)
plotDEG(data, metrics, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.05, lineSize = 0.3)
plotDEG(data, metrics, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.01, lineSize = 0.3)
plotDEG(data, metrics, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.001, lineSize = 0.3)
plotDEG(data, metrics, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.0001, lineSize = 0.3)
plotDEG(data, metrics, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.00001, lineSize = 0.3)
