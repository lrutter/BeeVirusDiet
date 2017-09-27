# This is the data from EDASeq vignette

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
load("../../data/bees_data_48.rda") # load just to read information, then read data_12 later
data = bees_data
treatments = unique(sapply(colnames(data[,-1]), function(x) unlist(strsplit(x,"[.]"))[1]))
colNames = sapply(colnames(data), function(x) unlist(strsplit(x,"[.]"))[1])
baseDir = "deg-12"
load("../../data/bees_metrics_12.rda")
metrics = metrics12
rm(data,bees_data,metrics12)

for (i in 1:(length(treatments)-1)){
  for (j in (i+1):(length(treatments))){
    outDir = paste0(baseDir, "/", treatments[i],"_",treatments[j])
    ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
    load(paste0("../../data/bees_", treatments[i], "_", treatments[j], "_12.rda"))
    data[,-1] = log(data[,-1]+1)
    
    #plotDEG(data, metrics, outDir=outDir, threshVar="FDR", threshVal=0.05)
    plotDEG(data, metrics, outDir=outDir, threshVar="FDR", threshVal=0.01)
    plotDEG(data, metrics, outDir=outDir, threshVar="FDR", threshVal=0.001)
    # plotDEG(data, metrics, outDir=outDir, pointSize=1, option = "volcano", threshVar="FDR", threshVal=0.05)
    # plotDEG(data, metrics, outDir=outDir, option="scatterOrthogonal", threshVar="FDR", threshOrth = 1, threshVal=0.05)
    # plotDEG(data, metrics, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.05, lineSize = 0.3)
    # plotDEG(data, metrics, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.01, lineSize = 0.3)
    # plotDEG(data, metrics, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.001, lineSize = 0.3)
    # plotDEG(data, metrics, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.0001, lineSize = 0.3)
    # plotDEG(data, metrics, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.00001, lineSize = 0.3)    
  }
}
