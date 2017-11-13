library(bigPint)
library(data.table)

outDir = getwd()
data <- readRDS("../data.Rds")

data <- as.data.frame(data)
setDT(data, keep.rownames = TRUE)[]
colnames(data)[1] = "ID"
data <- as.data.frame(data)
data[,-1] <- log(data[,-1]+1)

dataMetrics <- readRDS("../dataMetrics.Rds")

plotDEG(data = data, dataMetrics = dataMetrics, outDir = outDir, pointSize = 0.5, xbins = 10, threshVar = "FDR", threshVal = 0.05, option = "scatterHexagon")
