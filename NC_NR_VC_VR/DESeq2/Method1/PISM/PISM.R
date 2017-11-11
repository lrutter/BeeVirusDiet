# Points that are DEG but not outside PI are plotted blue, points that are DEG and outside PI are plotted red, and points that are not DEG but are outside PI are plotted grey.
library(bigPint)

outDir = getwd()
data <- readRDS("../dataRLD.Rds")

setDT(data, keep.rownames = TRUE)[]
colnames(data)[1] = "ID"
data <- as.data.frame(data)

dataMetrics <- readRDS("../dataMetrics.Rds")

plotDEG(data = data, dataMetrics = dataMetrics, outDir = outDir, threshVar = "padj", threshVal = 0.05, option = "scatterPrediction")
