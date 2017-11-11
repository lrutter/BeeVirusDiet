library(bigPint)

outDir = getwd()
data <- readRDS("../dataRLD.Rds")

setDT(data, keep.rownames = TRUE)[]
colnames(data)[1] = "ID"
data <- as.data.frame(data)

dataMetrics <- readRDS("../dataMetrics.Rds")

plotDEG(data = data, dataMetrics = dataMetrics, outDir = outDir, pointSize = 0.5, xbins = 10, threshVar = "padj", threshVal = 0.05, option = "scatterHexagon")
