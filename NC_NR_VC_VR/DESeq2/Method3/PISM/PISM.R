# Points that are DEG but not outside PI are plotted blue, points that are DEG and outside PI are plotted red, and points that are not DEG but are outside PI are plotted grey.
library(bigPint)

outDir = getwd()
dataMetrics <- readRDS("../dataFolder/dataMetrics.Rds")

for (i in 1:length(names(dataMetrics))){
  currPair = names(dataMetrics)[i]
  data <- readRDS(paste0("../dataFolder/dataRLD_", currPair,".Rds"))
  setDT(data, keep.rownames = TRUE)[]
  colnames(data)[1] = "ID"
  data <- as.data.frame(data)
  plotDEG(data = data, dataMetrics = dataMetrics, outDir = outDir, threshVar = "padj", threshVal = 0.05, option = "scatterPrediction")
}
