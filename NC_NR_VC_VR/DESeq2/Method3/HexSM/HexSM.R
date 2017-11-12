library(bigPint)

outDir = getwd()
dataMetrics <- readRDS("../dataFolder/dataMetrics.Rds")

for (i in 1:length(names(dataMetrics))){
  currPair = names(dataMetrics)[i]
  data <- readRDS(paste0("../dataFolder/dataRLD_", currPair,".Rds"))
  setDT(data, keep.rownames = TRUE)[]
  colnames(data)[1] = "ID"
  data <- as.data.frame(data)
  plotDEG(data = data, dataMetrics = dataMetrics, outDir = outDir, pointSize = 0.5, xbins = 10, threshVar = "padj", threshVal = 0.05, option = "scatterHexagon")
}
