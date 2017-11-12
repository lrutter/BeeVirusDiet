library(bigPint)

outDir = getwd()
dataMetrics <- readRDS("../dataFolder/dataMetrics.Rds")

for (i in 3:length(names(dataMetrics))){
  currPair = names(dataMetrics)[i]
  data <- readRDS(paste0("../dataFolder/dataRLD_", currPair,".Rds"))
  setDT(data, keep.rownames = TRUE)[]
  colnames(data)[1] = "ID"
  data <- as.data.frame(data)
  for (nC in 2:15){  
    plotClusters(data=data, dataMetrics=dataMetrics, nC, threshVar="padj", threshVal=0.05, outDir=outDir)  
  }
}
  