library(bigPint)

outDir = getwd()
data <- readRDS("../dataRLD.Rds")

setDT(data, keep.rownames = TRUE)[]
colnames(data)[1] = "ID"
data <- as.data.frame(data)

dataMetrics <- readRDS("../dataMetrics.Rds")

for (nC in 7:15){
  plotClusters(data=data, dataMetrics=dataMetrics, nC, threshVar="padj", threshVal=0.05, outDir=outDir)  
}
