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

for (nC in 2:6){
  plotClusters(data=data, dataMetrics=dataMetrics, nC=nC, threshVar="FDR", threshVal=0.05, outDir=outDir)  
}
