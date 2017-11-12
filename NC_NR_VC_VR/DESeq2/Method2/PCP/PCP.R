library(bigPint)

outDir = getwd()
dataMetrics <- readRDS("../dataMetrics.Rds")

data <- readRDS("../dataRLD.Rds")
setDT(data, keep.rownames = TRUE)[]
colnames(data)[1] = "ID"
data <- as.data.frame(data)

colnames(data) <- c("ID","N.1","N.2","N.3","N.4","N.5","N.6","N.7","N.8","N.9","N.10","N.11","N.12","V.1","V.2","V.3","V.4","V.5","V.6","V.7","V.8","V.9","V.10","V.11","V.12")
plotDEG(data = data, dataMetrics = dataMetrics, outDir = outDir, threshVar = "padj", threshVal = 0.05, lineSize=0.5, option = "parallelCoord")

colnames(data) <- c("ID","C.1","C.2","C.3","C.4","C.5","C.6","R.1","R.2","R.3","R.4","R.5","R.6","C.7","C.8","C.9","C.10","C.11","C.12","R.7","R.8","R.9","R.10","R.11","R.12")
data <- data[,c(1:7,14:19,8:13,20:25)]
plotDEG(data = data, dataMetrics = dataMetrics, outDir = outDir, threshVar = "padj", threshVal = 0.05, lineSize=0.1, option = "parallelCoord")
