# This is the data from EDASeq vignette

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
load("../../data/bees_data_48.rda")
data = bees_data
treatments = unique(sapply(colnames(data[,-1]), function(x) unlist(strsplit(x,"[.]"))[1]))
colNames = sapply(colnames(data), function(x) unlist(strsplit(x,"[.]"))[1])
baseDir = "48"

for (i in 1:(length(treatments)-1)){
  for (j in (i+1):(length(treatments))){
    outDir = paste0(baseDir, "/", treatments[i],"_",treatments[j])
    ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
    dataSel <- data[,which(colNames %in% c(treatments[i],treatments[j]))]
    # p <- plotMDS(dataSel, col = c("red","red","red","red","red","red","blue","blue","blue","blue","blue","blue"), cex=0.6)
    
    png(filename = paste0(outDir, "/",treatments[i],"_",treatments[j],"_MDS.png"))
    plotMDS(dataSel, col = c("red","red","red","red","red","red","blue","blue","blue","blue","blue","blue"), cex=0.6)
    dev.off()
    
    # dataSel = log(dataSel+1)
    # setDT(dataSel, keep.rownames = TRUE)[]
    # colnames(dataSel)[1] = "ID"
    # dataSel = as.data.frame(dataSel)
    # 
    # plotScatterStatic(dataSel, outDir = outDir)
    # boxSel = dataSel[,-1] %>% gather(Sample,Count)
    # bPlot = ggplot(boxSel, aes(x=Sample, y=Count)) + geom_boxplot()
    # png(filename=paste0(outDir = outDir,"/",treatments[i],"_",treatments[j],"_boxplot.jpg"))
    # bPlot
    # dev.off()
    # plotScatterStatic(dataSel, threshOrth = 0.5, outDir = outDir, option="orthogonal")
    # plotScatterStatic(dataSel, threshOrth = 1, outDir = outDir, option="orthogonal")
    # plotScatterStatic(dataSel, threshOrth = 2, outDir = outDir, option="orthogonal")
    # plotScatterStatic(dataSel, threshOrth = 3, outDir = outDir, option="orthogonal")
    # plotScatterStatic(dataSel, threshOrth = 4, outDir = outDir, option="orthogonal")
    # plotScatterStatic(dataSel, outDir = outDir, option="prediction")
    # plotScatterStatic(dataSel, piLevel=0.99, outDir = outDir, option="prediction")
    # plotScatterStatic(dataSel, piLevel=0.99999, outDir = outDir, option="prediction")
  }
}
