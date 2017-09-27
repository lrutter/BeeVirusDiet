# This is the data from EDASeq vignette

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
load("../../data/bees_data_48.rda") # load just to read information, then read data_12 later
data = bees_data
treatments = unique(sapply(colnames(data[,-1]), function(x) unlist(strsplit(x,"[.]"))[1]))
colNames = sapply(colnames(data), function(x) unlist(strsplit(x,"[.]"))[1])
baseDir = "12"
rm(data,bees_data)

for (i in 1:(length(treatments)-1)){
  for (j in (i+1):(length(treatments))){
    outDir = paste0(baseDir, "/", treatments[i],"_",treatments[j])
    ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
    load(paste0("../../data/bees_", treatments[i], "_", treatments[j], "_12.rda"))
    dataSel <- data
    
    png(filename = paste0(outDir, "/",treatments[i],"_",treatments[j],"_MDS.png"))
    plotMDS(dataSel[,-1], col = c("red","red","red","red","red","red","blue","blue","blue","blue","blue","blue"), cex=0.6)
    dev.off()
    
    dataSel[,-1] = log(dataSel[,-1]+1)

    boxSel = dataSel[,-1] %>% gather(Sample,Count)
    bPlot = ggplot(boxSel, aes(x=Sample, y=Count)) + geom_boxplot()
    png(filename=paste0(outDir = outDir,"/",treatments[i],"_",treatments[j],"_boxplot.jpg"))
    bPlot
    dev.off()
    plotScatterStatic(dataSel, threshOrth = 0.5, outDir = outDir, option="orthogonal")
    plotScatterStatic(dataSel, threshOrth = 1, outDir = outDir, option="orthogonal")
    plotScatterStatic(dataSel, threshOrth = 2, outDir = outDir, option="orthogonal")
    plotScatterStatic(dataSel, threshOrth = 3, outDir = outDir, option="orthogonal")
    plotScatterStatic(dataSel, threshOrth = 4, outDir = outDir, option="orthogonal")
    plotScatterStatic(dataSel, outDir = outDir, option="prediction")
    plotScatterStatic(dataSel, piLevel=0.99, outDir = outDir, option="prediction")
    plotScatterStatic(dataSel, piLevel=0.99999, outDir = outDir, option="prediction")
  }
}
