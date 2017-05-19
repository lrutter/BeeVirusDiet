library(plotly)
library(GGally)
library(hexbin)
library(htmlwidgets)
library(tidyr)
library(shiny)
library(edgeR)
library(EDASeq)
library(dplyr)
library(data.table)
library(ggplot2)

outDir = "/Users/lindz/BeeVirusDiet/scatMatStaticRLog"

my_fn <- function(data, mapping, ...){
  x = data[,c(as.character(mapping$x))]
  y = data[,c(as.character(mapping$y))]
  h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
  hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
  attr(hexdf, "cID") <- h@cID
  p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-1*buffer, maxRange[2]+buffer), ylim = c(-1*buffer, maxRange[2]+buffer))
  p
}

# Following DESeq2 vignette (using rlog)
beeCounts <-read.delim(file="../AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
colnames(beeCounts) <- c("NC.1", "NC.2", "NR.1", "VR.1", "NS.1", "VP.1", "NS.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VS.1", "VS.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VS.3", "NP.5", "VC.5", "VS.4", "NS.3", "VS.5", "VP.5", "NR.3", "NR.4", "VC.6", "NS.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NS.5", "VP.6", "NS.6", "VR.5", "VR.6", "VS.6")
countdata <- beeCounts[ , order(names(beeCounts))]
countdata <- as.matrix(countdata)
coldata = data.frame(row.names = colnames(countdata), virus = unlist(lapply(colnames(countdata), function (x) substring(unlist(strsplit(x, "[.]"))[1],1,1))), diet = unlist(lapply(colnames(countdata), function (x) substring(unlist(strsplit(x, "[.]"))[1],2,2))), treatment = unlist(lapply(colnames(countdata), function (x) unlist(strsplit(x, "[.]"))[1])))
dds = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ treatment)
dds <- DESeq(dds)
rld <- rlog(dds)


bindata <- as.data.frame(assay(rld))
setDT(bindata, keep.rownames = TRUE)[]
colnames(bindata)[1] <- "ID"
bindata$ID <- as.character(bindata$ID)
colNames <- colnames(bindata)
myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
myPairs <- myPairs[-which(myPairs=="ID")]
colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
bindata <- as.data.frame(bindata)

for (i in 1:(length(myPairs)-1)){
  for (j in (i+1):length(myPairs)){
    group1 = myPairs[i]
    group2 = myPairs[j]
    bindataSel <- cbind(ID=bindata$ID, bindata[,which(colGroups %in% c(group1, group2))])
    maxVal = max(bindataSel[,-1])
    minVal = min(bindataSel[,-1])
    maxRange = c(minVal, maxVal)
    xbins=10
    buffer = maxRange[2]/xbins
    p <- ggpairs(bindataSel[,-1], lower = list(continuous = my_fn))
    jpeg(filename=paste0(outDir, "/", group1, "_", group2, ".jpg"), height=700, width=700)
    print(p)
    dev.off()
  }
}




