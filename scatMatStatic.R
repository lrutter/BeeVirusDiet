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

outDir = "/Users/lindz/BeeVirusDiet/scatMatStatic"

my_fn <- function(data, mapping, ...){
  x = data[,c(as.character(mapping$x))]
  y = data[,c(as.character(mapping$y))]
  h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
  hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
  attr(hexdf, "cID") <- h@cID
  p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-1*buffer, maxRange[2]+buffer), ylim = c(-1*buffer, maxRange[2]+buffer))
  p
}

bindata <- readRDS("beeDataFN.rds")
bindata[,2:ncol(bindata)] <- log(bindata[,2:ncol(bindata)])
# Log will cause negative values, so set them to zero
bindata[,2:ncol(bindata)] <- apply(bindata[,2:ncol(bindata)], 2, function(x) ifelse(x < 0, 0, x))
bindata$ID <- as.character(bindata$ID)
colNames <- colnames(bindata)
myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
myPairs <- myPairs[-which(myPairs=="ID")]
colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])

for (i in 1:(length(myPairs)-1)){
  for (j in (i+1):length(myPairs)){
    group1 = myPairs[i]
    group2 = myPairs[j]
    bindataSel <- cbind(ID=bindata$ID, bindata[,which(colGroups %in% c(group1, group2))])
    maxVal = max(abs(bindataSel[,-1]))
    maxRange = c(0, maxVal)
    xbins=10
    buffer = maxRange[2]/xbins
    p <- ggpairs(bindataSel[,-1], lower = list(continuous = my_fn))
    jpeg(filename=paste0(outDir, "/", group1, "_", group2, ".jpg"), height=700, width=700)
    print(p)
    dev.off()
  }
}




