library(plotly)
library(GGally)
library(hexbin)
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(limma)
library(Glimma)
library(edgeR)

rm(list=ls())

outDir = "/Users/lindz/bigPint/tblshoot/CheckReps"

my_fn <- function(data, mapping, ...){
  xChar = as.character(mapping$x)
  yChar = as.character(mapping$y)
  x = data[,c(xChar)]
  y = data[,c(yChar)]
  h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
  hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
  attr(hexdf, "cID") <- h@cID
  p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-1*buffer, maxRange[2]+buffer), ylim = c(-1*buffer, maxRange[2]+buffer))
  p
}

dds <- readRDS("/Users/lindz/bigPint/tblshoot/AllPairs/data_limma.Rds")

# Change group1 and group2 as needed
group1 ="VC"

sampleIndex <- which(sapply(colnames(dds[[1]]), function(x) unlist(strsplit(x,"[.]"))[1]) %in% group1)

bindataSel <- as.data.frame(dds[[1]])[, sampleIndex]
setDT(bindataSel, keep.rownames = TRUE)[]
colnames(bindataSel)[1] <- "ID"
bindataSel$ID <- as.character(bindataSel$ID)
bindataSel <- as.data.frame(bindataSel)

x <- DGEList(counts=bindataSel[,2:ncol(bindataSel)])
bindataSelX <- rlog(x[[1]])

# Take logs
bindataSelX <- as.data.frame(bindataSelX)
bindataSelX$ID <- bindataSel$ID
bindataSelX2 <- bindataSelX[,c(7,1,2,3,4,5,6)]
bindataSel <- bindataSelX2

maxVal = max(bindataSel[,-1])
minVal = min(bindataSel[,-1])
maxRange = c(minVal, maxVal)
xbins=10
buffer = maxRange[2]/xbins
p <- ggpairs(bindataSel[,-1], lower = list(continuous = my_fn)) + theme(strip.text = element_text(size = 20))
jpeg(filename=paste0(outDir, "/", group1, "_rlog_", xbins, ".jpg"), height=1400, width=1400)
print(p)
dev.off()

