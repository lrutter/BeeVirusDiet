library(plotly)
library(GGally)
library(hexbin)
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(DESeq2)

# Look at VS and NP normalized/filtered to entire dataset

outDir = "/Users/lindz/BeeVirusDiet/tblshoot-VSNP/Method3/HexSM"

my_fn <- function(data, mapping, ...){
  xChar = as.character(mapping$x)
  yChar = as.character(mapping$y)
  x = data[,c(xChar)]
  y = data[,c(yChar)]
  h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
  hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
  attr(hexdf, "cID") <- h@cID
  p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-1*buffer, maxRange[2]+buffer), ylim = c(-1*buffer, maxRange[2]+buffer)) + geom_point(data = degData, aes_string(x=xChar, y=yChar), inherit.aes = FALSE, color = "orange", size = 4)
  p
}

# Change group1 and group2 as needed
group1 ="NC"
group2 ="NS"

dds <- readRDS(paste0("/Users/lindz/BeeVirusDiet/tblshoot-VSNP/Method3/comparisonData/", group1, "_", group2, ".rds"))[[1]]
rld <- readRDS(paste0("/Users/lindz/BeeVirusDiet/tblshoot-VSNP/Method3/comparisonData/", group1, "_", group2, ".rds"))[[3]]

sampleIndex <- which(sapply(colnames(assay(rld)), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(group1, group2))

bindataSel <- as.data.frame(assay(rld))[, sampleIndex]
setDT(bindataSel, keep.rownames = TRUE)[]
colnames(bindataSel)[1] <- "ID"
bindataSel$ID <- as.character(bindataSel$ID)
bindataSel <- as.data.frame(bindataSel)

res <- results(dds, contrast=c("treatment",group1,group2))
degIndex <- which(res@listData$padj<0.05) 
degData <- bindataSel[degIndex,]

maxVal = max(bindataSel[,-1])
minVal = min(bindataSel[,-1])
maxRange = c(minVal, maxVal)
xbins=10
buffer = maxRange[2]/xbins
p <- ggpairs(bindataSel[,-1], lower = list(continuous = my_fn))
jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_larger.jpg"), height=1400, width=1400)
print(p)
dev.off()

