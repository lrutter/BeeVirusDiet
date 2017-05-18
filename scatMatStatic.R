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

bindata <- readRDS("beeDataFN.rds")
bindata[,2:ncol(bindata)] <- log(bindata[,2:ncol(bindata)])
# Log will cause negative values, so set them to zero
bindata[,2:ncol(bindata)] <- apply(bindata[,2:ncol(bindata)], 2, function(x) ifelse(x < 0, 0, x))

bindata$ID <- as.character(bindata$ID)
colNames <- colnames(bindata)
myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
myPairs <- myPairs[-which(myPairs=="ID")]
colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])

bindataSel <- cbind(ID=bindata$ID, bindata[,which(colGroups %in% c("NC","NP"))])

