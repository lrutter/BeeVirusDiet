library(edgeR)
library(ggplot2)
library(GGally)
library(EDASeq)
library(utils)

# In this script, I calculate the number of DEGs between all pairwise combinations of treatments. However, instead of feeding in the whole data frame into the DEG analysis, I only feed in the subset of data that contains the treatment pair of interest.

findPairs <- function(data){
  y <- DGEList(counts=data)
  minLib <- min(y$samples$lib.size)
  keep <- rowSums(cpm(y)>3) >= 6
  y <- y[keep, , keep.lib.sizes=FALSE]
  y[[1]] <- betweenLaneNormalization(y[[1]], which="full", round=FALSE)
  Group = sapply(colnames(y), function(x) unlist(strsplit(x,"[.]"))[1])
  design <- model.matrix(~0+Group, data=y$samples)
  colnames(design) <- levels(Group)
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  contrast=rep(0,ncol(fit))
  contrast[1]=1
  contrast[2]=-1
  lrt <- glmLRT(fit, contrast=contrast)
  lrt <- topTags(lrt, n = nrow(y[[1]]))[[1]]
  setDT(lrt, keep.rownames = TRUE)[]
  colnames(lrt)[1] = "ID"
  lrt <- as.data.frame(lrt)
  lrt
}

thisPath <- "/Users/lindz/BeeVirusDiet/edgeR-btwnLane"
beeCounts <-read.delim(file="../../data/AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
colnames(beeCounts) <- c("NC.1", "NC.2", "NR.1", "VR.1", "NS.1", "VP.1", "NS.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VS.1", "VS.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VS.3", "NP.5", "VC.5", "VS.4", "NS.3", "VS.5", "VP.5", "NR.3", "NR.4", "VC.6", "NS.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NS.5", "VP.6", "NS.6", "VR.5", "VR.6", "VS.6")
beeCounts <- beeCounts[ , order(names(beeCounts))]

treatments = unique(sapply(colnames(beeCounts), function(x) unlist(strsplit(x,"[.]"))[1]))
colNames = sapply(colnames(beeCounts), function(x) unlist(strsplit(x,"[.]"))[1])

allPairs = data.frame(Treatment1 = factor(), Treatment2 = factor(), NumberDEG = numeric())
metrics=list()
for (i in 1:(length(treatments)-1)){
  for (j in (i+1):length(treatments)){
    keep = which(colNames %in% c(treatments[i], treatments[j]))
    beeCounts2 = beeCounts[, keep]
    lrt <- findPairs(beeCounts2)
    metrics[[paste0(treatments[i], "_", treatments[2])]] <- lrt
    
    lrt <- lrt[which(lrt$FDR<0.05),]
    lrtLength <- nrow(lrt)
    allPairs = rbind(allPairs, data.frame(Treatment1 = factor(treatments[i]), Treatment2 = factor(treatments[j]), NumberDEG = lrtLength))
  }
}







X = NULL
remVals = c(1,2,3,4,5,6,7,8,9,10,11,12)

for (j in 1:12){
  for (i in 1:ncol(combn(remVals,j))){
    beeCounts2 = beeCounts[,-combn(remVals,j)[,i]]
    x = findPairs(beeCounts2)
    X = rbind(X, x)
  }
}

R = NULL
for (j in 1:12){
  for (i in 1:ncol(combn(remVals,j))){
    combn(remVals,j)[,i]
    rowName=colnames(beeCounts)[combn(remVals,j)[,i]]
    rowName=paste(rowName, sep="", collapse="_")
    R = rbind(R, rowName)
  }
}

rownames(X) = R
colnames(X) = paste0(allPairs$Treatment1, "_", allPairs$Treatment2)
remTable = X
remTable2 = remTable[order(rowSds(remTable)),]
saveRDS(remTable2,"remTable.rds")
