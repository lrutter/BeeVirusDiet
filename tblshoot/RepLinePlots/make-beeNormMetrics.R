library(edgeR)
library(ggplot2)
library(GGally)
library(EDASeq)
library(utils)
library(data.table)

beeCounts <-read.delim(file="../../data/AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
colnames(beeCounts) <- c("NC.1", "NC.2", "NR.1", "VR.1", "NS.1", "VP.1", "NS.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VS.1", "VS.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VS.3", "NP.5", "VC.5", "VS.4", "NS.3", "VS.5", "VP.5", "NR.3", "NR.4", "VC.6", "NS.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NS.5", "VP.6", "NS.6", "VR.5", "VR.6", "VS.6")
beeCounts <- beeCounts[ , order(names(beeCounts))]
y <- DGEList(counts=beeCounts)
minLib <- min(y$samples$lib.size)
keep <- rowSums(cpm(y)>3) >= 24
# Number of genes 15,314--> 8,581
y <- y[keep, , keep.lib.sizes=FALSE]
y[[1]] <- betweenLaneNormalization(y[[1]], which="full", round=FALSE)
Group = sapply(colnames(y), function(x) unlist(strsplit(x,"[.]"))[1])
design <- model.matrix(~0+Group, data=y$samples)
colnames(design) <- levels(Group)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)


bees_metrics <- list()

for (i in 1:(ncol(fit)-1)){
  for (j in (i+1):ncol(fit)){
    contrast=rep(0,ncol(fit))
    contrast[i]=1
    contrast[j]=-1
    lrt <- glmLRT(fit, contrast=contrast)
    lrt <- topTags(lrt, n = nrow(y[[1]]))[[1]]
    
    setDT(lrt, keep.rownames = TRUE)[]
    colnames(lrt)[1] = "ID"
    lrt <- as.data.frame(lrt)
    
    bees_metrics[[paste0(unique(Group)[i], "_", unique(Group)[j])]] <- lrt
  }
}

save(bees_metrics, file = "../../data/bees_metrics.rda")
bees_data = data.frame(y[[1]])
setDT(bees_data, keep.rownames = TRUE)[]
colnames(bees_data)[1] = "ID"
bees_data = as.data.frame(bees_data)

save(bees_data, file = "../../data/bees_data.rda")
