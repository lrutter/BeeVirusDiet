library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(DESeq2)

outDir = "/Users/lindz/BeeVirusDiet/tblshoot-VSNP/Method3/PCP"

# Change group1 and group2 as needed
group1 ="NC"
group2 ="NS"

dds <- readRDS(paste0("/Users/lindz/BeeVirusDiet/tblshoot-VSNP/Method3/comparisonData/", group1, "_", group2, ".rds"))[[1]]
rld <- readRDS(paste0("/Users/lindz/BeeVirusDiet/tblshoot-VSNP/Method3/comparisonData/", group1, "_", group2, ".rds"))[[3]]

sampleIndex <- which(sapply(colnames(assay(rld)), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(group1, group2))

bindataSel <- as.data.frame(assay(rld))[, sampleIndex]
setDT(bindataSel, keep.rownames = TRUE)[]
colnames(bindataSel)[1] <- "ID"
bindataSel$ID <- as.factor(bindataSel$ID)
bindataSel <- as.data.frame(bindataSel)

res <- results(dds, contrast=c("treatment",group1,group2))
degIndex <- which(res@listData$padj<0.05) 
pcpDat <- bindataSel[degIndex,]

nVar = ncol(bindataSel)
colNms <- colnames(bindataSel[, c(2:nVar)])
    
boxDat <- bindataSel[, c(1:nVar)] %>% gather(key, val, -c(ID))
colnames(boxDat) <- c("ID", "Sample", "Count")
pcpDat2 <- pcpDat[, c(1:nVar)] %>% gather(key, val, -c(ID))
colnames(pcpDat2) <- c("ID", "Sample", "Count")

p <- ggplot(boxDat, aes(x = Sample, y = Count)) + geom_boxplot() + geom_line(data=pcpDat2, aes(x = Sample, y = Count, group = ID), size = 0.2, color = "red") + theme(axis.text.x=element_text(angle=90, hjust=1, size = 12))

jpeg(filename=paste0(outDir, "/", group1, "_", group2, ".jpg"), height=650, width=1250)
print(p)
dev.off()
