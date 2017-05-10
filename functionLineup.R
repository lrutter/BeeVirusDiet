# This function outputs several items into a folder called PermLineup:
# 1) Triangle lineups of the top 100 DEGs
# 2) Correct locations in Correct.csv
# 3) TopDEG1.csv (real data), TopDEG2.csv (permuted), TopDEG3.csv (permuted)
# 4) Permutations.csv tells the permutation orders

# It is different from sigDiffPerm.R because it does not include FDR value in plot outputs.

library(ggplot2)
library(edgeR)
library(dplyr)
library(gtools)
library(tibble)
library(readr)
library(EDASeq)
library(data.table)

# countTable should be a dataframe with first column (type "chr") named "ID" and rest of columns (type "num") named sample names. Have no extra columns. The sample names must be in the format groupName.repNumber (ex: DU.1).

getLineups <- function(countTable, group1, group2, nRep, nPerm, outDir, indScale=FALSE){
  groupNames = unlist(lapply(colnames(countTable), function (x) unlist(strsplit(x, "[.]"))[1]))
  listCond = groupNames[c(which(unlist(lapply(colnames(countTable), function (x) unlist(strsplit(x, "[.]"))[1]))==group1), which(unlist(lapply(colnames(countTable), function (x) unlist(strsplit(x, "[.]"))[1]))==group2))]
  permCols = colnames(countTable)[which(groupNames %in% c(group1, group2))]
  
  allComb <- getPerms(length(listCond))
  allCombLab <- allComb
  for (i in 1:(length(listCond))){
    allCombLab[which(allCombLab == i)] = permCols[i]
  }

  permList <- list()
  correctPlace <- list()

  if (!dir.exists(paste0(getwd(), "/", outDir))){
    dir.create(paste0(getwd(), "/", outDir))
  }
  
  countTable2 = countTable[,which(groupNames %in% c(group1, group2))]
  randPerm = c(1, sample(2:nrow(allCombLab), (nPerm-1)))

  colnames(countTable2) <- allCombLab[randPerm[1],]
  groupNames2 <- unlist(lapply(colnames(countTable2), function (x) unlist(strsplit(x, "[.]"))[1]))
  y <- DGEList(counts=countTable2, group=groupNames2)
  keep <- rowSums(cpm(y)>3) >= 6
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y = estimateCommonDisp(y)
  y = estimateTagwiseDisp(y)
  de = exactTest(y, pair=c(group1,group2))
  tt = topTags(de, n=nrow(y))
  tt <- merge( y[[1]], tt[[1]], by="row.names")
  permList[[1]] = arrange(tt, FDR)
  colnames(permList[[1]])[1] <- "ID"
  write.csv(permList[[1]], file= paste(getwd(), "/", outDir, "/TopDEG1.csv", sep=""))
    
  for(i in 2:nPerm){
    colnames(y[[1]]) = allCombLab[randPerm[i],]
    y[[2]]$group <- as.factor(unlist(lapply(colnames(y[[1]]), function (x) unlist(strsplit(x, "[.]"))[1])))
    attributes(y[[4]])$dimnames[[2]] = allCombLab[randPerm[i],]
    de = exactTest(y, pair=c(group1,group2))
    tt = topTags(de, n=nrow(y))
    
    tt <- merge( y[[1]][,which(groupNames2 %in% c(group1, group2))], tt[[1]], by="row.names")
    permList[[i]] = arrange(tt, FDR)
    colnames(permList[[i]])[1] <- "ID"
    write.csv(permList[[i]], file= paste(getwd(), "/", outDir, "/TopDEG", i, ".csv", sep=""))
  }

  for (i in 1:100){
    fullDat <- data.frame()
    lineup <- permute(seq(1:nPerm))
    correctPlace[i] <- which(lineup==1)
    nRep=6
    for (j in 1:nPerm){
      gene = permList[[j]][i,2:(2*nRep+1)]
      x = unlist(lapply(colnames(gene), function (x) unlist(strsplit(x, "[.]"))[1]))
      x[x==group1] <- 1
      x[x==group2] <- 2
      dat = data.frame(x=x,y=t(gene),z=which(lineup==j))
      colnames(dat)=c("x","y","z")
      dat$x=as.factor(dat$x)
      levels(dat$x)=c(group1,group2)
      dat$meanG1 = mean(filter(dat, x==group1)$y)
      dat$meanG2 = mean(filter(dat, x==group2)$y)
      fullDat <- rbind(fullDat, dat)
    }
    if (indScale){
      allPlot = ggplot(fullDat, aes(x, y)) + geom_point(aes(colour = factor(x)), shape = 20, size=5, alpha = 0.5) + scale_shape(solid = FALSE) + ggtitle(paste("Transcript: ", i)) + ylab("Read Count") + theme(axis.title.x = element_blank(), legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12), plot.title=element_text(hjust=0.5)) + labs(colour = "Group", size=12) + geom_segment(aes(x = 1, y = meanG1, xend = 2, yend = meanG2), colour="gray25", size = 0.1) + facet_wrap(~ z, ncol = 5, scales = "free_y")
    }
    else{
      allPlot = ggplot(fullDat, aes(x, y)) + geom_point(aes(colour = factor(x)), shape = 20, size=5, alpha = 0.5) + scale_shape(solid = FALSE) + ggtitle(paste("Transcript: ", i)) + ylab("Read Count") + scale_y_continuous(limits=c(0, max(fullDat$y))) + theme(axis.title.x = element_blank(), legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12), plot.title=element_text(hjust=0.5)) + labs(colour = "Group", size=12) + geom_segment(aes(x = 1, y = meanG1, xend = 2, yend = meanG2), colour="gray25", size = 0.1) + facet_wrap(~ z, ncol = 5)
    }
    jpeg(file = paste0(getwd(), "/", outDir, "/", "Gene", i, ".jpg"), height = ceiling(nPerm/5)*175, width = 700)
    print(allPlot)
    dev.off()
  }

  correctPlace = data.frame(correctPlace)
  colnames(correctPlace) = 1:100
  correctPlace = t(correctPlace)
  colnames(correctPlace) = "DataPlot"
  correctPlace = data.frame(correctPlace)
  correctPlace <- rownames_to_column(correctPlace, "Gene")
  write.csv(correctPlace, row.names = FALSE, file = paste0(getwd(), "/", outDir, "/Correct.csv"))

  # permInfo --> allComb[1,]
  permInfo <- data.frame(matrix(unlist(allCombLab), nrow = nPerm))
  rownames(permInfo) <- c("Data", paste0(rep("Permute",(nPerm-1)),1:(nPerm-1)))
  colnames(permInfo) <- c(1:(2*nRep))
  write.table(permInfo, sep=",",  col.names=FALSE, file = paste0(getwd(), "/", outDir, "/Permutations.csv"))
}

# Use counts straight from lanes (unflitered and unnormalized)
beeCounts <-read.delim(file="AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
colnames(beeCounts) <- c("NC.1", "NC.2", "NR.1", "VR.1", "NS.1", "VP.1", "NS.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VS.1", "VS.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VS.3", "NP.5", "VC.5", "VS.4", "NS.3", "VS.5", "VP.5", "NR.3", "NR.4", "VC.6", "NS.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NS.5", "VP.6", "NS.6", "VR.5", "VR.6", "VS.6")
countTable <- beeCounts[ , order(names(beeCounts))]

groups <- unique(unlist(lapply(colnames(countTable), function (x) unlist(strsplit(x, "[.]"))[1])))

getLineups(countTable = countTable, group1="NS", group2="NP", nRep=6, nPerm=20, outDir=paste0("/HBPerm","/NS-NP"), indScale=FALSE)
getLineups(countTable = countTable, group1="NC", group2="NP", nRep=6, nPerm=20, outDir=paste0("/HBPerm","/NC-NP"), indScale=FALSE)
