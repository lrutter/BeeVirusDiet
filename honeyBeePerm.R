library(edgeR)

# Generate original count table (beeDataFN.rds)
beeCounts <-read.delim(file="AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
colnames(beeCounts) <- c("NC.1", "NC.2", "NR.1", "VR.1", "NS.1", "VP.1", "NS.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VS.1", "VS.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VS.3", "NP.5", "VC.5", "VS.4", "NS.3", "VS.5", "VP.5", "NR.3", "NR.4", "VC.6", "NS.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NS.5", "VP.6", "NS.6", "VR.5", "VR.6", "VS.6")
beeCounts <- beeCounts[ , order(names(beeCounts))]
nRep=6
group = rep(c("NC","NP","NR","NS","VC","VP","VR","VS"), each = nRep)
y <- DGEList(counts=beeCounts, group=group)
keep <- rowSums(cpm(y)>3) >= 24
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y = estimateCommonDisp(y)
y = estimateTagwiseDisp(y)
de = exactTest(y, pair=c("NR","VP"))
tt = topTags(de, n=nrow(y))

# GLM of NR-VP had 113 DEGs (with FDR <0.05). Now, exactTest of NR-VP had 65 DEGs.
# NRVP is the GLM result of comparing NR and VP
NRVP <- readRDS("/Users/lindz/BeeVirusDiet/edgeR/DEGPairs/NRVP.rds")
#Slight differences in FDR values and orders if do exactTest() or GLM for NR and VP.
# 93 of the 113 DEGs from GLM that had FDR<0.05 are in the top 113 rows from exactTest tt
sum(rownames(NRVP) %in% rownames(tt[[1]][1:nrow(NRVP),]))
#55 of the 65 DEGs from exactTest tt are in the top 65 DEGs from GLM
sum(rownames(tt[[1]][1:65,]) %in% rownames(NRVP[1:65,]))

# Now check to see if permutation changes values in each row. If not, then permutation function does not need to rerun calcNormFactors(), estimateCommonDisp(), and estimateTagwiseDisp() for each permutation.
tt <- merge( y[[1]][,which(group %in% c("NR", "VP"))], tt[[1]], by="row.names")
# Save the non-permuted normalized and filtered data
NRVP_NoPerm <- tt


set.seed(10)
beeCounts <-read.delim(file="AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
colNames <- c("NC.1", "NC.2", "NR.1", "VR.1", "NS.1", "VP.1", "NS.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VS.1", "VS.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VS.3", "NP.5", "VC.5", "VS.4", "NS.3", "VS.5", "VP.5", "NR.3", "NR.4", "VC.6", "NS.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NS.5", "VP.6", "NS.6", "VR.5", "VR.6", "VS.6")
colnames(beeCounts) <- sample(colNames)

group2 = unlist(lapply(colnames(beeCounts), function (x) unlist(strsplit(x, "[.]"))[1]))
y2 <- DGEList(counts=beeCounts, group=group2)
keep2 <- rowSums(cpm(y2)>3) >= 24
y2 <- y2[keep2, , keep.lib.sizes=FALSE]
y2 <- calcNormFactors(y2)
y2 = estimateCommonDisp(y2)
y2 = estimateTagwiseDisp(y2)
de2 = exactTest(y2, pair=c("NR","VP"))
tt2 = topTags(de2, n=nrow(y2))
tt2 <- merge( y2[[1]][,which(group2 %in% c("NR", "VP"))], tt2[[1]], by="row.names")
NRVP_Perm <- tt2

# Can see that NRVP and NRVP_Perm does not change count values in each column, even though p-values/FDR change
# See PermCountSame.xls
# Original 6=VP.1; 19=VP.4; 33=VP.5; 44=VP.6
# Permuted 6=VP.2; 19-NR.5; 33=NR.1; 44=VP.6