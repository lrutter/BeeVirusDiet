# This is a filtering and normalization script following edgeR

library(edgeR)
library(ggplot2)
library(GGally)
library(EDASeq)

beeCounts <-read.delim(file="AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)

colnames(beeCounts) <- c("NC.1", "NC.2", "NR.1", "VR.1", "NN.1", "VP.1", "NN.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VN.1", "VN.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VN.3", "NP.5", "VC.5", "VN.4", "NN.3", "VN.5", "VP.5", "NR.3", "NR.4", "VC.6", "NN.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NN.5", "VP.6", "NN.6", "VR.5", "VR.6", "VN.6")

beeCounts <- beeCounts[ , order(names(beeCounts))]

y <- DGEList(counts=beeCounts)

################## Strict edgeR Filtering ################## 
# Smallest library has 3,044,259 reads
minLib <- min(y$samples$lib.size)

# CPM of 3 corresponds to a count of ~9 in the minimum number of samples in a group (24). Usually a gene is required to have a count of 5-10 in a library to be considered expressed in that library.
keep <- rowSums(cpm(y)>3) >= 24
# Number of genes 15,314--> 8,581
y <- y[keep, , keep.lib.sizes=FALSE]

################ edgeR Normalization ################
yERN <- calcNormFactors(y)

################ edgeR Normalization Visuals ##############

# Cannot see much in this boxplot. Many large outliers.
ggparcoord(data.frame(yERN[[1]]), columns=1:48, alphaLines=0, boxplot=TRUE, scale="globalminmax") + coord_flip()

allGroups <- c(rep("NC",6), rep("NN",6), rep("NP",6), rep("NR",6), rep("VC",6), rep("VN",6), rep("VP",6), rep("VR",6))
yERN$samples$group <- allGroups
plotMDS(y, col = c("red","deeppink","darkorange","gold","green2", "green4","blue", "purple")[factor(allGroups)], cex=0.6)

plotMDS(y, col = c("blue","blue","blue","blue","red","red","red","red")[factor(allGroups)], cex=0.6)

plotMDS(y, col = c("white","red","blue","white","white","white","white","white")[factor(allGroups)], cex=0.6)

plotMDS(y, col = c("white","white","white","white","white","red","blue","white")[factor(allGroups)], cex=0.6)

plotMDS(y, col = c("white","white","blue","white","white","red","white","white")[factor(allGroups)], cex=0.6)



######## edgeR GLM approach (3.2.3) ######## 
################################################

group <- c(rep("NC",6), rep("NN",6), rep("NP",6), rep("NR",6), rep("VC",6), rep("VN",6), rep("VP",6), rep("VR",6))
y$samples$group <- c(rep("NC",6), rep("NN",6), rep("NP",6), rep("NR",6), rep("VC",6), rep("VN",6), rep("VP",6), rep("VR",6))

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
fit <- glmFit(y, design)

my.contrasts <- makeContrasts(
  Drug.1vs0 = Drug.1h-Drug.0h,
  Drug.2vs0 = Drug.2h-Drug.0h,
  Placebo.1vs0 = Placebo.1h-Placebo.0h,
  Placebo.2vs0 = Placebo.2h-Placebo.0h,
  DrugvsPlacebo.0h = Drug.0h-Placebo.0h,
  DrugvsPlacebo.1h = (Drug.1h-Drug.0h)-(Placebo.1h-Placebo.0h),
  DrugvsPlacebo.2h = (Drug.2h-Drug.0h)-(Placebo.2h-Placebo.0h),
+ levels=design)


