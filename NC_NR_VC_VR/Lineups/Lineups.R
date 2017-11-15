library(edgeR)
library(ggplot2)
library(GGally)
library(EDASeq)
library(utils)

thisPath <- getwd()

beeCounts <- readRDS("../data/data.Rds")
beeCounts <- cbind(ID=rownames(beeCounts), beeCounts)
beeCounts$ID <- as.character(beeCounts$ID)

plotPermutations(beeCounts, nPerm = 20, topThresh = 50, outDir = getwd())
