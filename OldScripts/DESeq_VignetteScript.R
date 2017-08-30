source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("parathyroidSE")
library("DESeq2")
library("parathyroidSE")

data("parathyroidGenesSE")
se <- parathyroidGenesSE
colnames(se) <- se$run
colData(se)

#### Starting from count tables #### 
countdata <- assay(parathyroidGenesSE)
coldata <- colData(parathyroidGenesSE)
rownames(coldata) <- coldata$run
colnames(countdata) <- coldata$run
# We now have all the ingredients to prepare our data object in a form that is suitable for analysis, namely: 1) countdata: a table with the read counts 2) coldata: a table with metadata on the count table’s columns. To now construct the data object from the matrix of counts and the metadata table, we use:
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ patient + treatment)
ddsFull <- ddsFullCountTable

# Take care of technical replicates (libraries derived from same sample). Here, sample SRS308873 is run twice.
as.data.frame(colData(ddsFull)[ ,c("sample","patient","treatment","time")])

# We recommend to first add together technical replicates (i.e., libraries derived from the same samples), such that we have one column per sample. We have implemented a convenience function for this, which can take am object, either SummarizedExperiment or DESeqDataSet, and a grouping factor, in this case the sample name, and return the object with the counts summed up for each unique sample.
ddsCollapsed <- collapseReplicates(ddsFull, groupby = ddsFull$sample, run = ddsFull$run )
head(as.data.frame(colData(ddsCollapsed)[ ,c("sample","runsCollapsed")]), 12)

# We can confirm that the counts for the new object are equal to the summed up counts of the columns that had the same value for the grouping factor
original <- rowSums(counts(ddsFull)[, ddsFull$sample == "SRS308873"])
all(original == counts(ddsCollapsed)[, "SRS308873"])


#####################################################################################
# Running the DESeq pipeline
dds <- ddsCollapsed[, ddsCollapsed$time == "48h"]

# Make sure Control is the first level in the treatment factor, so  the default log2 fold changes are treatment over control
dds$treatment <- relevel(dds$treatment, "Control")
as.data.frame(colData(dds))

# Run the main command
dds <- DESeq(dds)

# Calling results without any arguments extracts the estimated log2 fold changes and p values for the last variable in the design formula.
res <- results(dds)

# In general, the results for a comparison of any two levels of a variable can be extracted using the contrast argument to results. The user should specify three values: the name of the variable, the name of the level in the numerator, and the name of the level in the denominator. Here we extract results for the log2 of the fold change of DPN / Control.
res <- results(dds, contrast = c("treatment", "DPN", "Control"))

# We subset the results table to these genes and then sort it by the log2 fold change estimate to get the significant genes with the strongest down-regulation
resSig <- res[which(res$padj < 0.1 ), ]
head(resSig[order(resSig$log2FoldChange), ])


plotMA(res, ylim = c(-1, 1))
plotDispEsts(dds, ylim = c(1e-6, 1e1))
hist(res$pvalue, breaks=20, col="grey")


# create bins using the quantile function
qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
# "cut" the genes into the bins
bins <- cut( res$baseMean, qs )
# rename the levels of the bins using the middle point
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
# calculate the ratio of £p£ values less than .01 for each bin
ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) ) # plot these ratios
barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")

# Independent filtering (Removal of low-count genes) can result in more p-values. This function tells us the number of rejections for various mean normalized count cut-offs.
attr(res,"filterThreshold")

# Performing rlog transformation for exploratory analysis (visualization)
rld <- rlog(dds)
head(assay(rld))

# Compare first two samples. With rlog on right, there is no longer excessive variation in low read counts, like there was in ordinary logarithm on left.
par(mfrow = c(1, 2))
plot(log2(1+counts(dds, normalized=TRUE)[, 1:2]), col="#00000020", pch=20, cex=0.3)
plot(assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3)

# Examine distance btwn samples and see if it fits experimental design
sampleDists <- dist(t(assay(rld)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$treatment, rld$patient, sep="-")
colnames(sampleDistMatrix) <- NULL
library("gplots")
library("RColorBrewer")
colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
heatmap.2(sampleDistMatrix, trace="none", col=colours)

# Another way to visualize sample-to-sample distances is a principal-components analysis (PCA)
ramp <- 1:3/3
cols <- c(rgb(ramp, 0, 0),
          rgb(0, ramp, 0),
          rgb(0, 0, ramp),
          rgb(ramp, 0, ramp))
print(plotPCA(rld, intgroup = c("patient", "treatment"), col=cols))

# We can also do hiearchical clustering of genes (instead of samples). We usually only do this for the most variable genes. Here, we choose the 35 genes with the most variance across samples.
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))



