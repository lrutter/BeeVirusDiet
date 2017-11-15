require(EDASeq)
require(yeastRNASeq)
require(leeBamViews)

# Unaligned reads
files <- list.files(file.path(system.file(package = "yeastRNASeq"), "reads"), pattern = "fastq", full.names = TRUE)
names(files) <- gsub("\\.fastq.*", "", basename(files))
met <- DataFrame(conditions=c(rep("mut",2), rep("wt",2)), row.names=names(files))
fastq <- FastqFileList(files)
elementMetadata(fastq) <- met

# Aligned reads
files <- list.files(file.path(system.file(package = "leeBamViews"), "bam"), pattern = "bam$", full.names = TRUE)
names(files) <- gsub("\\.bam", "", basename(files))
gt <- gsub(".*/", "", files)
gt <- gsub("_.*", "", gt)
lane <- gsub(".*(.)$", "\\1", gt)
geno <- gsub(".$", "", gt)
pd <- DataFrame(geno=geno, lane=lane, row.names=paste(geno,lane,sep="."))
bfs <- BamFileList(files)
elementMetadata(bfs) <- pd

# Read level EDA
colors <- c(rep(rgb(1,0,0,alpha=0.7),2), rep(rgb(0,0,1,alpha=0.7),2), rep(rgb(0,1,0,alpha=0.7),2), rep(rgb(0,1,1,alpha=0.7),2))
barplot(bfs,las=2,col=colors)

plotQuality(bfs,col=colors,lty=1)
legend("topright",unique(elementMetadata(bfs)[,1]), fill=unique(colors))

plotQuality(bfs[[1]],cex.axis=.8)
barplot(bfs[[1]],las=2)
plotNtFrequency(bfs[[1]])

# Gene-level EDA
