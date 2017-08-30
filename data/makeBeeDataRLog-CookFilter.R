---
  title: "Honey Bees - Count Table Pictures"
author: Lindsay Rutter
output:
  packagedocs::package_docs:
  toc: true
toc_collapse: true
vignette: >
  %\VignetteEngine{packagedocs::redirect}
---
  
  <meta http-equiv="content-type" content="text/html;charset=utf-8" />
  
  ```{r global_options, include=FALSE}
# R output pre blocks are styled by default to indicate output
knitr::opts_chunk$set(
  comment = NA,
  cache = TRUE,
  fig.height = 8,
  fig.width = 10
)
```

```{r}
beeCounts <-read.delim(file="../AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
colnames(beeCounts) <- c("NC.1", "NC.2", "NR.1", "VR.1", "NS.1", "VP.1", "NS.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VS.1", "VS.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VS.3", "NP.5", "VC.5", "VS.4", "NS.3", "VS.5", "VP.5", "NR.3", "NR.4", "VC.6", "NS.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NS.5", "VP.6", "NS.6", "VR.5", "VR.6", "VS.6")
countdata <- beeCounts[ , order(names(beeCounts))]
countdata <- as.matrix(countdata)
coldata = data.frame(row.names = colnames(countdata), virus = unlist(lapply(colnames(countdata), function (x) substring(unlist(strsplit(x, "[.]"))[1],1,1))), diet = unlist(lapply(colnames(countdata), function (x) substring(unlist(strsplit(x, "[.]"))[1],2,2))), treatment = unlist(lapply(colnames(countdata), function (x) unlist(strsplit(x, "[.]"))[1])))
dds = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ treatment)
dds <- DESeq(dds)

# 1189 genes (7.7% of total) have all zero counts
length(which(mcols(dds)$allZero))

# Cooks are stored
cookValues = assays(dds)[["cooks"]]
# Each column has 1189 NA values for Cook Values (because all zeros)
summary(cookValues)




rld <- rlog(dds)
dat <- as.data.frame(assay(rld))
setDT(dat, keep.rownames = TRUE)[]
colnames(dat)[1] <- "ID"
dat$ID <- as.character(dat$ID)
dat <- as.data.frame(dat)
saveRDS(dat, "beeDataRLog.rds")
