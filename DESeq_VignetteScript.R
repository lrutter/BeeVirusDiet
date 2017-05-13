source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("parathyroidSE")
library("DESeq2")
library("parathyroidSE")

data("parathyroidGenesSE")
se <- parathyroidGenesSE

