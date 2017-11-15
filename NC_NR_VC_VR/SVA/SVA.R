thisPath <- getwd()

beeCounts <- readRDS("../data/data.Rds")
e <- as.matrix(beeCounts)

exVars <- read.csv("../data/extraVarClean.csv")
exVars <- exVars[c(1:6,13:18,25:30,37:42),]
colnames(exVars)[1] <- "sample"
exVars$Day <- as.factor(exVars$Day)
exVars$sample <- c(rep("NC",6),rep("NR",6),rep("VC",6),rep("VR",6))
p <- exVars

mod = model.matrix(~as.factor(sample) + as.factor(Lane) + as.factor(Day) + Mortality + SBV + IAPV + rnaConc, data=p)
mod0 = model.matrix(~as.factor(Lane) + as.factor(Day) + Mortality + SBV + IAPV + rnaConc, data=p)

n.sv = num.sv(e,mod,method="leek")
svobj = sva(e,mod,mod0,n.sv=n.sv)

pValues = f.pvalue(e,mod,mod0) #440 <0.05
qValues = p.adjust(pValues,method="BH") #0 <0.05

modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(edata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

###################
library(sva)
library(bladderbatch)
library(pamr)
library(limma)
data(bladderdata)

###### USE SVA FOR UNKNOWN EXTRA VARIABLES ###### 

# Set up data (57 samples; 22,283 genes)
pheno = pData(bladderEset) #ExtraVariables
edata = exprs(bladderEset) #CountTable
mod = model.matrix(~as.factor(cancer), data=pheno)
mod0 = model.matrix(~1,data=pheno)

n.sv = num.sv(edata,mod,method="leek")
svobj = sva(edata,mod,mod0,n.sv=n.sv)

pValues = f.pvalue(edata,mod,mod0)
qValues = p.adjust(pValues,method="BH")

# Now these p/q values account for surrogate vars
modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(edata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

# Now use limma package to get DEGs
fit = lmFit(edata,modSv)
contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)),"C3"=c(-1,0,1,rep(0,svobj$n.sv)))
fitContrasts = contrasts.fit(fit,contrast.matrix)
eb = eBayes(fitContrasts)
topTableF(eb, adjust="BH")

###### USE COMBAT FOR KNOWN BATCH EFFECTS ###### 
