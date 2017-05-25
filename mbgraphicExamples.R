library(mbgraphic)
data(Election2013)
dim(Election2013)
iaunivariate(Election2013)

election_num <- Election2013[,sapply(Election2013,is.numeric)] 
dim(election_num)
iacorrgram(election_num) # allows user to zoom in on individual subplots

# clustering based on optimal leaf ordering of variables
vc <- varclust(Election2013,mincor=0.8)
summary(vc)
# the reduced data set
election_reduced <- vc$dfclusrep
dim(election_reduced)

# sdf calculates and measures 9 metrics from scagnostics package
scagdf <- sdf(election_reduced)
summary(scagdf)

# Add custom/additional measures to scagnostics object
addscag <- scag2sdf(election_reduced,scagfun.list=list(dcor2d=dcor2d,splines2d=splines2d))
# merge 'addscag' and 'scagsdf'
scagdf2 <- mergesdfdata(scagdf,addscag)
# merged list contains  11 scagnostics
names(scagdf2$sdf)

# Interactive PCP
iascagpcp(scagdf2)

# Interactive scaggram
iascaggram(scagdf) # has special shiny input box

