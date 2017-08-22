results <- readRDS("RemSampleComp.rds")

tResults <- t(results[,-1])
colnames(tResults) <- results[,1]
tResults <- as.data.frame(tResults)
setDT(tResults, keep.rownames = TRUE)[]
colnames(tResults)[1] <- "Pair"
tResults$ID <- as.factor(tResults$Pair)
tResults <- as.data.frame(tResults)

boxDat <- tResults %>% gather(key, val, -c(Pair))
boxDat$val <- as.numeric(boxDat$val)
p <- ggplot(boxDat, aes(key, val)) + geom_boxplot()
ggplotly(p)

#works
getOutliers <- function(input){
  rowInfo = list()
  for (i in 1:28){
    row = unlist(input[i,])
    hiThresh <- median(row) + 1.5*IQR(row)
    loThresh <- median(row) - 1.5*IQR(row)
    hiOutlier <- names(which(row > hiThresh))[-1]
    loOutlier <- names(which(row < loThresh))[-1]
    rowInfo[paste0(colnames(input)[i+1], "_loOutlier")] <-list(loOutlier)
    rowInfo[paste0(colnames(input)[i+1], "_hiOutlier")] <-list(hiOutlier)
  }
  rowInfo
}



getOutliers <- function(input){
  rowInfo = list()
  for (i in 1:28){
    row = unlist(input[i,])
    hiThresh <- median(row) + 1.5*IQR(row)
    loThresh <- median(row) - 1.5*IQR(row)
    hiOutlier <- names(sort(row[2:49][which(row[2:49] > hiThresh)]))
    loOutlier <- names(sort(row[2:49][which(row[2:49] < loThresh)]))
    rowInfo[paste0(input[i,1], "_loOutlier")] <-list(loOutlier)
    rowInfo[paste0(input[i,1], "_hiOutlier")] <-list(hiOutlier)
  }
  rowInfo
}

remOutliers <- getOutliers(results)
