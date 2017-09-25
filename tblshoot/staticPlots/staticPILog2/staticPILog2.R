library(ggplot2) 

outDir = "/Users/lindz/BeeVirusDiet/staticPILog2"
predLevel = 0.999999
pointSize = 0.1

my_fn <- function(data, mapping, ...){
  xChar = as.character(mapping$x)
  yChar = as.character(mapping$y)
  x = data[,c(xChar)]
  y = data[,c(yChar)]
  m <- lm(y ~ x, data = data) # changed right side to "data" from "dat"
  mpi <- cbind(dat, predict(m, interval = "prediction", level=predLevel))
  # Keep only points that are outside the prediction interval
  plotPoints <- mpi[which(!(mpi[yChar] > mpi$lwr & mpi[yChar] < mpi$upr)),]
  pred_interval <- predict(m, newdata=data.frame(x=newx), interval="prediction", level = predLevel)
  pred_interval <- as.data.frame(pred_interval)
  pred_interval[xChar] = newx
  p <- ggplot(data = plotPoints, aes_string(x = xChar)) + geom_point(aes_string(y = yChar), size=pointSize) + geom_ribbon(data= pred_interval, aes(ymin = lwr, ymax = upr), fill = "blue", alpha = 0.2) + coord_cartesian(xlim = c(minX, maxX), ylim = c(minX, maxX))
  p
}

dat <- readRDS("beeDataFN.rds")
dat[,2:ncol(dat)] <- log2(dat[,2:ncol(dat)]+1)
# Log will cause negative values, so set them to zero
dat[,2:ncol(dat)] <- apply(dat[,2:ncol(dat)], 2, function(x) ifelse(x < 0, 0, x))
dat$ID <- as.character(dat$ID)
colNames <- colnames(dat)
myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
myPairs <- myPairs[-which(myPairs=="ID")]
colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])

# Create prediction interval data frame with upper and lower lines corresponding to sequence covering minimum and maximum of x values in original dataset (don't consider last column because it is not on the x axis in any of the individual scatterplots)
minX <- min(dat[,c(2:(ncol(dat)-1))])
maxX <- max(dat[,c(2:(ncol(dat)-1))])
newx <- seq(minX - (maxX-minX), maxX + (maxX-minX), by=0.05)

for (i in 1:(length(myPairs)-1)){
  for (j in (i+1):length(myPairs)){
    group1 = myPairs[i]
    group2 = myPairs[j]
    datSel <- cbind(ID=dat$ID, dat[,which(colGroups %in% c(group1, group2))])
    maxVal = max(abs(datSel[,-1]))
    maxRange = c(0, maxVal)
    #xbins=10
    #buffer = maxRange[2]/xbins
    p <- ggpairs(datSel[,-1], lower = list(continuous = my_fn))
    jpeg(filename=paste0(outDir, "/", group1, "_", group2, ".jpg"), height=700, width=700)
    print(p)
    dev.off()
  }
}
