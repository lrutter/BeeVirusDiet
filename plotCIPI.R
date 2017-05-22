library(ggplot2) 

dat <- readRDS("/Users/lindz/BeeVirusDiet/beeDataRLog.rds")
dat <- dat[,c(1, 8:11)]
predLevel = 0.999999

# Create prediction interval data frame with upper and lower lines corresponding to sequence covering minimum and maximum of x values in original dataset (don't consider last column because it is not on the x axis in any of the individual scatterplots)
minX <- min(dat[,c(2:(ncol(dat)-1))])
maxX <- max(dat[,c(2:(ncol(dat)-1))])
newx <- seq(minX - (maxX-minX), maxX + (maxX-minX), by=0.05)

#y=wt; x=qsec

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
  p <- ggplot(data = plotPoints, aes_string(x = xChar)) + geom_point(aes_string(y = yChar), size=0.5) + geom_ribbon(data= pred_interval, aes(ymin = lwr, ymax = upr), fill = "blue", alpha = 0.2) + coord_cartesian(xlim = c(minX, maxX))
  p
}

#[,-1]
p <- ggpairs(dat[,c(2:4)], lower = list(continuous = my_fn))

