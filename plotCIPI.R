library(ggplot2) 

dat <- data.frame(qsec=mtcars$qsec, wt=mtcars$wt)
m <- lm(wt ~ qsec, data = dat) 
mpi <- cbind(dat, predict(m, interval = "prediction"))
# Keep only points that are outside the prediction interval
plotPoints <- mpi[which(!(mpi$wt > mpi$lwr & mpi$wt < mpi$upr)),]

# Create prediction interval data frame with upper and lower lines corresponding to sequence covering minimum and maximum of x values in original dataset
minX <- min(mpi$qsec)
maxX <- max(mpi$qsec)
newx <- seq(minX - (maxX-minX), maxX + (maxX-minX), by=0.05)
pred_interval <- predict(m, newdata=data.frame(qsec=newx), interval="prediction", level = 0.95)
pred_interval <- as.data.frame(pred_interval)

pred_interval$qsec = newx
ggplot(data=plotPoints, aes(x = qsec)) + geom_point(aes(y=wt)) + geom_ribbon(data=pred_interval, aes(ymin = lwr, ymax = upr), fill = "blue", alpha = 0.2) + coord_cartesian(xlim = c(minX, maxX))

