dat <- readRDS("/Users/lindz/BeeVirusDiet/beeDataFN.rds")
dat <- dat[,c(1, 8:13, 44:49)]
dat[,-1] <- log2(dat[,-1]+1)

minVal = min(dat[,-1])
maxVal = max(dat[,-1])

my_fn <- function(data, mapping, ...){
  x = data[,c(as.character(mapping$x))]
  y = data[,c(as.character(mapping$y))]
  p <- ggplot(data = dat, aes(x=x, y=y)) + coord_cartesian(xlim = c(minVal, maxVal), ylim = c(minVal, maxVal))
  p
}

p <- ggpairs(dat[,-1], lower = list(continuous = my_fn))

ggpairs(dat[,c(2,3)], lower = list(continuous = my_fn))

nCol = ncol(dat)


attach(faithful) 
eruption.lm = lm(eruptions ~ waiting) 
newdata = data.frame(waiting=80)
predict(eruption.lm, newdata, interval="predict")

