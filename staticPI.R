dat <- readRDS("/Users/lindz/BeeVirusDiet/beeDataRLog.rds")
dat <- dat[,c(1, 8:13, 44:49)]

attach(faithful) 
eruption.lm = lm(eruptions ~ waiting) 
newdata = data.frame(waiting=80) 

predict(eruption.lm, newdata, interval="predict") 
nCol = ncol(dat)
