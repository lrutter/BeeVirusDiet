set.seed(1)
dat <- data.frame(ID = paste0("ID",1:100), A=rnorm(100,0,1), B=rnorm(100,0,1), C=rnorm(100,0,1))
subDat <- dat[sample(1:100, 5),]

boxDat <- dat[, c(1:4)] %>% gather(key, val, -c(ID))
colnames(boxDat) <- c("ID", "Sample", "Count")
pcpDat <- subDat[, c(1:4)] %>% gather(key, val, -c(ID))
colnames(pcpDat) <- c("ID", "Sample", "Count")

ggplot(boxDat, aes(x = Sample, y = Count)) + geom_boxplot()
ggplot(pcpDat, aes(x = Sample, y = Count, group = ID)) + geom_line(size = 0.1)

