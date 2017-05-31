library(ggplot2)
library(plotly)
library(tidyr)

set.seed(1)
dat <- data.frame(Name = paste0("Name",1:100), A=abs(rnorm(100)), B=abs(rnorm(100)), C=abs(rnorm(100)), D=abs(rnorm(100)), E=abs(rnorm(100)))
dat$Name <- as.character(dat$Name)
nVar = ncol(dat)
boxDat <- dat[, c(1:nVar)] %>% gather(key, val, -c(Name))
colnames(boxDat) <- c("Name", "Group", "Count")

BP <- ggplot(boxDat, aes(x = Group, y = Count)) + geom_boxplot()
ggBP <- ggplotly(BP)
