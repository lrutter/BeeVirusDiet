library(ggplot2)
library(hexbin)

set.seed(1)
bDat <- data.frame(Group1 = rnorm(100,0,1), Group2 = rnorm(100,0,1))
points <- data.frame(Group1 = rnorm(10,0.5,1), Group2 = rnorm(10,0.5,1))
maxVal = max(max(bDat), max(points))
minVal = min(min(bDat), min(points))
maxRange = c(minVal, maxVal)
xbins=10
buffer = maxRange[2]/xbins
xChar = "Group1"
yChar = "Group2"
x = bDat[,c(xChar)]
y = bDat[,c(yChar)]
h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
attr(hexdf, "cID") <- h@cID

ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(maxRange[1], maxRange[2]), ylim = c(maxRange[1], maxRange[2])) + geom_point(data = points, aes(x=Group1, y=Group2), inherit.aes = FALSE, color = "orange", size = 0.1)

