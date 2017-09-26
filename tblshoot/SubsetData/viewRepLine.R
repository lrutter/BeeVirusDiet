load("../../data/bees_metrics.rda")
load("../../data/bees_data.rda")

bees_data[,-1] = log(bees_data[,-1]+1)

plotRepLine(bees_data, bees_metrics)
