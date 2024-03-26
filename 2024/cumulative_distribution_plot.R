####### make cumulative distribution plot 

setwd("C:/Users/queenie.tsang/Desktop/CSDE1/")

cumulative_curve <- read.csv("cumulative_curve.csv")

plot(ecdf(cumulative_curve))
