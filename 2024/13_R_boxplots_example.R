###### R boxplots example

# Convert the variable dose from a numeric to a factor variable
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
head(ToothGrowth)

#     len supp dose
# 1  4.2   VC  0.5
# 2 11.5   VC  0.5
# 3  7.3   VC  0.5
# 4  5.8   VC  0.5
# 5  6.4   VC  0.5
# 6 10.0   VC  0.5

library(ggplot2)
# Basic box plot
p <- ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_boxplot()
p

# Box plot with jittered points
# 0.2 : degree of jitter in x direction
p + geom_jitter(shape=16, position=position_jitter(0.2))

# Change box plot line colors by groups
p<-ggplot(ToothGrowth, aes(x=dose, y=len, color=dose)) +
  geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.2))
p
