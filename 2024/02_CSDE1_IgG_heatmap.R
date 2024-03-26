#### CSDE1 vs IgG heatmap
#### March 26 2024
### Queenie Tsang

setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/results/")

#read in the expression values in TPM for the selected genes
CSDE1_IgG<- read.csv("CSDE1_IgG_expression_genes_set2.csv")

#get rid of the first column which is just row number from original matrix file
CSDE1_IgG$X <- NULL

#https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
Heatmap()

