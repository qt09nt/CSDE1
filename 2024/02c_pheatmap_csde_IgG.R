##CSDE1 vs IgG heatmap

#https://biostatsquid.com/step-by-step-heatmap-tutorial-with-pheatmap/
#March 26 2024


setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/results/")

#read in the expression values in TPM for the selected genes
CSDE1_IgG<- read.csv("CSDE1_IgG_expression_genes_set2.csv")

#get rid of the first column which is just row number from original matrix file
CSDE1_IgG$X <- NULL

row.names(CSDE1_IgG)<-CSDE1_IgG$target_id

#keep just the sample columns
CSDE1_IgG_keep <- CSDE1_IgG[,2:7]

###create an annotation dataframe to annotate sample type
condition <- c("CSDE1", "CSDE1", "CSDE1", "IgG", "IgG", "IgG")
