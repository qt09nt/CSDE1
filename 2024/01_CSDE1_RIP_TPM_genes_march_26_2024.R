### CSDE1 and IgG RIP in TPM
### March 26 2024
### Queenie Tsang

library(tidyr)
library(tidyverse)
library(ggplot2)
library(dplyr)

###read in functions 
setwd("C:/Users/queenie.tsang/Desktop/CSDE1/")
source("C:/Users/queenie.tsang/Desktop/CSDE1/new_Celf2_polysome_seq_functions_feb_2024.R")

#read in the metadata file 
setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/")

s2c <- read.csv("CSDE1_RIP_metadata.csv")

setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_output/")

setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP")

#-------------------------- Defining the Annotation -------------------------------------------
t2g <- read.table("refseq2gene_mouse.txt", colClasses =rep("character", 2), header=TRUE)
t2g <- dplyr::distinct(t2g, target_id, .keep_all= TRUE)

# ------------------------- Protein Coding Genes List from ENSEMBL BIOMART --------------------
protein_coding <- read.csv("protein_coding_genes_biomart_ensembl.txt")
protein_coding$gene_caps <- toupper(protein_coding$Gene.name)

##################### 



##Keep only the sample, State and path columns of the metadata file
s2c <- s2c[,c("Sample", "State", "path")]

#subset the metadata table for just IgG and CSDE1 RIP samples 
s2c_CSDE1_IgG <- s2c[s2c$State %in% c("Csde", "IgG"),]

#rename the metadata "state" column to "condition"
colnames(s2c_CSDE1_IgG)[2]<-"condition"
colnames(s2c_CSDE1_IgG)[1]<-"sample"

# set the condition you want as the reference condition, in this case it's IgG
s2c_CSDE1_IgG<- set_metadata_control_condition(metadata_table = s2c_CSDE1_IgG, reference = "IgG")

#initialize the sleuth object for comparing CSDE1 with IgG
so_CSDE1_IgG <- sleuth_prep(s2c_CSDE1_IgG, aggregation_column = "gene", gene_mode= TRUE, extra_bootstrap_summary = TRUE, target_mapping = t2g, transform_fun_counts=function (x) log2(x+0.5))

##fit the full model 
so_CSDE1_IgG <- sleuth_fit(so_CSDE1_IgG, ~condition, 'full')

####### get the expression table
CSDE1_IgG_kallisto_table <- kallisto_table(so_CSDE1_IgG, use_filtered=TRUE, normalized=TRUE)

######### reformat the expression table so that samples are columns
CSDE1_IgG_kallisto_table_keep <- CSDE1_IgG_kallisto_table[,c("target_id", "sample", "tpm")]

CSDE1_IgG_expression <- CSDE1_IgG_kallisto_table_keep %>% spread(key = sample, value = tpm)

#save the expression matrix 
write.csv(CSDE1_IgG_expression, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/results/CSDE1_IgG_all_genes_expression_TPM.csv")


CSDE1_IgG_expression$target_id <- toupper(CSDE1_IgG_expression$target_id) 
  
#### keep just the genes requested: Pabpc1, Vim, FABP7, Ybx1
CSDE1_IgG_expression_genes <- CSDE1_IgG_expression[CSDE1_IgG_expression$target_id %in% c("PABPC1", "VIM", "FABP7", "YBX1"),]

#### save as a csv file 
write.csv(CSDE1_IgG_expression_genes, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/results/CSDE1_IgG_expression_genes.csv")

###### get the expression values in TPM for the second set of genes requested:

# BCL11B
# CRIM1
# CUX2
# FOXO1
# GABRA5
# LHX2
# NR4A2
# OMA1
# PCDH17
# Pou3f1
# RELN
# RPRM
# SATB2
# SOX5
# TBR1
# TLE3
# UNC5D

#### keep just the second set of genes requested: 
CSDE1_IgG_expression_genes_set2 <- CSDE1_IgG_expression[CSDE1_IgG_expression$target_id %in% c("BCL11B","CRIM1", "CUX2", "FOXO1",
                                                                                         "GABRA5",
                                                                                         "LHX2",
                                                                                         "NR4A2",
                                                                                         "OMA1",
                                                                                         "PCDH17",
                                                                                         "POU3F1",
                                                                                         "RELN",
                                                                                         "RPRM",
                                                                                         "SATB2",
                                                                                         "SOX5",
                                                                                         "TBR1",
                                                                                         "TLE3",
                                                                                         "UNC5D"),]

### save to csv
write.csv(CSDE1_IgG_expression_genes_set2, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/results/CSDE1_IgG_expression_genes_set2.csv")


### read in the TPM layer genes file
TMP_layer_genes <-read.csv("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/results/TMP_layer.csv")

TMP_layer_genes <- na.omit(TMP_layer_genes)

TMP_layer_genes_copy <- TMP_layer_genes
TMP_layer_genes_copy$symbol <- toupper(TMP_layer_genes_copy$symbol)

#keep just the TMP layer genes from the expression matrix of all genes
CSDE1_IgG_expression_copy <- CSDE1_IgG_expression
CSDE1_IgG_expression_copy$target_id <- toupper(CSDE1_IgG_expression_copy$target_id)

CSDE1_IgG_expression_TMP_layer <-CSDE1_IgG_expression_copy[CSDE1_IgG_expression_copy$target_id %in% TMP_layer_genes_copy$symbol, ]

#combine the original TMP_layer table with the TPM expression 
merged_new <- merge(TMP_layer_genes_copy, CSDE1_IgG_expression_TMP_layer, by.x="symbol", by.y = "target_id", all = FALSE)

write.csv(merged_new, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/results/TMP_layer_genes_TPM_merged_new.csv")
