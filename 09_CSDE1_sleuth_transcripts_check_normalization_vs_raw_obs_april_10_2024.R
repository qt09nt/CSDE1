####### CSDE1 RIP dataset - try running the Sleuth tests with 3 samples for creating the Sleuth Object 
####### and try normalization step using: 

#norm_fun_tpm: a function to perform between sample normalization on the TPM. The default is the DESeq method.
#See norm_factors for details.

library(tidyr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(sleuth)

###read in functions 
setwd("C:/Users/queenie.tsang/Desktop/CSDE1/")
source("C:/Users/queenie.tsang/Desktop/CSDE1/new_Celf2_polysome_seq_functions_feb_2024.R")

#read in the metadata file 
setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/")

s2c <- read.csv("CSDE1_RIP_metadata.csv")

########## #-------------------------- Defining the Annotation -------------------------------------------
t2g <- read.table("refseq2gene_mouse.txt", colClasses =rep("character", 2), header=TRUE)
t2g <- dplyr::distinct(t2g, target_id, .keep_all= TRUE)

# ------------------------- Protein Coding Genes List from ENSEMBL BIOMART --------------------
protein_coding <- read.csv("protein_coding_genes_biomart_ensembl.txt") 
protein_coding$gene_caps <- toupper(protein_coding$Gene.name)

############

s2c <-s2c[,c("Sample", "State", "path")]
colnames(s2c) <- c("sample", "condition", "path")

######## use transform_fun_tpm() parameter since we want the expression values in TPM
######## norm_fun_tpm: a function to perform between sample normalization on the TPM. The default is the DESeq method.
#### keep the default settings for the Sleuth prep 
so_all_samples <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

so_all_samples <- sleuth_fit(so_all_samples, ~condition, 'full')

########### get the expression table 
#normalized = TRUE, return the normalized data
#normalized = FALSE, returns the raw data
CSDE1_all_samples_kallisto_table <- kallisto_table(so_all_samples, use_filtered=TRUE, normalized=TRUE)

######## reformat the expression table so that samples are columns
CSDE1_all_samples_kallisto_table_keep <- CSDE1_all_samples_kallisto_table[,c("target_id", "sample", "tpm")]

CSDE1_all_samples_expression_TPM <- CSDE1_all_samples_kallisto_table_keep %>% spread(key = sample, value = tpm)

write.csv(CSDE1_all_samples_expression_TPM, "results/CSDE1_all_samples_expression_TPM_sleuth_prep_default_settings_normalized_true.csv")


############################## try NOT normalizing and generate the kallisto table 
CSDE1_all_samples_expression_raw_data_kallisto_table <- kallisto_table(so_all_samples, use_filtered=TRUE, normalized = FALSE)
CSDE1_all_samples_expression_raw_data_keep <- CSDE1_all_samples_expression_raw_data[,c("target_id", "sample", "tpm")]

############ 
CSDE1_all_samples_expression_raw_data <- CSDE1_all_samples_expression_raw_data_keep  %>% spread(key = sample, value = tpm )
write.csv(CSDE1_all_samples_expression_raw_data, "results/CSDE1_all_samples_expression_TPM_sleuth_prep_default_settings_raw_data_normalize_false.csv")


################# try out the sleuth matrix obs_norm() vs obs_raw()
#https://rdrr.io/github/pachterlab/sleuth/man/sleuth_to_matrix.html

sleuth_matrix_obs_norm_tpm <- sleuth_to_matrix(so_all_samples, 'obs_norm', 'tpm')
######## this table matches exactly the table output using kallisto_table(normalized=TRUE)

############### 
write.csv(sleuth_matrix_obs_norm_tpm, "CSDE1_all_samples_sleuth_matrix_obs_norm_tpm.csv")
sleuth_matrix_obs_raw_tpm <- sleuth_to_matrix(so_all_samples, 'obs_raw', 'tpm')

###### 
write.csv(sleuth_matrix_obs_raw_tpm, "results/April 10 2024/CSDE1_all_samples_sleuth_matrix_obs_raw_tpm.csv")

########### This sleuth matrix obs_raw TPM is exactly the same as the one generated from Kallisto_table(normalized=FALSE)
##########


