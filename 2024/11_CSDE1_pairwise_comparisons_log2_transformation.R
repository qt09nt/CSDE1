### Sleuth RNA seq just do pairwise comparison for generating the TPM values, log 2 transformation

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

############
#subset the s2c metadata table for CSDE1 and IgG 
s2c_CSDE1_IGG <- s2c[s2c$condition %in% c("IgG", "Csde"),]

####### set the IgG condition as the base line/reference condition:
s2c_CSDE1_IGG <- set_metadata_control_condition(metadata_table=s2c_CSDE1_IGG, reference="IgG")

so_CSDE1_IGG <- sleuth_prep(s2c_CSDE1_IGG, extra_bootstrap_summary = TRUE, transform_fun_tpm = function(x) log2(x + 0.5))

so_CSDE1_IGG <- sleuth_fit(so_CSDE1_IGG, ~condition, 'full')
so_CSDE1_IGG <-sleuth_fit(so_CSDE1_IGG, ~1, 'reduced')

####### get the kallisto table expression table from this CSDE1 and IGG sleuth object:
CSDE1_IGG_kallisto_table <- kallisto_table(so_CSDE1_IGG, use_filtered=TRUE, normalized=TRUE)

######## reformat the expression table so that samples are columns
CSDE1_IGG_samples_kallisto_table_keep <- CSDE1_IGG_kallisto_table[,c("target_id", "sample", "tpm")]

CSDE1_IGG_samples_expression_TPM <- CSDE1_IGG_samples_kallisto_table_keep %>% spread(key = sample, value = tpm)

#### log 2 transform the TPM table 
CSDE1_IGG_samples_expression_TPM_copy <- CSDE1_IGG_samples_expression_TPM 
CSDE1_IGG_samples_expression_TPM_copy[,2:7] <- log2(CSDE1_IGG_samples_expression_TPM_copy[,2:7]+0.5)

### 
write.csv(CSDE1_IGG_samples_expression_TPM_copy, "results/April 10 2024/CSDE1_IGG_samples_kallisto_table_expression_TPM_normalized_log2.csv")


#generate CSDE1 vs IgG wald test results

##################### Sleuth wald test
csde1_igg_wald_test <- sleuth_wt(so_CSDE1_IGG, which_beta='conditionCsde')

models(so_CSDE1_IGG)

#output results
sleuth_results_csde1_igg_wald_test <- sleuth_results(csde1_igg_wald_test, test="conditionCsde", show_all = TRUE)

#remove the rows where it's just NA:
sleuth_wald_test_csde1_igg <- na.omit(sleuth_results_csde1_igg_wald_test)

#### 
write.csv(sleuth_wald_test_csde1_igg, "results/April 10 2024/sleuth_wald_test_csde1_vs_igg_log2.csv")

############# 

######### try merging the 2 tables together:
CSDE1_vs_IgG_combined <- merge(CSDE1_IGG_samples_expression_TPM_copy, sleuth_wald_test_csde1_igg, by.x="target_id", by.y="target_id")

write.csv(CSDE1_vs_IgG_combined, "results/April 10 2024/CSDE_vs_IgG_combined_log2.csv")
