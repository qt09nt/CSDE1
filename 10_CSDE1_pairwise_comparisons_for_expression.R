### Sleuth RNA seq just do pairwise comparison for generating the TPM values, default transformation

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

so_CSDE1_IGG <- sleuth_prep(s2c_CSDE1_IGG, extra_bootstrap_summary = TRUE)


so_CSDE1_IGG <- sleuth_fit(so_CSDE1_IGG, ~condition, 'full')
so_CSDE1_IGG <-sleuth_fit(so_CSDE1_IGG, ~1, 'reduced')

####### get the kallisto table expression table from this CSDE1 and IGG sleuth object:
CSDE1_IGG_kallisto_table <- kallisto_table(so_CSDE1_IGG, use_filtered=TRUE, normalized=TRUE)

######## reformat the expression table so that samples are columns
CSDE1_IGG_samples_kallisto_table_keep <- CSDE1_IGG_kallisto_table[,c("target_id", "sample", "tpm")]

CSDE1_IGG_samples_expression_TPM <- CSDE1_IGG_samples_kallisto_table_keep %>% spread(key = sample, value = tpm)

### 
write.csv(CSDE1_IGG_samples_expression_TPM, "results/April 10 2024/CSDE1_IGG_samples_kallisto_table_expression_TPM_normalized.csv")


### read in the previously generated CSDE1 vs IgG wald test results using default settings

CSDE1_vs_IgG_wald_test<-read.csv("results/April 10 2024/sleuth_wald_test_csde1_vs_igg_default_settings_April_9.csv")

