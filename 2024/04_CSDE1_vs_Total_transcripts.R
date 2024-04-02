########## Analysis of CSDE1 vs Total RIP in transcripts mode 
######### April 2 2024

library(sleuth)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(dplyr)


##########
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


#############

#subset the metadata table for just CSDE1 and Capture/Total RIP samples 
s2c_CSDE1_total <- s2c[s2c$condition %in% c("Csde", "Capture"),]

#get the kallisto expression matrix in TPM
# set the condition you want as the reference condition, in this case it's IgG
s2c_CSDE1_total<- set_metadata_control_condition(metadata_table = s2c_CSDE1_total, reference = "Capture")

#initialize the sleuth object for comparing CSDE1 with Capture(AKA Input)
so_CSDE1_total <- sleuth_prep(s2c_CSDE1_total, extra_bootstrap_summary = TRUE, transform_fun_counts=function (x) log2(x+0.5))

##fit the full model 
so_CSDE1_total <- sleuth_fit(so_CSDE1_total, ~condition, 'full')

###### Do the sleuth wald test CSDE1 vs Total
models(so_CSDE1_total)


conditionCsde

#wald test for Csde1 versus total differential expression analysis 
wald_test_csde_vs_total <-sleuth_wt(so_CSDE1_total, which_beta = "conditionCsde")

#output results for csde1 vs IgG wald test
sleuth_wald_test_csde_vs_total <- sleuth_results(wald_test_csde_vs_total, test = "conditionCsde", show_all = TRUE)

#remove the genes with NA for all samples 
new_sleuth_wald_test_csde_vs_total <- na.omit(sleuth_wald_test_csde_vs_total)

################# write the CSDE1 vs Total/Input/Capture sleuth results
write.csv(new_sleuth_wald_test_csde_vs_total, "results/sleuth_wald_test_csde_vs_total_all_transcripts.csv")














################## FOR THE CSDE1 DATASET TPM MATRIX, include all 3 conditions (CSDE1, Input and IgG for the normalization)
###################

###
####### get the expression table
CSDE1_total_kallisto_table <- kallisto_table(so_CSDE1_total, use_filtered=TRUE, normalized=TRUE)

######### reformat the expression table so that samples are columns
CSDE1_total_kallisto_table_keep <- CSDE1_total_kallisto_table[,c("target_id", "sample", "tpm")]

CSDE1_total_expression <- CSDE1_total_kallisto_table_keep %>% spread(key = sample, value = tpm)

write.csv(CSDE1_total_expression, "results/CSDE1_total_transcripts_mode_expression_in_TPM.csv")





