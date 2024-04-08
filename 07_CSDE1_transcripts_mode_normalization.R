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
so_all_samples <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, transform_fun_tpm=function(x) log2(x+0.5))

so_all_samples <- sleuth_fit(so_all_samples, ~condition, 'full')

########### get the expression table 

CSDE1_all_samples_kallisto_table <- kallisto_table(so_all_samples, use_filtered=TRUE, normalized=TRUE)

######## reformat the expression table so that samples are columns
CSDE1_all_samples_kallisto_table_keep <- CSDE1_all_samples_kallisto_table[,c("target_id", "sample", "tpm")]

CSDE1_all_samples_expression_TPM <- CSDE1_all_samples_kallisto_table_keep %>% spread(key = sample, value = tpm)

setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP")

write.csv(CSDE1_all_samples_expression_TPM, "results/CSDE1_all_sampes_expression_TPM_normalize_TPM.csv")


#################
#subset the metadata table for just CSDE1 and Capture/Total RIP samples 
s2c_CSDE1_total <- s2c[s2c$condition %in% c("Csde", "Capture"),]

#get the kallisto expression matrix in TPM
# set the condition you want as the reference condition, in this case it's IgG
s2c_CSDE1_total<- set_metadata_control_condition(metadata_table = s2c_CSDE1_total, reference = "Capture")

#initialize the sleuth object for comparing CSDE1 with Capture(AKA Input)
so_CSDE1_total <- sleuth_prep(s2c_CSDE1_total, extra_bootstrap_summary = TRUE, transform_fun_tpm=function (x) log2(x+0.5))

##fit the full model 
so_CSDE1_total <- sleuth_fit(so_CSDE1_total, ~condition, 'full')

###### Do the sleuth wald test CSDE1 vs Total
models(so_CSDE1_total)


#conditionCsde

#wald test for Csde1 versus total differential expression analysis 
wald_test_csde_vs_total <-sleuth_wt(so_CSDE1_total, which_beta = "conditionCsde")

#output results for csde1 vs IgG wald test
sleuth_wald_test_csde_vs_total <- sleuth_results(wald_test_csde_vs_total, test = "conditionCsde", show_all = TRUE)

#remove the genes with NA for all samples 
new_sleuth_wald_test_csde_vs_total <- na.omit(sleuth_wald_test_csde_vs_total)

################# write the CSDE1 vs Total/Input/Capture sleuth results
write.csv(new_sleuth_wald_test_csde_vs_total, "results/sleuth_wald_test_csde_vs_total_all_transcript_new_April_8.csv")


######################## Wald Test for CSDE1 vs IgG:

#subset the metadata table with the conditions being tested:
s2c_CSDE1_IgG <- s2c[,c("sample", "condition", "path")]

#set the condition you want as the reference condition, in this case it's IgG
s2c_CSDE1_IgG<- set_metadata_control_condition(metadata_table = s2c_CSDE1_IgG, reference = "IgG")

#initialize the sleuth object for comparing CSDE1 with IgG
so_CSDE1_IgG <- sleuth_prep(s2c_CSDE1_IgG, extra_bootstrap_summary = TRUE, transform_fun_tpm=function (x) log2(x+0.5))

#fit the full model 
so_CSDE1_IgG <- sleuth_fit(so_CSDE1_IgG, ~condition, 'full')

models(so_CSDE1_IgG)

#wald test for Csde1 versus IgG differential expression analysis 
wald_test_csde_vs_igg <-sleuth_wt(so_CSDE1_IgG, which_beta = "conditionCsde")

#output results for csde1 vs IgG wald test
sleuth_wald_test_csde_vs_igg <- sleuth_results(wald_test_csde_vs_igg, test = "conditionCsde", show_all = TRUE)

#remove the genes with NA for all samples 
new_sleuth_wald_test_csde_vs_igg <- na.omit(sleuth_wald_test_csde_vs_igg)

write.csv(new_sleuth_wald_test_csde_vs_igg, "results/sleuth_wald_test_CSDE1_IgG_newApril8_2024.csv") 


############## merge the tables together
#keep only the b value and q value of CSDE vs IgG table
csde1_vs_total_keep <-new_sleuth_wald_test_csde_vs_total[,c("target_id","qval","b")]
colnames(csde1_vs_total_keep)<-c("target_id","qval_csde1_vs_total","b_csde1_vs_total")
  
csde1_vs_igg_keep <-new_sleuth_wald_test_csde_vs_igg[,c("target_id", "qval", "b")] 
colnames(csde1_vs_igg_keep)<- c("target_id", "qval_csde1_vs_igg", "b_csde1_vs_igg")

combined <-merge(CSDE1_all_samples_expression_TPM, csde1_vs_total_keep, by.x="target_id", by.y="target_id")

combined2 <- merge(combined, csde1_vs_igg_keep, by.x="target_id", by.y="target_id")


write.csv(combined2, "results/CSDE1_TPM_expression_transform_fun_tpm_log2_April_8.csv")

####### log 2 transform the TPM values: and calculate log2FC separately
combined2[,2:10] <- log2(combined2[,2:10] + 0.5) 

write.csv(combined2, "results/CSDE1_dataset_log2_transformed_TPM_April_8_2024_new_4PM.csv")
