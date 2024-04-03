
########## Analysis of CSDE1 RIP in transcripts mode - Expression table in TPM
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


################## FOR THE EXAMPLE TPM MATRIX, include all 3 conditions (CSDE1, Input and IgG for the normalization)

s2c <- s2c[,colnames(s2c) %in% c("Sample", "State", "path")]
colnames(s2c)<- c("sample", "condition", "path")

so_all_samples <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, transform_fun_counts=function (x) log2(x+0.5))

##fit the full model 
so_all_samples <- sleuth_fit(so_all_samples, ~condition, 'full')

CSDE1_all_samples_kallisto_table <- kallisto_table(so_all_samples, use_filtered=TRUE, normalized=TRUE)

############ reformat the expression table so that samples are columns
CSDE1_all_samples_kallisto_table_keep <- CSDE1_all_samples_kallisto_table[,c("target_id", "sample", "tpm")]

CSDE1_all_samples_expression <- CSDE1_all_samples_kallisto_table_keep %>% spread(key = sample, value = tpm)

write.csv(CSDE1_all_samples_expression, "results/CSDE1_all_samples_expression_transcripts_TPM.csv")

########## combine what we have so far for the merged table:

setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/results/")

sleuth_wald_test_csde_vs_IgG <- read.csv("sleuth_wald_test_csde_vs_IgG_transcripts_all_genes.csv")

##keep only target_id, b value and qvalue for CSDE1 vs IGG
sleuth_wald_test_csde_vs_IgG_keep <- sleuth_wald_test_csde_vs_IgG[,c("target_id","qval","b")]
colnames(sleuth_wald_test_csde_vs_IgG_keep)<- c("target_id", "qval_csde1_vs_igg","b_csde1_vs_igg")

sleuth_wald_test_csde_vs_total<-read.csv("sleuth_wald_test_csde_vs_total_all_transcripts.csv")
sleuth_wald_test_csde_vs_total_keep <- sleuth_wald_test_csde_vs_total[,c("target_id","qval","b")]
colnames(sleuth_wald_test_csde_vs_total_keep)<- c("target_id","qval_csde1_vs_total", "b_csde1_vs_total")

merge1<- merge(CSDE1_all_samples_expression, sleuth_wald_test_csde_vs_IgG_keep, by.x="target_id", by.y="target_id")

#reorder columns so bvalue CSDE1 vs IGG b value comes before qvalue CSDE1 vs IgG: 
merge2<- merge1[,c(1:10, 12, 11)]

#do the same for b value CSDE1 vs total: 
sleuth_wald_test_csde_vs_total_keep2 <- sleuth_wald_test_csde_vs_total_keep[, c(1,3,2)]

### merge the rest of the dataframe with the sleuth wald test CSDE1 vs total df:
merge3 <- merge(merge2, sleuth_wald_test_csde_vs_total_keep2)


############### Add the gene column as the last column
t2g$gene <- toupper(t2g$gene)

merge4 <- merge(merge3, t2g, by.x="target_id", by.y="target_id")

### save this result
write.csv(merge4, "CSDE1_dataset_all_transcripts.csv")

CSDE1_all_transcripts<-read.csv("results/CSDE1_dataset_all_transcripts.csv")

#filter for just protein coding genes 
CSDE1_protein_coding <-CSDE1_all_transcripts[CSDE1_all_transcripts$gene %in% protein_coding$gene_caps,]

## save
write.csv(CSDE1_protein_coding, "CSDE1_protein_coding_transcripts.csv")

############## try merging the ENSEMBL transcript ID into the table:

ENSEMBL_df <-read.csv("results/GRCm38_ENSEMBL_ID_Refseq_mRNA_ID_mart_export.txt")


######## problem with merging the ENSEMBL transcript IDs is that the RefSeq IDs in the CSDE1 target_id column
#include a decimal number at the end ie.NM_001001178.1, whereas the ENSEMBL Biomart table, the RefSeq mRNA IDs and
# RefSeq mRNA predicted IDs do NOT contain a decimal at the end

#### can try merging the ENSEMBL ids but likely won't work:
test1 <- merge(ENSEMBL_df, CSDE1_all_transcripts, by.x = "RefSeq.mRNA.ID", by.y="target_id")
#0 observations

test2 <- merge(ENSEMBL_df, CSDE1_all_transcripts, by.x = "RefSeq.mRNA.predicted.ID", by.y="target_id")
### 0 observations

# It looks like it may not be possible to merge the ENSEMBL transcript IDS via the method of merging ENSEMBL Biomart table 
# with the CSDE1 table, because in the CSDE1 table, the target_id/accession number ends in a decimal (ie. NM_001001178.1) but 
# the RefSeq mRNA ID, and RefSeq mRNA predicted IDs from the ENSEMBL Biomart exported table DON't contain decimals at the end 
# (ie. NM_001272030)

### may need to re-do via Sleuth and use the ENSEMBL ID directly