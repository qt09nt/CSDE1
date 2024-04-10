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

###subset for CSDE1 and total samples
s2c_CSDE1_total 

####### set the total condition as the base line/reference condition:
s2c_CSDE1_total <- set_metadata_control_condition(metadata_table=s2c_CSDE1_total, reference="Capture")

so_CSDE1_total <- sleuth_prep(s2c_CSDE1_total, extra_bootstrap_summary = TRUE, transform_fun_tpm = function(x) log2(x + 0.5))

so_CSDE1_total <- sleuth_fit(so_CSDE1_total, ~condition, 'full')
so_CSDE1_total <-sleuth_fit(so_CSDE1_total, ~1, 'reduced')