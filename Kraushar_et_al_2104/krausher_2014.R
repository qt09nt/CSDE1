#Written by Queenie Tsang

#November 2023

#script for processing the Krausher et al (2014) polysome profiling dataset
#this dataset comes from the supplementary information section of this paper:
# 
# "Kraushar, M. L., Thompson, K., Wijeratne, H. R. S., Viljetic, B., Sakers, K., Marson, J. W., Kontoyiannis, D. L., Buyske, S., Hart, R. P., & Rasin, M. R. (2014). 
#   Temporally defined neocortical translation and polysome assembly are determined by the RNA-binding protein Hu antigen R. Proceedings of the National Academy of Sciences 
#   of the United States of America, 111(36). https://doi.org/10.1073/pnas.1408305111"
# 
# https://www.pnas.org/doi/full/10.1073/pnas.1408305111

# The dataset was filtered with Ribodetector, then aligned with Kallisto pseudoalignment. This part is to do the differential expression
#analysis with Sleuth.

library(sleuth)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(dplyr)

#set directory containing Ribodetector filtered Kallisto files:
setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Kraushar et al 2014 dataset/kallisto/kallisto")

#read in metadata table
krausher_meta <-read.csv("Krausher_ribodetector_kallisto_metadata.csv")

#subset the metadata mono and poly samples
unique(krausher_meta$condition)

krausher_mono_poly_meta <- krausher_meta[krausher_meta$condition %in% c("M_P0_WT", "P_P0_WT"),]


#-------------------------- Defining the Annotation -------------------------------------------
t2g <- read.table("refseq2gene_mouse", colClasses =rep("character", 2), header=TRUE)
t2g <- dplyr::distinct(t2g, target_id, .keep_all= TRUE)

# ------------------------- Protein Coding Genes List from ENSEMBL BIOMART --------------------
protein_coding <- read.csv("protein_coding_genes_biomart_ensembl.txt")
protein_coding$gene_caps <- toupper(protein_coding$Gene.name)

# -------------------------- filtering function ------------------
#filter out any features that do not have at least 5 estimated counts in at least 47
design_filter <- design_filter <- function(design, row, min_reads=5, min_prop = 0.47){
  sum(apply(design, 2, function(x){
    y <- as.factor(x);
    return(max(tapply(row, y, function(f){sum(f >= min_reads)})/
                 tapply(row, y, length)) == 1 
           || basic_filter(row, min_reads, min_prop)
    )
  })) > 0}
# --------------------------------------

##### ----------------------------- Relevel the Condition column with selected reference condition
krausher_mono_poly_meta$condition <- as.factor(krausher_mono_poly_meta$condition)
krausher_mono_poly_meta$condition <- relevel(krausher_mono_poly_meta$condition, ref="M_P0_WT")
#selected reference condition is the WT monosome 

krausher_mono_poly_meta$condition

#####wild type poly vs mono - create sleuth object
so_WTpoly_vs_WTmono <- sleuth_prep(krausher_mono_poly_meta, extra_bootstrap_summary = TRUE, aggregation_column="gene", gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))

#fit the full model
so_WTpoly_vs_WTmono <- sleuth_fit(so_WTpoly_vs_WTmono, ~condition, 'full')

models(so_WTpoly_vs_WTmono)

#Wald Test for differential expression of WT mono versus WT total 
wald_test_WTpoly_vs_WTmono <- sleuth_wt(so_WTpoly_vs_WTmono, which_beta = 'conditionP_P0_WT')

#output wald test results
sleuth_wald_test_WTpoly_vs_WTmono <- sleuth_results(wald_test_WTpoly_vs_WTmono, test = 'conditionP_P0_WT', show_all=TRUE)

#remove genes which are just NA for all samples
sleuth_wald_test_WTpoly_vs_WTmono_filtered <- na.omit(sleuth_wald_test_WTpoly_vs_WTmono)

#filter for protein coding genes only
sleuth_wald_test_WTpoly_vs_WTmono_filtered$target_id <- toupper(sleuth_wald_test_WTpoly_vs_WTmono_filtered$target_id)

sleuth_wald_test_WTpoly_vs_WTmono_protein_coding <- sleuth_wald_test_WTpoly_vs_WTmono_filtered[sleuth_wald_test_WTpoly_vs_WTmono_filtered$target_id %in% protein_coding$gene_caps,]

write.csv(sleuth_wald_test_WTpoly_vs_WTmono_protein_coding, "C:/Users/queenie.tsang/Desktop/CSDE1/Kraushar et al 2014 dataset/results/sleuth_wald_test_WTpoly_vs_WTmono_protein_coding.csv")