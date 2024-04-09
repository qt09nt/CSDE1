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

CSDE1_all_samples_kallisto_table <- kallisto_table(so_all_samples, use_filtered=TRUE, normalized=TRUE)

######## reformat the expression table so that samples are columns
CSDE1_all_samples_kallisto_table_keep <- CSDE1_all_samples_kallisto_table[,c("target_id", "sample", "tpm")]

CSDE1_all_samples_expression_TPM <- CSDE1_all_samples_kallisto_table_keep %>% spread(key = sample, value = tpm)

setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP")

write.csv(CSDE1_all_samples_expression_TPM, "results/CSDE1_all_samples_expression_TPM_sleuth_prep_default_settings.csv")

###########





########## fit the sleuth reduced model
so_all_samples <-sleuth_fit(so_all_samples, ~1, 'reduced')

####### results of the test can be examined with
so_all_samples <- sleuth_lrt(so_all_samples, 'reduced', 'full')

# examine the models that have been fit with the models() function
models(so_all_samples)

sleuth_table <- sleuth_results(so_all_samples, 'reduced:full', 'lrt', show_all = FALSE)

sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

head(sleuth_significant, 20)

######## Test significant differences between conditions using Wald Test 
#https://hbctraining.github.io/DGE_workshop_salmon/lessons/09_sleuth.html

#subset the s2c metadata table for CSDE1 and IgG 
s2c_CSDE1_IGG <- s2c[s2c$condition %in% c("IgG", "Csde"),]

####### set the IgG condition as the base line/reference condition:
s2c_CSDE1_IGG <- set_metadata_control_condition(metadata_table=s2c_CSDE1_IGG, reference="IgG")

so_CSDE1_IGG <- sleuth_prep(s2c_CSDE1_IGG, extra_bootstrap_summary = TRUE)

so_CSDE1_IGG <- sleuth_fit(so_CSDE1_IGG, ~condition, 'full')
so_CSDE1_IGG <-sleuth_fit(so_CSDE1_IGG, ~1, 'reduced')
so_CSDE1_IGG <- sleuth_lrt(so_CSDE1_IGG, 'reduced', 'full')

models(so_CSDE1_IGG)
#conditionCsde   is the test condition


######## results of the test can be examine with
sleuth_table <- sleuth_results(so_CSDE1_IGG, 'reduced:full', 'lrt', show_all=FALSE)

##### conditionIgG
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

write.csv(sleuth_table, "results/CSDE1_vs_IGG_LRT_table.csv")

###################### Sleuth wald test
csde1_igg_wald_test <- sleuth_wt(so_CSDE1_IGG, which_beta='conditionCsde')

#output results
sleuth_results_csde1_igg_wald_test <- sleuth_results(csde1_igg_wald_test, test="conditionCsde", show_all = TRUE)

#remove the rows where it's just NA:
sleuth_wald_test_csde1_igg <- na.omit(sleuth_results_csde1_igg_wald_test)

#### 
write.csv(sleuth_wald_test_csde1_igg, "results/sleuth_wald_test_csde1_vs_igg_default_settings_April_9.csv")


######################### Sleuth wald test
s2c_CSDE1_total <- s2c[s2c$condition %in% c("Csde", "Capture"),]

####### set the IgG condition as the base line/reference condition:
s2c_CSDE1_total <- set_metadata_control_condition(metadata_table=s2c_CSDE1_total, reference="Capture")
so_CSDE1_total <- sleuth_prep(s2c_CSDE1_total, extra_bootstrap_summary = TRUE)

#fit the full model
so_CSDE1_total <- sleuth_fit(so_CSDE1_total, ~condition, 'full')
so_CSDE1_total <- sleuth_fit(so_CSDE1_total, ~condition, 'reduced')

sleuth_wald_test_CSDE1_total <- sleuth_wt(so_CSDE1_total, which_beta='conditionCsde')
so_CSDE1_total <- sleuth_lrt(so_CSDE1_total, 'reduced', 'full')


sleuth_live(so_CSDE1_total)

############ try plotting the distributions using Sleuth plot_group_density()

plot_group_density(so_CSDE1_total, use_filtered = TRUE, units = "est_counts", 
                   trans ="log", grouping = setdiff(colnames(so_CSDE1_total$sample_to_covariates),
                                                    "sample"), offset = 1)

plot_group_density(so_CSDE1_IGG, use_filtered = TRUE, units = "est_counts", 
                   trans ="log", grouping = setdiff(colnames(so_CSDE1_IGG$sample_to_covariates),
                                                    "sample"), offset = 1)

### look at distribution for log 2 transformation
plot_group_density(so_CSDE1_total, use_filtered = TRUE, units = "est_counts", 
                   trans ="log2", grouping = setdiff(colnames(so_CSDE1_total$sample_to_covariates),
                                                    "sample"), offset = 1)

### look at distribution for log 2 transformation
plot_group_density(so_CSDE1_IGG, use_filtered = TRUE, units = "est_counts", 
                   trans ="log2", grouping = setdiff(colnames(so_CSDE1_IGG$sample_to_covariates),
                                                    "sample"), offset = 1)

######## check the distribution for TPM values 
plot_group_density(so_CSDE1_IGG, use_filtered = TRUE, units = "tpm", 
trans ="log", grouping = setdiff(colnames(so_CSDE1_total$sample_to_covariates), "sample"), offset = 1)

plot_group_density(so_CSDE1_total, use_filtered = TRUE, units = "tpm", 
trans ="log", grouping = setdiff(colnames(so_CSDE1_total$sample_to_covariates), "sample"), offset = 1)

plot_group_density(so_CSDE1_IGG, use_filtered = TRUE, units = "tpm", 
                   trans ="log2", grouping = setdiff(colnames(so_CSDE1_total$sample_to_covariates), "sample"), offset = 1)

plot_group_density(so_CSDE1_total, use_filtered = TRUE, units = "tpm", 
                   trans ="log2", grouping = setdiff(colnames(so_CSDE1_total$sample_to_covariates), "sample"), offset = 1)


############# double check the normalization of feature counts data by plotting the log transformed 
############ normalized est counts of CSDE1 and IgG, and CSDE1 and Capture samples
#### make boxplots of the genes - expression in TPM - should be normalized already during the sleuth_prep() step

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(dplyr)

#### pick a random transcript
#NM_001001130.2


CSDE1_NM_001001130.2 <- CSDE1_all_samples_expression_TPM[CSDE1_all_samples_expression_TPM$target_id == "NM_001001130.2",]

CSDE1_NM_001001130.2_t <- t(CSDE1_NM_001001130.2)
CSDE1_NM_001001130.2_t <-data.frame(CSDE1_NM_001001130.2_t)

CSDE1_NM_001001130.2_t$sample <- row.names(CSDE1_NM_001001130.2_t)

colnames(CSDE1_NM_001001130.2_t)[1] <- "normalized_TPM_value"

#rename this table to the name of transcript NM_001001130.2
CSDE1_NM_001001130.2_t <- CSDE1_NM_001001130.2_t[2:10,]

CSDE1_NM_001001130.2_t %>% 
  ggplot(aes(x=sample, y=normalized_TPM_value, fill=sample))+
  geom_boxplot() +
  scale_fill_viridis(discrete=TRUE, alpha = 0.6)+
  geom_jitter(color="black", size=0.4, alpha = 0.9)+
  theme_ipsum() +
  theme(
    plot.title = element_text(size = 11)
  ) +
  ggtitle("CSDE1 RIP dataset - TPM values with default settings - NM_001001130.2 transcript") +
  xlab("sample")
