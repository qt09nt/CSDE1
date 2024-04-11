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
s2c_CSDE1_total<- s2c[s2c$condition %in% c("Csde", "Capture"), ]

####### set the total condition as the base line/reference condition:
s2c_CSDE1_total <- set_metadata_control_condition(metadata_table=s2c_CSDE1_total, reference="Capture")

so_CSDE1_total <- sleuth_prep(s2c_CSDE1_total, extra_bootstrap_summary = TRUE, transform_fun_tpm = function(x) log2(x + 0.5))

so_CSDE1_total <- sleuth_fit(so_CSDE1_total, ~condition, 'full')
so_CSDE1_total <-sleuth_fit(so_CSDE1_total, ~1, 'reduced')


####### get the kallisto table expression table from this CSDE1 and IGG sleuth object:
CSDE1_total_kallisto_table <- kallisto_table(so_CSDE1_total, use_filtered=TRUE, normalized=TRUE)

####### keep only the relevant columns
CSDE1_total_kallisto_table_keep <- CSDE1_total_kallisto_table[,c("target_id", "sample", "tpm")]

########## reformat the table so columns are the samples and transcripts are the rows
CSDE1_total_expression_TPM <- CSDE1_total_kallisto_table_keep %>% spread(key = sample, value = tpm)

CSDE1_total_expression_TPM_copy <-CSDE1_total_expression_TPM
 
CSDE1_total_expression_TPM_copy[,2:7] <- log2(CSDE1_total_expression_TPM_copy[,2:7]+0.5)

################### save the results
write.csv(CSDE1_total_expression_TPM_copy, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/results/April 10 2024/CSDE1_total_expression_TPM_log2.csv")


######### get the CSDE total kallisto table for est counts (normalized):
CSDE1_total_kallisto_table_est_counts <- CSDE1_total_kallisto_table[,c("target_id", "sample", "est_counts")]

CSDE1_total_expression_est_counts<-CSDE1_total_kallisto_table_est_counts %>% spread(key = sample, value= est_counts)
CSDE1_total_expression_est_counts_copy <- CSDE1_total_expression_est_counts
CSDE1_total_expression_est_counts_copy[,2:7]<-log2(CSDE1_total_expression_est_counts_copy[,2:7]+0.5)

write.csv(CSDE1_total_expression_est_counts_copy, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/results/April 10 2024/CSDE1_total_expression_est_counts_normalized_log2.csv")

################################### check the raw data for est_counts
CSDE1_total_kallisto_table_est_counts_raw <- kallisto_table(so_CSDE1_total, use_filtered=TRUE, normalized=FALSE)
CSDE1_total_kallisto_table_est_counts_raw_keep <-CSDE1_total_kallisto_table_est_counts_raw[,c("target_id", "sample", "est_counts")]

###reformat the table so columns are samples and transcripts are rows
CSDE1_total_kallisto_est_counts_raw <- CSDE1_total_kallisto_table_est_counts_raw_keep %>% spread(key = sample, value= est_counts)
CSDE1_total_kallisto_est_counts_raw_copy <-CSDE1_total_kallisto_est_counts_raw

#log2 transform these 
CSDE1_total_kallisto_est_counts_raw_copy[,2:7]<-log2(CSDE1_total_kallisto_est_counts_raw_copy[,2:7]+0.5)

write.csv(CSDE1_total_kallisto_est_counts_raw_copy, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/results/April 10 2024/CSDE1_total_kallisto_est_counts_raw_log2.csv")


######### do the sleuth wald test for differential expression analysis 
models(so_CSDE1_total)

csde1_total_wald_test <- sleuth_wt(so_CSDE1_total, which_beta='conditionCsde')

#output results
sleuth_results_csde1_total_wald_test <- sleuth_results(csde1_total_wald_test, test="conditionCsde", show_all = TRUE)

#remove the rows where it's just NA:
sleuth_wald_test_csde1_total <- na.omit(sleuth_results_csde1_total_wald_test)


### filter for just transcripts with q value < 0.05
sleuth_wald_test_csde1_total_qval_0.05 <- dplyr::filter(sleuth_wald_test_csde1_total, qval <= 0.05)

#### 
write.csv(sleuth_wald_test_csde1_total, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/results/April 10 2024/sleuth_wald_test_csde1_vs_total_log2.csv")

########## merge the TPM table and the Sleuth wald test results together
combined_CSDE1_total <- merge(CSDE1_total_expression_TPM_copy, sleuth_wald_test_csde1_total, by.x="target_id", by.y="target_id")

combined_CSDE1_total <- combined_CSDE1_total[,1:10]
combined_CSDE1_total$pval <- NULL

write.csv(combined_CSDE1_total, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/results/CSDE_total_TPM_sleuth_wald_test_combined_log2.csv")

#######################
plot_group_density(so_CSDE1_total, use_filtered = TRUE, units = "est_counts", 
                   trans ="log2", grouping = setdiff(colnames(so_CSDE1_total$sample_to_covariates),
                                                    "sample"), offset =0.5)

################################# Quality Control Checking Normalization

########## check the distribution of the raw counts data versus normalized counts 

#Visualize the genes expression of CSDE1 and total samples


# library
library(ggplot2)
library(dplyr)
library(hrbrthemes)

# Build dataset with different distributions
data <- data.frame(
  type = c( rep("variable 1", 1000), rep("variable 2", 1000) ),
  value = c( rnorm(1000), rnorm(1000, mean=4) )
)

# Represent it
p <- data %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(fill="")

plot(p)

#### plot log2 est counts for one of the most highly differentially expressd transcripts 
plot_bootstrap(so_CSDE1_total, "XM_006531241.3", units="est_counts", color_by="condition")


library(ggplot2)

######### make a distribution plot of the Sleuth wald test b values CSDE2 vs total: 
ggplot(sleuth_wald_test_csde1_total, aes(x=b)) +
  geom_histogram(binwidth=1, fill="#404080", color="#e9ecef", alpha=0.9) + 
  theme_ipsum() +
  ggtitle("Distribution plot of Sleuth Wald Test - CSDE1 vs Total b values")+
  theme_ipsum()+
  theme(
    plot.title = element_text(size = 15)
  )


################## boxplots of raw counts vs normalized counts:
########## try plotting all the transcripts raw counts:
p <- ggplot(CSDE1_total_kallisto_table_est_counts_raw_keep, aes(x=sample, y=est_counts)) + 
  geom_boxplot()
p


##plot one of the transcripts (random)
est_counts_raw_XM_011243456.2 <- CSDE1_total_kallisto_est_counts_raw_copy[CSDE1_total_kallisto_est_counts_raw_copy$target_id == "XM_011243456.2", ]

est_counts_raw_XM_011243456.2_t <- t(est_counts_raw_XM_011243456.2)
sample_name <-c("Csde1_1",   "Csde1_2",   "Csde1_3",   "Input1",    "Input2",    "Input3")
est_counts_log2 <-c(3.589267, 5.653282, 5.039427, 9.21577, 9.036536, 8.726252)
condition_est_counts_raw <-c("Csde1",   "Csde1",   "Csde1",   "Input",    "Input",    "Input")
  
est_counts_raw_XM_011243456.2 <- data.frame(sample_name, est_counts_log2, condition_est_counts_raw)                   

p <- ggplot(est_counts_raw_XM_011243456.2, aes(x=sample_name, y=est_counts_log2, colour=condition_est_counts_raw)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2))
p
              

################# check the same transcript for the normalized log2 est_counts values:
est_counts_normalized_XM_011243456.2 <- CSDE1_total_expression_est_counts_copy[CSDE1_total_expression_est_counts_copy$target_id == "XM_011243456.2",]
est_counts_normalized_XM_011243456.2_t <- t(est_counts_normalized_XM_011243456.2)

est_counts_log2_normalized <- c(4.58307, 5.157492, 6.248255, 8.616663, 8.541677, 8.078101)
sample_est_counts_normalized <- c("Csde1_1", "Csde1_2", "Csde1_3", "Input1", "Input2", "Input3")
condition_est_counts_normalized <- c("Csde", "Csde", "Csde", "Input", "Input", "Input")

est_counts_normalized_XM_011243456.2 <- data.frame(est_counts_log2_normalized, sample_est_counts_normalized, condition_est_counts_normalized)

p <- ggplot(est_counts_normalized_XM_011243456.2, aes(x=sample_name, y=est_counts_log2, colour=condition_est_counts_normalized)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2))
p



############







#https://pachterlab.github.io/sleuth/docs/plot_volcano.html
plot_volcano

