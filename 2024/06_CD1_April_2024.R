########## CD1 E17 dataset analysis - transcripts version

######### Queenie Tsang
##### April 3 2024

######## get the expression matrix for the CD1 dataset for set 2 and set 3 only

setwd("C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024")

#use only set 2 and set 3 of total, poly, mono from CD1 dataset
meta<- read.csv("C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/results/CD1_metadata_set2_set3_only.csv")

### read in the ENSEMBL Biomart table


#load in functions for processing data
source("new_Celf2_polysome_seq_functions_feb_2024.R")

library(sleuth)
library(dplyr)
library(ggplot2)
library(readxl)
library(ggrepel)
library(stringr)
library(tidyverse)
library(tidyr)

#-------------------------- Defining the Annotation -------------------------------------------
t2g <- read.table("refseq2gene_mouse.txt", colClasses =rep("character", 2), header=TRUE)
t2g <- dplyr::distinct(t2g, target_id, .keep_all= TRUE)

# ------------------------- Protein Coding Genes List from ENSEMBL BIOMART --------------------
protein_coding <- read.csv("protein_coding_genes_biomart_ensembl.txt")
protein_coding$gene_caps <- toupper(protein_coding$Gene.name)

ENSEMBL <-read.csv("GRCm38_ENSEMBL_ID_Refseq_mRNA_ID_mart_export.txt", sep="\t")

##################### First get the Expression matrix in TPM  for transcripts
#############normalized different conditions separately ie. only mono, only
#############poly, only total

#subset the metadata table for mono samples only
meta_mono <- meta[meta$condition == "Mono",]
so_mono <- sleuth_prep(meta_mono, extra_bootstrap_summary = TRUE, transform_fun_counts=function (x) log2(x+0.5))
CD1_mono_kallisto_table <- kallisto_table(so_mono, use_filtered=TRUE, normalized=TRUE)

#reformat the table, keeping only target_id, sample, and TPM columns
CD1_mono_kallisto_table_keep <- CD1_mono_kallisto_table[,c("target_id", "sample", "tpm")]

#reformat the table so that the samples are columns
CD1_mono_kallisto_table_keep <- CD1_mono_kallisto_table_keep %>% spread(key = sample, value = tpm)

#save the table 
write.csv(CD1_mono_kallisto_table_keep, "results/CD1_mono_kallisto_table_all_transcripts_TPM.csv")

########### do the same thing for Poly samples
meta_poly <- meta[meta$condition == "Poly",]
so_poly <- sleuth_prep(meta_poly, extra_bootstrap_summary = TRUE, transform_fun_counts=function (x) log2(x+0.5))
###### 
CD1_poly_kallisto <- kallisto_table(so_poly, use_filtered=TRUE, normalized=TRUE)
CD1_poly_kallisto_keep <- CD1_poly_kallisto[, c("target_id", "sample", "tpm")]

#reformat this table 
CD1_poly_kallisto_keep <- CD1_poly_kallisto_keep %>% spread(key = sample, value = tpm)

######### save this
write.csv(CD1_poly_kallisto_keep, "results/CD1_poly_expression_table_all_transcripts_TPM.csv")

################# do the same thing for Total
meta_total <- meta[meta$condition == "Total", ]
so_total <- sleuth_prep(meta_total, extra_bootstrap_summary = TRUE, transform_fun_counts=function (x) log2(x+0.5))
CD1_total_kallisto <- kallisto_table(so_total, use_filtered=TRUE, normalized=TRUE)
CD1_total_kallisto_keep <- CD1_total_kallisto[,c("target_id", "sample", "tpm")]
CD1_total_kallisto_keep <- CD1_total_kallisto_keep %>% spread(key = sample, value = tpm)

## save this 
write.csv(CD1_total_kallisto_keep, "results/CD1_total_kallisto_expression_all_transcripts_TPM.csv")

############# do Sleuth Wald test for Poly vs Mono


########### subset the metadata table for Poly and Mono samples
meta_poly_mono <- meta[meta$condition %in% c("Poly", "Mono"),]

###set the condition column as factor 
meta_poly_mono$condition <- as.factor(meta_poly_mono$condition)

meta_poly_mono$condition <- relevel(meta_poly_mono$condition, ref="Mono")

#sleuth wald test in transcripts mode
so_poly_mono <- sleuth_prep(meta_poly_mono, extra_bootstrap_summary = TRUE, transform_fun_counts = function(x) log2(x+0.5))

so_poly_mono <- sleuth_fit(so_poly_mono, ~condition, 'full')
models(so_poly_mono)
#test condition = conditionPoly

wald_test_poly_vs_mono <- sleuth_wt(so_poly_mono, which_beta='conditionPoly')
wald_test_poly_vs_mono_results <- sleuth_results(wald_test_poly_vs_mono, test = 'conditionPoly', show_all = TRUE)

#remove NA rows
wald_test_poly_vs_mono_results <- na.omit(wald_test_poly_vs_mono_results)


#### save results
write.csv(wald_test_poly_vs_mono_results, "results/CD1_wald_test_poly_vs_mono_all_transcripts.csv")


########################### do the same test for Poly vs Total
meta_poly_total <- meta[meta$condition %in% c("Poly", "Total"),]
meta_poly_total$condition <- as.factor(meta_poly_total$condition)
meta_poly_total$condition <- relevel(meta_poly_total$condition, ref="Total")

####### sleuth wald test in transcripts mode 
so_poly_total <- sleuth_prep(meta_poly_total, extra_bootstrap_summary = TRUE, transform_fun_counts = function(x) log2(x+0.5))
so_poly_total <-sleuth_fit(so_poly_total, ~condition, 'full')
models(so_poly_total)
#conditionPoly

wald_test_poly_vs_total <- sleuth_wt(so_poly_total, which_beta='conditionPoly')
wald_test_poly_total_results <- sleuth_results(wald_test_poly_vs_total,  test = 'conditionPoly', show_all = TRUE)

###remove rows with NA:
wald_test_poly_total_results_no_NA <- na.omit(wald_test_poly_total_results)

##save 
write.csv(wald_test_poly_total_results_no_NA, "results/CD1_wald_test_poly_total_results_all_transcripts.csv")

############################### merge the different tables into a single table as requested by Guang


combined<- merge(CD1_poly_kallisto_keep, CD1_mono_kallisto_table_keep, by.x="target_id", by.y="target_id")

combined1 <- merge(combined, CD1_total_kallisto_keep, by.x="target_id", by.y="target_id")

###merge the sleuth results for Poly vs Mono
poly_vs_mono_keep <-wald_test_poly_vs_mono_results[,c("target_id","qval","b")]

#rearrange the columns 
poly_vs_mono_keep <- poly_vs_mono_keep[,c(1,3,2)]

#rename the columns
colnames(poly_vs_mono_keep) <- c("target_id", "poly_vs_mono_b", "poly_vs_mono_qval") 

########## merge these
combined2 <- merge(combined1, poly_vs_mono_keep, by.x="target_id", by.y="target_id")

################## merge the Poly vs Total wald test results

poly_vs_total_keep <- wald_test_poly_total_results_no_NA[,c("target_id", "qval","b")]

########### reorder the columns 
poly_vs_total_keep <- poly_vs_total_keep[,c(1,3,2)]
colnames(poly_vs_total_keep)<-c("target_id", "poly_vs_total_b", "poly_vs_total_qval")


###### merge this with the rest
combined3 <- merge(combined2, poly_vs_total_keep, by.x="target_id", by.y="target_id")

##### add the gene names column
combined4 <- merge(combined3, t2g, by.x="target_id", by.y="target_id")


######## save this table with all transcripts
write.csv(combined4, "results/CD1_all_transcripts.csv")

combined4 <- read.csv("results/CD1_all_transcripts.csv")
  

############# filter for protein coding genes only
############
combined4$gene <- toupper(combined4$gene)

combined_protein_coding <- combined4[combined4$gene %in% protein_coding$gene_caps,]

#save table
write.csv(combined_protein_coding, "results/CD1_transcripts_filtered_for_protein_coding_genes.csv")


############ the decimals at the end of the target_id refers to the version number:

#https://genome.ucsc.edu/FAQ/FAQgenes.html
#https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly

# What is a gene or transcript accession?
#   
# Gene symbols such as BRCA1 are easy to remember but sometimes change and are not specific to an organism. Therefore most databases
# internally use unique identifiers to refer to sequences and some journals require authors to use these in manuscripts.
# 
# The most common accession numbers encountered by users are either from Ensembl, GENCODE or RefSeq. Human Ensembl/GENCODE gene accession numbers
# start with ENSG followed by a number and version number separated by a dot, e.g. "ENSG00000012048.21" for latest BRCA1. Every ENSG-gene has at least one transcript
# assigned to it. The transcript identifiers start with with ENST and are likewise followed by a version number, e.g. "ENST00000619216.1". Additional details on 
# Ensembl IDs can be found on the Ensembl FAQ page.
# 
# NCBI refers to genes with plain numbers, e.g. 672 for BRCA1. Manually curated RefSeq transcript identifiers start with NM_ (coding) or NR_ (non-coding), 
# followed by a number and version number separated by a dot, e.g. "NR_046018.2". If the transcript was predicted by the NCBI Gnomon software, the prefix is XM_ but
# these are rare in human. A table of these and other RefSeq prefixes can be found on the NCBI website. 


######### try removing the decimal ending (version number) for the CD1 table (all transcripts) - contains 
#########  47,654 transcripts originally

########### what to do about the transcripts which have both a "RefSeq.mRNA.ID" and ""RefSeq.mRNA.predicted.ID" in the 
#### ENSEMBL table???  keep just the row with RefSeq mRNA ID?


combined5 <- combined4 
combined5_new <- combined5 %>% separate_wider_delim(target_id, delim=".", names = c("target_id", "version"))

###### try combining the ENSEMBL IDs 
ENSEMBL_keep <- ENSEMBL[,c("Transcript.stable.ID.version", "RefSeq.mRNA.ID","RefSeq.mRNA.predicted.ID","RefSeq.ncRNA.ID")]

merged_by_RefSeq_mRNA_ID <- merge(combined5_new, ENSEMBL_keep, by.x="target_id", by.y="RefSeq.mRNA.ID", keep_all=FALSE)
#34,455 rows

write.csv(merged_by_RefSeq_mRNA_ID, "results/CD1_merged_by_RefSeq_mRNA_ID.csv")

##### for the transcripts in the ENSEMBL table which have both a "RefSeq.mRNA.ID" AND "RefSeq.mRNA.predicted.ID"
#### keep just the transcript rows where "RefSeq.mRNA.ID" is empty AND there is something in the "RefSeq.mRNA.predicted.ID" column
#### which will prevent redundancy

##### subset the ENSEMBL_keep table for only rows where "RefSeq.mRNA.predicted.ID" is NOT empty and "RefSeq.mRNA.ID" IS empty
test_ensembl_mRNA_predicted_ID <- ENSEMBL_keep[ENSEMBL_keep$RefSeq.mRNA.ID == "" & ENSEMBL_keep$RefSeq.mRNA.predicted.ID != "",]

###### merge the CD1 table with ENSEMBL table by the RefSeq.mRNA.predicted.ID
merged_by_RefSeq_predicted_ID <- merge(combined5_new, test_ensembl_mRNA_predicted_ID, by.x="target_id", by.y="RefSeq.mRNA.predicted.ID", keep_all=FALSE)
#3153 rows

write.csv(merged_by_RefSeq_predicted_ID, "results/CD1_merged_by_RefSeq_predicted_ID.csv")
################ 



############ subset the ENSEMBL_keep table to keep only transcripts that have a RefSeq.ncRNA.ID only and DO NOT have a 
#"RefSeq.mRNA.ID" and DO NOT have a "RefSeq.mRNA.predicted.ID"

ENSEMBL_RefSeq_nCRNA_ID <- ENSEMBL_keep[ENSEMBL_keep$RefSeq.mRNA.ID == "",]
ENSEMBL_RefSeq_nCRNA_ID <- ENSEMBL_RefSeq_nCRNA_ID[ENSEMBL_RefSeq_nCRNA_ID$RefSeq.mRNA.predicted.ID == "",]
ENSEMBL_RefSeq_nCRNA_ID <- ENSEMBL_RefSeq_nCRNA_ID[!(ENSEMBL_RefSeq_nCRNA_ID$RefSeq.ncRNA.ID == ""),]

#### merge the combined table with the ENSEMBL_RefSeq_ncRNA_ID
merged_by_RefSeq_ncRNA_ID <- merge(combined5_new, ENSEMBL_RefSeq_nCRNA_ID, by.x="target_id", by.y="RefSeq.ncRNA.ID")
#1374 rows


###### combine the 3 merged tables: merged_by_RefSeq_mRNA_ID, merged_by_RefSeq_predicted_ID, merged_by_RefSeq_ncRNA_ID
###keep only the same columns for each table for merging using rbind()
merged_by_RefSeq_predicted_ID<-merged_by_RefSeq_predicted_ID[,c("target_id", "Transcript.stable.ID.version","CD1.E17.2.Poly_S5", "CD1.E17.3.Poly_S8",  "CD1.E17.2.Mono_S6",           
                                 "CD1.E17.3.Mono_S9", "CD1.E17.2.Total_S4", "CD1.E17.3.Total_S7",          
                                 "poly_vs_mono_b", "poly_vs_mono_qval", "poly_vs_total_b",             
                                 "poly_vs_total_qval", "gene")]

merged_by_RefSeq_mRNA_ID<-merged_by_RefSeq_mRNA_ID[,c("target_id", "Transcript.stable.ID.version", "CD1.E17.2.Poly_S5", "CD1.E17.3.Poly_S8",  "CD1.E17.2.Mono_S6",           
                                                      "CD1.E17.3.Mono_S9", "CD1.E17.2.Total_S4", "CD1.E17.3.Total_S7",          
                                                      "poly_vs_mono_b", "poly_vs_mono_qval", "poly_vs_total_b",             
                                                      "poly_vs_total_qval", "gene")]

merged_by_RefSeq_ncRNA_ID <- merged_by_RefSeq_ncRNA_ID[,c("target_id", "Transcript.stable.ID.version", "CD1.E17.2.Poly_S5", "CD1.E17.3.Poly_S8",  "CD1.E17.2.Mono_S6",           
                                                          "CD1.E17.3.Mono_S9", "CD1.E17.2.Total_S4", "CD1.E17.3.Total_S7",          
                                                          "poly_vs_mono_b", "poly_vs_mono_qval", "poly_vs_total_b",             
                                                          "poly_vs_total_qval", "gene")]

#combine the 3 dataframes together
stacked <-rbind(merged_by_RefSeq_mRNA_ID, merged_by_RefSeq_predicted_ID)

stacked2 <- rbind(stacked, merged_by_RefSeq_ncRNA_ID)

########### it looks the resulting table still has several rows (transcripts) which are duplicated (where the values for all columns)
########### are exactly the same
########### try keeping just the unique rows from the dataframe
stacked3 <-unique(stacked2)

write.csv(stacked3, "results/CD1_E17_cortex_polysome_all_transcripts_ENSEMBL_ID_by_merging_ENSEMBL_Biomart.csv")
#23540 transcripts which correspond to 15245 unique genes 

#### the original dataframe containing all transcripts  contains 47654 rows
combined4_unique <- unique(combined4)

#23540/47654 = 0.4939774

# So using the merging ENSEMBL ID Biomart exported table method, only around 50% of the transcripts are left 
# because only around 50% of the transcripts have both the RefSeq ID and a corresponding ENSEMBL ID


