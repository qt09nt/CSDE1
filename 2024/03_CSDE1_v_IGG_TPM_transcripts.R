########## Analysis of CSDE1 vs IgG RIP in transcripts mode 
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
##Keep only the sample, State and path columns of the metadata file
s2c <- s2c[,c("Sample", "State", "path")]

#subset the metadata table for just IgG and CSDE1 RIP samples 
s2c_CSDE1_IgG <- s2c[s2c$State %in% c("Csde", "IgG"),]

#rename the metadata "state" column to "condition"
colnames(s2c_CSDE1_IgG)[2]<-"condition"
colnames(s2c_CSDE1_IgG)[1]<-"sample"

# set the condition you want as the reference condition, in this case it's IgG
s2c_CSDE1_IgG<- set_metadata_control_condition(metadata_table = s2c_CSDE1_IgG, reference = "IgG")

#initialize the sleuth object for comparing CSDE1 with IgG
so_CSDE1_IgG <- sleuth_prep(s2c_CSDE1_IgG, extra_bootstrap_summary = TRUE, transform_fun_counts=function (x) log2(x+0.5))

##fit the full model 
so_CSDE1_IgG <- sleuth_fit(so_CSDE1_IgG, ~condition, 'full')

####### get the expression table
CSDE1_IgG_kallisto_table <- kallisto_table(so_CSDE1_IgG, use_filtered=TRUE, normalized=TRUE)

######### reformat the expression table so that samples are columns
CSDE1_IgG_kallisto_table_keep <- CSDE1_IgG_kallisto_table[,c("target_id", "sample", "tpm")]

CSDE1_IgG_expression <- CSDE1_IgG_kallisto_table_keep %>% spread(key = sample, value = tpm)


####### create a combined table which has both the CSDE1 & IgG samples in TPM, with transcripts and gene names
CSDE1_IgG_transcripts_merged <- merge(CSDE1_IgG_expression, t2g, by.x="target_id", by.y="target_id")

############ save the expression matrix which contains all genes (both protein coding and non protein coding)
write.csv(CSDE1_IgG_transcripts_merged, "results/CSDE1_IgG_transcripts_TPM_merged_all_transcripts.csv")


############ filter for transcripts which match to protein coding genes
##### convert the gene column to capitals to better match the protein_coding table for filtering
CSDE1_IgG_transcripts_merged$gene <- toupper(CSDE1_IgG_transcripts_merged$gene)

####keep just the protein coding genes
CSDE1_IgG_transcripts_merged_protein_coding <- CSDE1_IgG_transcripts_merged[CSDE1_IgG_transcripts_merged$gene %in% protein_coding$gene_caps,]

####save this as a file:
write.csv(CSDE1_IgG_transcripts_merged_protein_coding,"results/CSDE1_IgG_transcripts_TPM_protein_coding.csv")

############## for the non-protein coding genes, see how many transcripts there are and how many genes
CSDE1_IgG_transcripts_merged_non_protein_coding <- CSDE1_IgG_transcripts_merged[!(CSDE1_IgG_transcripts_merged$gene %in% protein_coding$gene_caps),]

### get the number of unique genes in the non-protein coding genes list
length(unique(CSDE1_IgG_transcripts_merged_non_protein_coding$gene))

colnames(s2c)<-c("sample","condition", "path")  

######### do the sleuth wald test for CSDE1 vs IgG
#subset s2c for CSDE1 and IgG
s2c_csde_igg <-s2c[s2c$condition %in% c("IgG", "Csde"),]

##### first relevel the factors so that IgG is the base level 
s2c_csde_igg$condition <- as.factor(s2c_csde_igg$condition)
s2c_csde_igg$condition <- relevel(s2c_csde_igg$condition, ref= "IgG")

so_csde_igg <- sleuth_prep(s2c_csde_igg, extra_bootstrap_summary = TRUE, transform_fun_counts= function(x) log2(x+0.5))
so_csde_igg <- sleuth_fit(so_csde_igg, ~condition, 'full')

models(so_csde_igg) 

#conditionCsde
#wald test for Csde1 versus IgG differential expression analysis 
csde_vs_igg_wt <-sleuth_wt(so_csde_igg, which_beta = "conditionCsde")

#output results for csde1 vs IgG wald test
sleuth_wald_test_csde_vs_IgG <- sleuth_results(csde_vs_igg_wt, test = "conditionCsde", show_all = TRUE)

#remove the genes with NA for all samples 
new_sleuth_wald_test_csde_vs_IgG <- na.omit(sleuth_wald_test_csde_vs_IgG)

### write the Sleuth wald test CSDE1 vs IgG results
write.csv(new_sleuth_wald_test_csde_vs_IgG, "results/sleuth_wald_test_csde_vs_IgG_transcripts_all_genes.csv")
