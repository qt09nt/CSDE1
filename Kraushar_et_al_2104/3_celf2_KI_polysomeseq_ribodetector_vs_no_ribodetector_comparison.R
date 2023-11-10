#comparison of Ribodetector filtering normalization versus no Ribodetector normalization:
#Celf2-KI-Polysome-seq
#Nov 8 2023

#for the non-Ribodetector filtering pipeline: Kallisto-Sleuth-filtering for protein coding genes only with ENSEMBL-Biomart list
#for the Ribodetector filtering pipeline: Ribodetector-Kallisto-Sleuth-filtering for protein coding genes only with ENSEMBL-Biomart list

library("sleuth")
library(tidyr)
library(ggplot2)

#set directory to kallisto files with no Ribodetector filtering:

setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/kallisto/kallisto_output/kallisto_output")

#read in the metadata file for no Ribodetector filtering:
metadata_no_ribodetector <-read.csv("metadata_kallisto_Celf2_KI_polysomeseq.csv", header = T)

#subset for WT-mono and WT-total
metadata_no_ribodetector_WTmono_WTtotal <- metadata_no_ribodetector[metadata_no_ribodetector$condition %in% c("WT-mono", "WT-total"),]

#subset for WT-mono and WT-poly 
meta_no_ribodetector_WTmono_WTpoly <- metadata_no_ribodetector[metadata_no_ribodetector$condition %in% c("WT-mono", "WT-poly"),]

#subset for WT-poly and WT-total
metadata_no_ribodetector_WTpoly_WTtotal <- metadata_no_ribodetector[metadata_no_ribodetector$condition %in% c("WT-poly", "WT-total"),]

#-------------------------- Defining the Annotation -------------------------------------------
t2g <- read.table("refseq2gene_mouse", colClasses =rep("character", 2), header=TRUE)
t2g <- dplyr::distinct(t2g, target_id, .keep_all= TRUE)

# ------------------------- Protein Coding Genes List from ENSEMBL BIOMART --------------------
setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq")

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

############################ Relevel the Condition column with selected reference condition
metadata_no_ribodetector_WTmono_WTtotal$condition <- as.factor(metadata_no_ribodetector_WTmono_WTtotal$condition)
metadata_no_ribodetector_WTmono_WTtotal$condition <- relevel(metadata_no_ribodetector_WTmono_WTtotal$condition, ref="WT-total")

#WT-mono vs WT-poly: relevel the conditions to set "WT-mono" as the reference
meta_no_ribodetector_WTmono_WTpoly$condition <- as.factor(meta_no_ribodetector_WTmono_WTpoly$condition)
meta_no_ribodetector_WTmono_WTpoly$condition <- relevel(meta_no_ribodetector_WTmono_WTpoly$condition, ref="WT-mono")


### WTpoly vs WTtotal
metadata_no_ribodetector_WTpoly_WTtotal$condition <- as.factor(metadata_no_ribodetector_WTpoly_WTtotal$condition)
metadata_no_ribodetector_WTpoly_WTtotal$condition <- relevel(metadata_no_ribodetector_WTpoly_WTtotal$condition, ref= "WT-total")

#####wild type mono vs total - create sleuth object
so_WT_mono_total_no_ribodetector <- sleuth_prep(metadata_no_ribodetector_WTmono_WTtotal, extra_bootstrap_summary = TRUE, aggregation_column="gene", gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))

# Wt-mono vs WT-poly - create sleuth object
so_WTmono_WTpoly_no_ribodetector <- sleuth_prep(meta_no_ribodetector_WTmono_WTpoly, extra_bootstrap_summary = TRUE, aggregation_column="gene", gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))

# WTpoly vs WT-total create sleuth object 
so_WTpoly_WTtotal_no_ribodetector <- sleuth_prep(metadata_no_ribodetector_WTpoly_WTtotal, extra_bootstrap_summary = TRUE, aggregation_column="gene", gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))

#fit the full model
so_WT_mono_total_no_ribodetector <- sleuth_fit(so_WT_mono_total_no_ribodetector, ~condition, 'full')
so_WTmono_WTpoly_no_ribodetector <- sleuth_fit(so_WTmono_WTpoly_no_ribodetector, ~condition, 'full')
so_WTpoly_WTtotal_no_ribodetector <- sleuth_fit(so_WTpoly_WTtotal_no_ribodetector, ~condition, 'full')

#check the models that have been fit
models(so_WT_mono_total_no_ribodetector)
models(so_WTmono_WTpoly_no_ribodetector)
models(so_WTpoly_WTtotal_no_ribodetector)

#Wald Test for differential expression of WT mono versus WT total - no Ribodetector filtered
wald_test_WT_mono_vs_WT_total_no_ribodetector <- sleuth_wt(so_WT_mono_total_no_ribodetector, which_beta = 'conditionWT-mono')
wald_test_WTmono_vs_WTpoly_no_ribodetector <- sleuth_wt(so_WTmono_WTpoly_no_ribodetector, which_beta = 'conditionWT-poly')
wald_test_WTpoly_vs_WTtotal_no_ribodetector<- sleuth_wt(so_WTpoly_WTtotal_no_ribodetector, which_beta = 'conditionWT-poly')

#output the results 
sleuth_wald_test_WT_mono_vs_WT_total_no_ribodetector <- sleuth_results(wald_test_WT_mono_vs_WT_total_no_ribodetector, test = 'conditionWT-mono',
                                                       show_all = TRUE)
sleuth_wald_test_WTmono_vs_WTpoly_no_ribodetector <- sleuth_results(wald_test_WTmono_vs_WTpoly_no_ribodetector, test = 'conditionWT-poly', show_all=TRUE)

sleuth_wald_test_WTpoly_vs_WTtotal_no_ribodetector <- sleuth_results(wald_test_WTpoly_vs_WTtotal_no_ribodetector, test = 'conditionWT-poly', show_all=TRUE)



#### remove the genes which are just NA for all samples                                                                    
sleuth_wald_test_WT_mono_vs_WT_total_no_ribodetector <- na.omit(sleuth_wald_test_WT_mono_vs_WT_total_no_ribodetector)
sleuth_wald_test_WTmono_vs_WTpoly_no_ribodetector <- na.omit(sleuth_wald_test_WTmono_vs_WTpoly_no_ribodetector)
sleuth_wald_test_WTpoly_vs_WTtotal_no_ribodetector <- na.omit(sleuth_wald_test_WTpoly_vs_WTtotal_no_ribodetector)

#########################  #filter for protein coding genes only
sleuth_wald_test_WT_mono_vs_WT_total_no_ribodetector$target_id <- toupper(sleuth_wald_test_WT_mono_vs_WT_total_no_ribodetector$target_id)
sleuth_wald_test_WT_mono_vs_WT_total_protein_coding_no_ribodetector <- sleuth_wald_test_WT_mono_vs_WT_total_no_ribodetector[sleuth_wald_test_WT_mono_vs_WT_total_no_ribodetector$target_id %in% protein_coding$gene_caps,]

#WTmono vs WTpoly no ribodetector - filter for protein coding genes only 
sleuth_wald_test_WTmono_vs_WTpoly_no_ribodetector$target_id <- toupper(sleuth_wald_test_WTmono_vs_WTpoly_no_ribodetector$target_id)
sleuth_wald_test_WTmono_vs_WTpoly_no_ribodetector_protein_coding <- sleuth_wald_test_WTmono_vs_WTpoly_no_ribodetector[sleuth_wald_test_WTmono_vs_WTpoly_no_ribodetector$target_id %in% protein_coding$gene_caps,]

#WTpoly vs WTtotal no ribodetector - filter for protein coding genes only
sleuth_wald_test_WTpoly_vs_WTtotal_no_ribodetector$target_id <- toupper(sleuth_wald_test_WTpoly_vs_WTtotal_no_ribodetector$target_id)
sleuth_wald_test_WTpoly_vs_WTtotal_no_ribodetector_protein_coding <- sleuth_wald_test_WTpoly_vs_WTtotal_no_ribodetector[sleuth_wald_test_WTpoly_vs_WTtotal_no_ribodetector$target_id %in% protein_coding$gene_caps,] 

#remove some of the ribosomal proteins manually when gene names contain "RPL"
sleuth_wald_test_WTpoly_WTtotal_no_ribodetector_protein_coding_remove_ribosomal_proteins <- sleuth_wald_test_WTpoly_vs_WTtotal_no_ribodetector_protein_coding[!grepl("RPL", sleuth_wald_test_WTpoly_vs_WTtotal_no_ribodetector_protein_coding$target_id),]
df[!grepl("REVERSE", df$Name),]

#save results:
write.csv(sleuth_wald_test_WT_mono_vs_WT_total_protein_coding_no_ribodetector, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/sleuth_wald_test_WT_mono_vs_WT_total_no_ribodetector_protein_coding.csv")
write.csv(sleuth_wald_test_WTmono_vs_WTpoly_no_ribodetector_protein_coding, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/sleuth_wald_test_WTmono_vs_WTpoly_no_ribodetector_protein_coding.csv")
write.csv(sleuth_wald_test_WTpoly_WTtotal_no_ribodetector_protein_coding_remove_ribosomal_proteins, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/sleuth_wald_test_WTpoly_WTtotal_no_ribodetector_protein_coding_remove_ribosomal_proteins.csv")


#looking at the Kallisto table containing the normalized counts
so_WT_mono_total_no_ribodetector_expression_matrix <- kallisto_table(so_WT_mono_total_no_ribodetector, normalized = TRUE)

#re-format the kallisto table to keep just select columns
WTmono_WTtotal_expression <- so_WT_mono_total_no_ribodetector_expression_matrix[,c("sample", "target_id", "scaled_reads_per_base")]
WTmono_WTtotal_expression_reformatted <-spread(WTmono_WTtotal_expression, sample, scaled_reads_per_base)

#filter the Kallisto table expression matrix for protein coding genes only
WTmono_WTtotal_expression_reformatted$target_id <- toupper(WTmono_WTtotal_expression_reformatted$target_id)
WTmono_WTtotal_expression_protein_coding <- WTmono_WTtotal_expression_reformatted[WTmono_WTtotal_expression_reformatted$target_id %in% protein_coding$gene_caps,]

#remove rows where expression is 0 across all samples
ans = WTmono_WTtotal_expression_protein_coding[rowSums(WTmono_WTtotal_expression_protein_coding[,2:7])>0,] 

ans2 <- ans
ans2[,2:7] <-log2(ans2[,2:7] + 0.5)

WTmono_WTtotal_expression_protein_coding_log2_no_ribodetector <-ans2

write.csv(WTmono_WTtotal_expression_protein_coding_log2_no_ribodetector, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/WTmono_WTtotal_expression_no_ribodetector_protein_coding_log2.csv")

### look at the expression matrix for Celf2-KI-Polysomeseq for WTmono vs WTtotal - list of genes after filtering for protein coding genes for 
# the pipeline where Ribodetector filtering has been applied and compare it to the list of genes remaining without Ribodetector 
setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector")

WTmono_WTtotal_ribodetector <- read.csv("WTmono_WTtotal_norm_counts_table_log2_protein_coding.csv")

#look at the genes which were filtered out by Ribodetector:
filtered_out_ribodetector <- WTmono_WTtotal_expression_protein_coding_log2_no_ribodetector[!(WTmono_WTtotal_expression_protein_coding_log2_no_ribodetector$target_id %in% WTmono_WTtotal_ribodetector$target_id), ]

#when you actually search up the genes which were filtered out by Ribodetector, some of these are actually protein
#coding genes
write.csv(filtered_out_ribodetector, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/filtered_out_genes_ribodetector_Celf2-KI-polysomeseq_WTmono_WTtotal.csv")


#### look at the genes filtered out by Ribodetector in the WT-mono vs WT-total results
ribodetector_sleuth_wald_test_WT_mono_vs_WT_total_protein_coding <-read.csv("sleuth_wald_test_WT_mono_vs_WT_total_protein_coding.csv")

sleuth_wald_test_WTmono_vs_WTtotal_genes_filtered_by_ribodetector <- sleuth_wald_test_WT_mono_vs_WT_total_protein_coding_no_ribodetector[!(sleuth_wald_test_WT_mono_vs_WT_total_protein_coding_no_ribodetector$target_id %in% ribodetector_sleuth_wald_test_WT_mono_vs_WT_total_protein_coding$target_id),]



#################### plot scatterplot of With Ribodetector (x-axis) and Without Ribodetector(y-axis)

############# load in the wald test results for WT-mono vs WT-total which have been Ribodetector filtered:

setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/")
wald_test_WTmono_WTtotal_ribodetector <-read.csv("sleuth_wald_test_WT_mono_vs_WT_total_protein_coding.csv", header = T)

#merge the Ribodetector filtered and non-Ribodetector filtered dataframes together before plotting
merged_WTmono_WT_total_ribodetector_vs_no_ribodetector<- merge(wald_test_WTmono_WTtotal_ribodetector, sleuth_wald_test_WT_mono_vs_WT_total_protein_coding_no_ribodetector, by.x = "target_id", by.y = "target_id")

#x-axis is with Ribodetector; y-axis is WITHOUT Ribodetector 
ggplot(merged_WTmono_WT_total_ribodetector_vs_no_ribodetector, aes(x = b.x, y = b.y)) +
  geom_point() +
  theme_bw() +
  ylab("Wald test b value  Celf2-KI Polysome seq WT-mono vs WT-total without Ribodetector") +
  xlab("Wald test b value Celf2-KI Polysome seq WT-mono vs WT-total with Ribodetector")

#### plot the adjusted q value; axis is with Ribodetector; y-axis is WITHOUT Ribodetector
ggplot(merged_WTmono_WT_total_ribodetector_vs_no_ribodetector, aes(x = qval.x, y = qval.y)) +
  geom_point() +
  theme_bw() +
  ylab("Wald test qval value  Celf2-KI Polysome seq WT-mono vs WT-total without Ribodetector") +
  xlab("Wald test qval value Celf2-KI Polysome seq WT-mono vs WT-total with Ribodetector")

########################## 

############################################################# Krausher et al (2014) dataset 

#read in the Kallisto pseudoalignment files (this one is Ribodetector filtered):
setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Kraushar et al 2014 dataset/kallisto/kallisto")

#read in the metadata file:
krausher_metadata <- read.csv("Krausher_ribodetector_kallisto_metadata.csv")  

#subset for just monosome and polysome
krausher_WTmono_WTpoly <- krausher_metadata[krausher_metadata$condition %in% c("M_P0_WT","P_P0_WT"),]



