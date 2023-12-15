### Celf2-KI-Polysomeseq 
### Wild Type versus Celf2-KI
#November 16 2023
#Queenie Tsang

library(sleuth)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(stringr)

#read in functions 
setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/R")
source("Cefl2_polysome_seq_functions.R")

#Ribodetector filtered Kallisto files: 
setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/kallisto/ribodetector/kallisto_ribodetector/")

#rename the samples in metadata file to correct for Celf2-Het-poly-2 and Celf2-Het-mono-2 label switch up (manually done in excel)
#from previous analysis (confirmed with Taylor who performed the experiment)

#metadata for Ribodetector-Kallisto output:
metadata <- read.csv("Celf2-KI-polysome_seq_ribodetector_metadata.csv")

### for WT-total versus KI-total
#subset the metadata table for WT-total and Celf2-KI-total
metadata_wttotal_vs_KItotal <- metadata[metadata$condition %in% c("Het-total", "WT-total"),]
metadata

meta_wtmono_vs_KImono <- metadata[metadata$condition %in% c("Het-mono", "WT-mono"),]

#KI poly vs KI total
metadata_KIpoly_KItotal <- metadata[metadata$condition %in% c("Het-poly", "Het-total"),]

## KI poly vs KI mono
#meta_KI-poly vs KI-mono
metadata_KIpoly_KImono <- metadata[metadata$condition %in% c("Het-mono", "Het-poly"),]

## WT poly vs WT mono
metadata_WTpoly_WTmono <- metadata[metadata$condition %in% c("WT-mono", "WT-poly"),]

## WTpoly vs WTtotal
metadata_WTpoly_WTtotal <- metadata[metadata$condition %in% c("WT-poly", "WT-total"),]

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

###-------------------- ------------------------------ initialize the sleuth object 

### set the reference condition 
metadata_wttotal_vs_KItotal$condition <- as.factor(metadata_wttotal_vs_KItotal$condition)
metadata_wttotal_vs_KItotal$condition <- relevel(metadata_wttotal_vs_KItotal$condition, ref = "WT-total")

meta_wtmono_vs_KImono$condition <- as.factor(meta_wtmono_vs_KImono$condition) 
meta_wtmono_vs_KImono$condition <- relevel(meta_wtmono_vs_KImono$condition, ref = "WT-mono")

metadata_KIpoly_KItotal$condition <- as.factor(metadata_KIpoly_KItotal$condition)
metadata_KIpoly_KItotal$condition <- relevel(metadata_KIpoly_KItotal$condition, ref = "Het-total")

metadata_KIpoly_KImono$condition <- as.factor(metadata_KIpoly_KImono$condition)
metadata_KIpoly_KImono$condition<-relevel(metadata_KIpoly_KImono$condition, ref= "Het-mono")

metadata_WTpoly_WTmono$condition <- as.factor(metadata_WTpoly_WTmono$condition)
metadata_WTpoly_WTmono$condition <- relevel(metadata_WTpoly_WTmono$condition, ref= "WT-mono")

metadata_WTpoly_WTtotal$condition <- as.factor(metadata_WTpoly_WTtotal$condition)
metadata_WTpoly_WTtotal$condition <- relevel(metadata_WTpoly_WTtotal$condition, ref = "WT-total")

metadata_WTpoly_KIpoly$condition <- as.factor(metadata_WTpoly_KIpoly$condition)
metadata_WTpoly_KIpoly$condition <- relevel(metadata_WTpoly_KIpoly$condition, ref = "WT-poly")

#initialize object for Sleuth object 
so_WTtotal_KItotal <- sleuth_prep(metadata_wttotal_vs_KItotal, extra_bootstrap_summary = TRUE, aggregation_column="gene",
                                  gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))
so_WTtotal_KItotal <- sleuth_fit(so_WTtotal_KItotal, ~condition, 'full')

so_all_samples <- sleuth_prep(metadata, extra_bootstrap_summary = TRUE, aggregation_column="gene",
                              gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))

so_WTpoly_KIpoly <- sleuth_prep(metadata_WTpoly_KIpoly, extra_bootstrap_summary = TRUE, aggregation_column = "gene",
                                gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))
so_WTpoly_KIpoly <- sleuth_fit(so_WTpoly_KIpoly, ~condition, 'full')

so_WTmono_KImono <- sleuth_prep(meta_wtmono_vs_KImono, extra_bootstrap_summary = TRUE, aggregation_column="gene",
                                gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))

so_WTmono_KImono <- sleuth_fit(so_WTmono_KImono, ~condition, 'full')


##### KI-poly vs KI-total
so_KIpoly_KItotal <- sleuth_prep(metadata_KIpoly_KItotal,  extra_bootstrap_summary = TRUE, aggregation_column="gene",
                                 gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))
so_KIpoly_KItotal <- sleuth_fit(so_KIpoly_KItotal, ~condition, 'full')


###### KI poly vs KI mono
so_KIpoly_KImono <- sleuth_prep(metadata_KIpoly_KImono,  extra_bootstrap_summary = TRUE, aggregation_column="gene",
                                gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))
so_KIpoly_KImono <- sleuth_fit(so_KIpoly_KImono, ~condition, 'full')

### WT poly vs WT mono
so_WTpoly_WTmono <- sleuth_prep(metadata_WTpoly_WTmono, extra_bootstrap_summary = TRUE, aggregation_column="gene",
                                gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))
so_WTpoly_WTmono <- sleuth_fit(so_WTpoly_WTmono, ~condition, 'full')

#WT poly vs WT total 
so_WTpoly_WTtotal <- sleuth_prep(metadata_WTpoly_WTtotal, extra_bootstrap_summary = TRUE, aggregation_column = "gene",
                                gene_mode = TRUE, extra_bootstrap_summary = TRUE, target_mapping= t2g, transform_fun_counts=function(x) log2(x + 0.5))
                                
so_WTpoly_WTtotal <- sleuth_fit(so_WTpoly_WTtotal, ~condition, 'full')                          

models(so_WTtotal_KItotal)
models(so_KIpoly_KItotal)
models(so_KIpoly_KImono)
models(so_WTpoly_WTmono)
models(so_WTpoly_WTtotal)
models(so_WTpoly_KIpoly)
models(so_WTmono_KImono)

#Wald test for WT-total vs KI-total
wald_test_wttotal_KItotal <- sleuth_wt(so_WTtotal_KItotal, which_beta="conditionHet-total")
sleuth_wald_test_WTtotal_vs_KItotal <- sleuth_results(wald_test_wttotal_KItotal, test = "conditionHet-total", show_all=TRUE)

wald_test_KIpoly_KItotal <- sleuth_wt(so_KIpoly_KItotal, which_beta="conditionHet-poly")
sleuth_wald_test_KIpoly_KItotal <- sleuth_results(wald_test_KIpoly_KItotal, test = "conditionHet-poly", show_all = TRUE)

wald_test_KIpoly_KImono <- sleuth_wt(so_KIpoly_KImono, which_beta ='conditionHet-poly')
sleuth_wald_test_KIpoly_KImono <- sleuth_results(wald_test_KIpoly_KImono, test = 'conditionHet-poly', show_all = TRUE)

wald_test_WTpoly_WTmono <- sleuth_wt(so_WTpoly_WTmono, which_beta = 'conditionWT-poly')
sleuth_wald_test_WTpoly_WTmono <- sleuth_results(wald_test_WTpoly_WTmono, test = 'conditionWT-poly', show_all = TRUE)

wald_test_WTpoly_WTtotal <- sleuth_wt(so_WTpoly_WTtotal, which_beta = 'conditionWT-poly')
sleuth_wald_test_WTpoly_WTtotal <- sleuth_results(wald_test_WTpoly_WTtotal, test = 'conditionWT-poly', show_all = TRUE)

wald_test_WTpoly_KIpoly <- sleuth_wt(so_WTpoly_KIpoly, which_beta = 'conditionHet-poly')
sleuth_wald_test_WTpoly_KIpoly <- sleuth_results(wald_test_WTpoly_KIpoly, test = 'conditionHet-poly', show_all = TRUE)

wald_test_WTmono_KImono <- sleuth_wt(so_WTmono_KImono, which_beta = 'conditionHet-mono')
sleuth_wald_test_WTmono_KImono <- sleuth_results(wald_test_WTmono_KImono, test = 'conditionHet-mono', show_all = TRUE)

#remove the NA values for genes with NA for all samples
wald_test_celf2_KI_WTtotal_vs_KItotal <- na.omit(sleuth_wald_test_WTtotal_vs_KItotal)
wald_test_celf2_KI_WTtotal_vs_KItotal$target_id <- toupper(wald_test_celf2_KI_WTtotal_vs_KItotal$target_id)

sleuth_wald_test_KIpoly_KItotal<- na.omit(sleuth_wald_test_KIpoly_KItotal)
sleuth_wald_test_KIpoly_KItotal$target_id <- toupper(sleuth_wald_test_KIpoly_KItotal$target_id)

wald_test_celf2_KIpoly_KImono <- na.omit(sleuth_wald_test_KIpoly_KImono)
wald_test_celf2_KIpoly_KImono$target_id <- toupper(wald_test_celf2_KIpoly_KImono$target_id)

wald_test_celf2_WTpoly_WTmono <- na.omit(sleuth_wald_test_WTpoly_WTmono)
wald_test_celf2_WTpoly_WTtotal <- na.omit(sleuth_wald_test_WTpoly_WTtotal)

wald_test_WTpoly_KIpoly <-na.omit(sleuth_wald_test_WTpoly_KIpoly)
wald_test_WTpoly_KIpoly$target_id<- toupper(wald_test_WTpoly_KIpoly$target_id)

sleuth_wald_test_WTmono_KImono <- na.omit(sleuth_wald_test_WTmono_KImono)
sleuth_wald_test_WTmono_KImono$target_id <- toupper(sleuth_wald_test_WTmono_KImono$target_id)

#save the Wald test for Celf 2 WTpoly vs WTmono (don't filter for protein coding genes)
write.csv(wald_test_celf2_WTpoly_WTmono, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/wald_test_celf2_WTpoly_WTmono_NOT_FILTERED_FOR_PROTEIN_CODING.csv")
write.csv(wald_test_celf2_WTpoly_WTtotal, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/wald_test_celf2_WTpoly_WTtotal_NOT_FILTERED_FOR_PROTEIN_CODING.csv")


#filter for protein coding genes only
wald_test_celf2_KI_WTtotal_vs_KItotal_protein_coding <- wald_test_celf2_KI_WTtotal_vs_KItotal[wald_test_celf2_KI_WTtotal_vs_KItotal$target_id %in% protein_coding$gene_caps,]
wald_test_celf2_KIpoly_vs_KItotal_protein_coding <- sleuth_wald_test_KIpoly_KItotal[sleuth_wald_test_KIpoly_KItotal$target_id %in% protein_coding$gene_caps,]
wald_test_celf2_KIpoly_KImono_protein_coding <-wald_test_celf2_KIpoly_KImono[wald_test_celf2_KIpoly_KImono$target_id %in% protein_coding$gene_caps,]

wald_test_WTpoly_KIpoly_protein_coding <- wald_test_WTpoly_KIpoly[wald_test_WTpoly_KIpoly$target_id %in% protein_coding$gene_caps, ]

sleuth_wald_test_WTmono_KImono_protein_coding <- sleuth_wald_test_WTmono_KImono[sleuth_wald_test_WTmono_KImono$target_id %in% protein_coding$gene_caps,]

setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/")

write.csv(wald_test_celf2_KI_WTtotal_vs_KItotal_protein_coding, "wald_test_celf2_KI_WTtotal_vs_KItotal_protein_coding.csv")
write.csv(wald_test_celf2_KIpoly_vs_KItotal_protein_coding, "wald_test_celf2_KIpoly_vs_KItotal_protein_coding.csv")
write.csv(wald_test_celf2_KIpoly_KImono_protein_coding, "wald_test_celf2_KIpoly_KImono_protein_coding.csv")
write.csv(wald_test_WTpoly_KIpoly_protein_coding, "wald_test_celf2_WTpoly_KIpoly_protein_coding.csv")

write.csv(sleuth_wald_test_WTmono_KImono_protein_coding, "sleuth_wald_test_WTmono_KImono_protein_coding.csv")

### get the normalized expression matrix
Celf2_WTtotal_KItotal_kallisto_table <- kallisto_table(so_WTtotal_KItotal, use_filtered = TRUE, normalized = TRUE)
celf2_KIpoly_KItotal_kallisto <- kallisto_table(so_KIpoly_KItotal,  use_filtered = TRUE, normalized = TRUE)
celf2_KIpoly_KImono_kallisto <- kallisto_table(so_KIpoly_KImono, use_filtered = TRUE, normalized = TRUE )
celf2_WTpoly_WTmono_kallisto <- kallisto_table(so_WTpoly_WTmono, use_filtered = TRUE, normalized = TRUE)

celf2_WT_KI_all_samples_kallisto <- kallisto_table(so_all_samples, use_filtered = TRUE, normalized = TRUE)

#keep just select columns from expression matrix
keep_celf2_wttotal_KItotal <- Celf2_WTtotal_KItotal_kallisto_table[,c("target_id", "sample", "scaled_reads_per_base")]
keep_KIpoly_KItotal <- celf2_KIpoly_KItotal_kallisto[,c("target_id", "sample", "scaled_reads_per_base")]
keep_KIpoly_KImono <- celf2_KIpoly_KImono_kallisto[,c("target_id", "sample", "scaled_reads_per_base")] 
keep_WTpoly_WTmono <- celf2_WTpoly_WTmono_kallisto[, c("target_id", "sample", "scaled_reads_per_base")]

keep_WT_KI_all_samples <- celf2_WT_KI_all_samples_kallisto[, c("target_id", "sample", "scaled_reads_per_base")]
  
###### re-format the table 
test <- keep_celf2_wttotal_KItotal  %>% spread(key = sample, value = scaled_reads_per_base)
test_KIpoly_KItotal <- keep_KIpoly_KItotal %>% spread(key = sample, value = scaled_reads_per_base)
test_KIpoly_KImono <- keep_KIpoly_KImono %>% spread(key = sample, value = scaled_reads_per_base)
test_WTpoly_WTmono <- keep_WTpoly_WTmono %>% spread(key = sample, value = scaled_reads_per_base)

test_WT_KI_all_samples <- keep_WT_KI_all_samples %>% spread(key = sample, value = scaled_reads_per_base)

#log2 transform the expression matrix
test[,2:7]<- log2(test[,2:7] + 0.5)
test_KIpoly_KItotal[,2:7] <- log2(test_KIpoly_KItotal[,2:7] + 0.5)
test_KIpoly_KImono[,2:7]<- log2(test_KIpoly_KImono[,2:7] + 0.5)
test_WTpoly_WTmono[,2:7]<- log2(test_WTpoly_WTmono[,2:7] + 0.5)

test_WT_KI_all_samples[,2:7]<-log2(test_WT_KI_all_samples[,2:7] + 0.5) 

#filter the normalized expression matrix for protein coding genes only
test_WT_KI_all_samples$target_id <- toupper(test_WT_KI_all_samples$target_id)
test_WT_KI_all_samples_protein_coding <- test_WT_KI_all_samples[test_WT_KI_all_samples$target_id %in% protein_coding$gene_caps,]

######## write the expression matrix into results directory
write.csv(test, "Celf2_WTtotal_KItotal_log2_normalized_expression_protein_coding.csv")
write.csv(test_KIpoly_KItotal, "celf2_KIpoly_KItotal_normalized_expr_protein_coding.csv")
write.csv(test_KIpoly_KImono, "celf2_KIpoly_KImono_normalized_expression_protein_codingo.csv")
write.csv(test_WTpoly_WTmono, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetectorcelf2_WTpoly_WTmono_normalized_expression.csv")

write.csv(test_WT_KI_all_samples, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/Celf2_WT_KI_all_samples_normalized_expression_no_filtering.csv" )
write.csv(test_WT_KI_all_samples_protein_coding, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/Celf2_WT_KI_all_samples_normalized_expression_WT_KI_all_samples_protein_coding.csv")

###################plot PCA plot for CELF2 WT and KI samples all

######## PCA plot for all wild type and KI samples mono, poly, total
metadata <- metadata[2:19,]

test_WT_KI_all_samples <-read.csv("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/Celf2_WT_KI_all_samples_normalized_expression_WT_KI_all_samples_protein_coding.csv")

test_WT_KI_all_samples <- test_WT_KI_all_samples[,2:20]


colnames(test_WT_KI_all_samples)[1]<- "Gene names"
test_WT_KI_all_samples$target_id <- test_WT_KI_all_samples[,1]

rownames(test_WT_KI_all_samples) = get_gene_names(test_WT_KI_all_samples)

celf2_all_samples <- select_columns(test_WT_KI_all_samples, "X.Celf2.Het.mono.2_S15_R1_001.fastq.norrna.fastq", "Celf2.WT.total.3_S7_R1_001.fastq.norrna.fastq")

#generate PCA plot with all CELF2 samples
celf2_all_samples_prcomp = prcomp(t(celf2_all_samples), scale = TRUE)$x

#for PC1 and PC2              
celf2_data_pc = data.frame(x=celf2_all_samples_prcomp[,1], y=celf2_all_samples_prcomp[,2], condition=metadata$condition, sample=row.names(celf2_all_samples_prcomp))                   

p <- (ggplot(celf2_data_pc, aes(x=x, y=y, color=condition, label = sample)) +
        geom_point(size=6) +
        scale_color_manual("Condition",
                           values = c("#A4D49C", "#4E80A5", "darkblue", "red", "gold", "purple")) +
        theme_bw() +
        #alpha = 0.1 +
        xlab("PC1") +
        ylab("PC2") +
        theme(text = element_text(size = 14)))

p + geom_text(check_overlap = TRUE, colour = "black", vjust = "inward", hjust = "inward", nudge_x = 0.05) + theme(aspect.ratio=1) 

#plot different combinations of principle components
pca_plot_combo(prcomp_df=celf2_all_samples_prcomp, first_pc = 2, second_pc = 3, first_pc_string = "2", second_pc_string = "3")
pca_plot_combo(prcomp_df=celf2_all_samples_prcomp, first_pc = 3, second_pc = 4, first_pc_string = "3", second_pc_string = "4")
pca_plot_combo(prcomp_df=celf2_all_samples_prcomp, first_pc = 4, second_pc = 5, first_pc_string = "4", second_pc_string = "5")

#################################################   Plot WT-Poly/WT-mono vs KI-Poly/KI-mono ############################

getwd()
setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/")

#read in the wald test results for Celf2-KI WT-Poly/WT-mono
wald_test_WTpoly_WTmono <- read.csv("CELF2_KI_polysome_sleuth_wald_test_WTpoly_vs_WTmono_protein_coding.csv")
#13528

#read in the table for CELF2 WT poly vs WT total 
wald_test_WTpoly_WTtotal<- read.csv("CELF2_KI_polysome_sleuth_wald_test_WTpoly_WTtotal_protein_coding.csv")

#read in the table for KI-poly vs KI-mono
wald_test_celf2_KIpoly_KImono_protein_coding<- read.csv("wald_test_celf2_KIpoly_KImono_protein_coding.csv")
#13484


#merge the Wild Type and Het (KI) Poly/Mono wald test tables 
merged_poly_vs_mono <- merge(wald_test_WTpoly_WTmono, wald_test_celf2_KIpoly_KImono_protein_coding, by.x = "target_id", by.y = "target_id", all = FALSE)

merged_WTpolyWTtotal_and_WTpolyWTmono <- merge(wald_test_WTpoly_WTtotal, wald_test_WTpoly_WTmono, by.x = "target_id", by.y = "target_id", all = FALSE )

 ggplot(merged_poly_vs_mono, aes(x=b.x, y=b.y)) +
  geom_point() +
  labs(x=("wald test b value WTpoly_WTmono"),
        y=("wald test b value KIpoly_KImono"),
       title=("CELF2-KI_Polysome-seq"))

 
 ggplot(merged_WTpolyWTtotal_and_WTpolyWTmono, aes(x=b.x, y=b.y)) +
   geom_point() +
   labs(x=("wald test b value WTpoly_vs_WTtotal"),
        y=("wald test b value WTpoly_vs_WTmono"),
        title=("CELF2-KI_Polysome-seq")) 



### overlay some of the CSDE1 binding proteins information on the CELF2 WTpoly vs WTmono  versus KIpoly_KImono plot
setwd("C:/Users/queenie.tsang/Desktop/CSDE1/ribodetector/results") 
 
wald_test_csde1_vs_input<-read.csv("wald_test_csde1_vs_input.csv")
wald_test_csde1_vs_IgG <- read.csv("wald_test_csde1_vs_IgG.csv")

#filter for genes which may be enriched in CSDE1 compared to input
# first filter for b (log2FC) > 1.33

wald_test_csde1_vs_input_b_1.33 <- wald_test_csde1_vs_input[wald_test_csde1_vs_input$b > 1.33,]
  
#filter for protein coding genes 
protein_coding <- read.csv("protein_coding_genes_biomart_ensembl.txt", header = T)
protein_coding$gene_caps <- toupper(protein_coding$Gene.name)
 
wald_test_csde1_vs_input_b_1.33$target_id <- toupper(wald_test_csde1_vs_input_b_1.33$target_id)
wald_test_csde1_vs_IgG$target_id <- toupper(wald_test_csde1_vs_IgG$target_id)

csde1_vs_input_b_1.33_protein_coding <- wald_test_csde1_vs_input_b_1.33[wald_test_csde1_vs_input_b_1.33$target_id %in% protein_coding$gene_caps, ]
csde1_vs_IgG_protein_coding <- wald_test_csde1_vs_IgG[wald_test_csde1_vs_IgG$target_id %in% protein_coding$gene_caps,]

#filter for qval < 0.05
csde1_vs_input_b_1.33_protein_coding_qval_0.05 <- csde1_vs_input_b_1.33_protein_coding[csde1_vs_input_b_1.33_protein_coding$qval < 0.05,]

write.csv(csde1_vs_input_b_1.33_protein_coding_qval_0.05, "csde1_vs_input_b_1.33_protein_coding_qval_0.05.csv")

csde1_vs_input_b_1.33_protein_coding_qval_0.05<- read.csv("csde1_vs_input_b_1.33_protein_coding_qval_0.05.csv")

###################### merge the wald test table for WTpoly-vs-WTtotal and WTpoly-vs-WTmono tables
merged_WTpolyWTtotal_and_WTpolyWTmono <- merge(wald_test_WTpoly_WTtotal, wald_test_WTpoly_WTmono, by.x = "target_id", by.y = "target_id", all = FALSE )

######### create a new column 
merged_poly_vs_mono$label <- "NA"
merged_WTpolyWTtotal_and_WTpolyWTmono$label <- "NA"

#genes upregulated in CSDE1
upregulated_csde1<- csde1_vs_input_b_1.33_protein_coding_qval_0.05$target_id
merged_poly_vs_mono<-data.frame(merged_poly_vs_mono)

merged_poly_vs_mono$label[merged_poly_vs_mono$target_id %in% upregulated_csde1]<- "up in CSDE1"

#re-plot with the genes upregulated in CSDE1 coloured
ggplot(merged_poly_vs_mono, aes(x=b.x, y=b.y, col=label, alpha= 0.3)) +
  geom_point() +
  scale_color_manual(values = c("grey30", "red" ))+
  labs(x=("wald test b value WTpoly_WTmono"),
       y=("wald test b value KIpoly_KImono"),
       title=("CELF2-KI_Polysome-seq"))

#"#A4D49C" = light green 
######### November 21 2023
# KI poly vs KI mono ratio 

setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector")
 
KIpoly_KImono <- read.csv("celf2_KIpoly_KImono_normalized_expression_protein_coding.csv")
 
#rename the column names to simplify them 
colnames(KIpoly_KImono)<-c("X", "target_id", "Celf2.Het.mono.2", "Celf2.Het.mono.1",  "Celf2.Het.mono.3",   "Celf2.Het.poly.1", "Celf2.Het.poly.2", "Celf2.Het.poly.3")


#function to calculate the poly to mono ratio for a given sample/condition
#input norm_expr_df is the normalized expression matrix which is log2 transformed from Sleuth
#poly_sample is the column name of the polysome sample
#mono_sample is the column name of the monosome sample
calculate_polymono_ratio <-function(norm_expr_df, poly_sample, mono_sample){
  norm_expr_df$new <- norm_expr_df[,poly_sample] / norm_expr_df[,mono_sample]
  norm_expr_df <<- norm_expr_df
  
  #rename the new calculations column
  #rename last column with conditions in fold change calculation
  colnames(norm_expr_df)[ncol(norm_expr_df)]<-paste(poly_sample, "_to_", mono_sample, "_ratio", sep="")
  return(norm_expr_df)
  
}

KIpoly_KImono <-calculate_polymono_ratio(norm_expr_df = KIpoly_KImono, poly_sample = "Celf2.Het.poly.1", 
                         mono_sample="Celf2.Het.mono.1")

KIpoly_KImono<-calculate_polymono_ratio(norm_expr_df = KIpoly_KImono, poly_sample = "Celf2.Het.poly.2", 
                                        mono_sample="Celf2.Het.mono.2")

KIpoly_KImono<-calculate_polymono_ratio(norm_expr_df = KIpoly_KImono, poly_sample = "Celf2.Het.poly.3", 
                                        mono_sample="Celf2.Het.mono.3")

#save the KI poly to mono calculations 
write.csv(KIpoly_KImono, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/celf2_KI_KIpoly_to_KImono_ratio_calculations.csv")

#read in the CELF2 WTpoly vs WTmono expression matrix (log2)
WTpoly_WTmono<-read.csv("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetectorcelf2_WTpoly_WTmono_normalized_expression.csv")

#simplify column names a bit
colnames(WTpoly_WTmono) <- c("X", "target_id",                               
                             "Celf2.WT.mono.1", "Celf2.WT.mono.2",
                             "Celf2.WT.mono.3", "Celf2.WT.poly.1",
                             "Celf2.WT.poly.2", "Celf2.WT.poly.3")
#calculate the poly to mono ratios for the Celf2 KI polysome wild type samples:
WTpoly_WTmono <- calculate_polymono_ratio(norm_expr_df = WTpoly_WTmono, poly_sample = "Celf2.WT.poly.1", mono_sample = "Celf2.WT.mono.1")

WTpoly_WTmono <- calculate_polymono_ratio(norm_expr_df = WTpoly_WTmono, poly_sample = "Celf2.WT.poly.2", mono_sample = "Celf2.WT.mono.2")

WTpoly_WTmono <- calculate_polymono_ratio(norm_expr_df = WTpoly_WTmono, poly_sample = "Celf2.WT.poly.3", mono_sample = "Celf2.WT.mono.3")

#save the WT poly to mono ratio calculations 
write.csv(WTpoly_WTmono, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/celf2_KI_WTpoly_to_WTmono_ratio_calculations.csv")

#merge the ratios columns for CElf2 WT poly vs mono and CELF2 KI poly vs mono together
data(mtcars)

WTpoly_WTmono_selected <- WTpoly_WTmono[,9:11]
WTpoly_WTmono_selected$target_id<-WTpoly_WTmono$target_id

KIpoly_KImono_selected <- KIpoly_KImono[,9:11]
KIpoly_KImono_selected$target_id<-KIpoly_KImono$target_id

new <- merge(WTpoly_WTmono_selected, KIpoly_KImono_selected, by.x = "target_id", by.y = "target_id")

write.csv(new, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/Celf2_KI_polysomeseq_WTpoly_WTmono_ratio_KIpoly_KImono_ratio.csv")

setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/")

celf2_poly_mono_ratios <- read.csv("Celf2_KI_polysomeseq_WTpoly_WTmono_ratio_KIpoly_KImono_ratio.csv")

##### Wald test for comparison of wild type versus KI
#https://www.biostars.org/p/316488/

row.names(celf2_poly_mono_ratios)<- celf2_poly_mono_ratios$target_id
celf2_poly_mono_ratios <- celf2_poly_mono_ratios[,3:8]

poly_mono_ratios_t <- t(celf2_poly_mono_ratios)
poly_mono_ratios_t <- as.data.frame(poly_mono_ratios_t)
condition <- data.frame(c("WTpoly_to_WTmono", "WTpoly_to_WTmono", "WTpoly_to_WTmono", "KIpoly_to_KImono", "KIpoly_to_KImono", "KIpoly_to_KImono"), poly_mono_ratios_t)
colnames(condition)[1] <- c("CaseControl")
## convert the condition column to categorical and ensure that 'WTpoly_to_WTmono' is the reference level
condition$CaseControl <- factor(condition$CaseControl, levels=c("WTpoly_to_WTmono", "KIpoly_to_KImono"))

#check the distribution of the data
hist(data.matrix(condition[,1:18354]))

#distribution of the data appears to be log-normal distribution 
#log-normal distribution is a statistical distribution of logarithmic
#values from a related normal distribution 
#?lognorm

model_condition <- glm(CaseControl ~ ZYX , data=condition, family=binomial(link="logit"))

#try a different data distribution for the glm model step
model_condition <- glm(CaseControl ~ ZYX , data=condition, family=gaussian(link="log"))

# Error in if (is.null(etastart) && is.null(start) && is.null(mustart) &&  : 
#              missing value where TRUE/FALSE needed
#              In addition: Warning message:
#                In Ops.factor(y, 0) : ‘<=’ not meaningful for factors

model_condition <- glm(CaseControl ~ ZYX , data=condition, family=gaussian(link="identity"))
# Error in glm.fit(x = c(1, 1, 1, 1, 1, 1, 1.00891305595112, 1.04836288127184,  : 
#                          NA/NaN/Inf in 'y'
#                        In addition: Warning messages:
#                          1: In Ops.factor(y, mu) : ‘-’ not meaningful for factors
#                        2: In Ops.factor(eta, offset) : ‘-’ not meaningful for factors
#                        3: In Ops.factor(y, mu) : ‘-’ not meaningful for factors
                       
#write this model part as a function so that it can be applied to all the genes in the poly:mono ratios
#dataframe 

#remove any NAs
condition <- condition[complete.cases(condition),]

library(aod)
new_pvalue_list <- list()
for (x in 2:ncol(condition)){
  gene = colnames(condition[x])
  #print(class(gene))
  #print(gene)
  #model_condition <- glm(CaseControl ~ gene, data=condition, family=binomial(link="logit"))
  model_condition <- glm(CaseControl ~ condition[,gene], data=condition, family=binomial(link="logit"))
  
  
  ratios_ZYX_wald_test <-wald.test(b=coef(model_condition), Sigma=vcov(model_condition), Terms=2)
  
  new <- data.frame(colnames(condition[,2:18355]))
  colnames(new)[1]<-"gene"
  new_t <- t(new)
  
  #save the wald test pvalue for that gene into a new column called "pvalue
  new_pvalue<-(ratios_ZYX_wald_test$result)[['chi2']][3]
  
  #append the pvalue to list
  
  new_pvalue_list <- c(new_pvalue_list, new_pvalue)
  new_pvalue_list <<- new_pvalue_list
}

#create a new dataframe with all the genes 
condition2 <- data.frame(colnames(condition[,2:ncol(condition)]))
#append the wald test pvalue to new column called "pvalue"
condition2$pvalue <- new_pvalue_list

colnames(condition2)[1]<- "gene"
colnames(condition2)[2]<- "wald_test_pvalue"

condition2$wald_test_pvalue <- as.numeric(condition2$wald_test_pvalue)

write.csv(condition2, "celf2_KI_polysomeseq_wtpoly_mono_ratios_vs_KI_poly_mono_ratios_wald_test.csv")

wald_test_celf2_wtpoly_mono_ratios_vs_KI_poly_mono_ratios<- read.csv("celf2_KI_polysomeseq_wtpoly_mono_ratios_vs_KI_poly_mono_ratios_wald_test.csv")

library(aod)


#### example from https://www.biostars.org/p/316488/
fakedata <- matrix(rbinom(10*20, 100, .1), ncol=10)
rownames(fakedata) <- paste0("sample", c(1:nrow(fakedata)))
fakedata <- data.frame(c(rep("control", 10), rep("case", 10)), fakedata)
colnames(fakedata) <- c("CaseControl", paste0("gene", c(1:10)))

head(fakedata, 14)

#check the distribution of the fake example dataset 
hist(data.matrix(fakedata[,2:ncol(fakedata)]))

fakedata$CaseControl <- factor(fakedata$CaseControl, levels=c("control","case"))
model <- glm(CaseControl ~ gene1, data=fakedata, family=binomial(link="logit"))

library(aod)
install.packages("aod")

coef(model_condition)
# (Intercept)      ZYX 
# 5.486773   -5.309659

#We then apply the Wald test on whichever model coefficient ('term')
#we want. Our gene is the second coefficient:

library(aod)
ratios_ZYX_wald_test <-wald.test(b=coef(model_condition), Sigma=vcov(model_condition), Terms=2)
# Wald test:
#   ----------
#   
#   Chi-squared test:
#   X2 = 0.0073, df = 1, P(> X2) = 0.93

(ratios_ZYX_wald_test$result)[['chi2']]
#       chi2          df           P 
# 0.007348734 1.000000000 0.931685220 

(ratios_ZYX_wald_test$result)[['chi2']][3]
#         P 
# 0.9316852 

################################################ try using Limma for the statistical test for the ratios

setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/")
celf2_poly_mono_ratios <- read.csv("Celf2_KI_polysomeseq_WTpoly_WTmono_ratio_KIpoly_KImono_ratio.csv")

library(limma)
library(dplyr)

row.names(celf2_poly_mono_ratios) <- celf2_poly_mono_ratios$target_id

#keep just the sample columns 
celf2_poly_mono_ratios <- select_columns(celf2_poly_mono_ratios, "Celf2.WT.poly.1_to_Celf2.WT.mono.1_ratio", "Celf2.Het.poly.3_to_Celf2.Het.mono.3_ratio")
celf2_poly_mono_ratios <- as.matrix(celf2_poly_mono_ratios)

#create a metadata table for the Celf2-WT-poly-to-mono ratios and the KI-poly-to-mono ratios
metadata_ratios = data.frame(sample = colnames(celf2_poly_mono_ratios))
metadata_ratios$condition <- c("WT", "WT", "WT", "Het", "Het", "Het")

write.csv(metadata_ratios, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/metadata_celf2_KI_polysome_seq_poly_mono_ratios_metadata.csv")
write.csv(celf2_poly_mono_ratios, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/celf2_poly_mono_ratios_matrix.csv")

celf2_poly_mono_ratios <- read.csv("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/celf2_poly_mono_ratios_matrix.csv")
metadata_ratios <- read.csv("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/metadata_celf2_KI_polysome_seq_poly_mono_ratios_metadata.csv")


##filter for just protein coding genes
colnames(celf2_poly_mono_ratios)[1]<- "gene"


celf2_poly_mono_ratios$gene <- toupper(celf2_poly_mono_ratios$gene)
celf2_poly_mono_ratios<- celf2_poly_mono_ratios[celf2_poly_mono_ratios$gene %in% protein_coding$gene_caps,] 
#13312 genes remaining after filtering for protein coding genes only\


######### differential expression analysis with Limma
### Differential expression analysis between cell types
difference = ifelse(metadata_ratios$condition == "WT", 1, 0)
design = cbind(control=1, difference=difference)
WT_KI_ratios_top_table = get_top_table(celf2_poly_mono_ratios, design)

write.table(WT_KI_ratios_top_table, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/Celf2_KI_polysomeseq_WTpolymono_ratios_vs_KIpolymono_ratios_top_table.txt", sep = "\t", quote = F, row.names = T, col.names = T)

write.table(WT_KI_ratios_top_table, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/Celf2_KI_polysomeseq_WTpolymono_ratios_vs_KIpolymono_ratios_top_table_protein_coding.txt", sep = "\t", quote = F, row.names = T, col.names = T)

print(volcano_plot(WT_KI_ratios_top_table, "Celf2_KI_polysomeseq_WTpolymono_ratios_vs_KIpolymono_ratios"))


#### read in limma top table results for stats test for WTpoly-mono ratios vs KIpoly-mono ratios
setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/")
WT_KI_ratios_top_table <- read.csv("Celf2_KI_polysomeseq_WTpolymono_ratios_vs_KIpolymono_ratios_top_table.txt", sep = "\t", header = T)


##### try plotting the average ratios of CELF2 WTpoly/mono vs KIpoly/mono with the CSDE1 upregulated genes overlayed
celf2_poly_mono_ratios$WT_poly_to_mono_avg_ratio <- rowMeans(celf2_poly_mono_ratios[,c(3:5)])
celf2_poly_mono_ratios$KI_poly_to_mono_avg_ratio <- rowMeans(celf2_poly_mono_ratios[,c(6:8)])

### add a labels column to identify the genes upregulated in CSDE1:
celf2_poly_mono_ratios$csde1_up <- "NA"

celf2_poly_mono_ratios$csde1_up[celf2_poly_mono_ratios$target_id %in% upregulated_csde1]<- "up in CSDE1"

write.csv(celf2_poly_mono_ratios, "celf2_poly_mono_ratios_csde1_labels.csv")


#filter the celf2_poly_mono_ratios table for protein coding genes only
celf2_poly_mono_ratios_protein_coding <- celf2_poly_mono_ratios[celf2_poly_mono_ratios$target_id %in% protein_coding$gene_caps,]


ggplot(celf2_poly_mono_ratios_protein_coding, aes(x=WT_poly_to_mono_avg_ratio, y=KI_poly_to_mono_avg_ratio, col=csde1_up, alpha= 0.5)) +
  geom_point() +
  scale_color_manual(values = c("grey30", "red" ))+
  labs(x=("average WTpoly_WTmono ratio"),
       y=("average KIpoly_KImono ratio"),
       title=("CELF2-KI_Polysome-seq"))

#overlay transcript length on this plot
setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq")

#read in the table for ENSEMBL protein coding genes and which have the transcript length info
protein_coding <- read.csv("mouse_ensembl_biomart_mouse_genes_GRCm39.txt")

protein_coding$gene_caps <- toupper(protein_coding$Gene.name)

merge_celf2_genes <- merge(celf2_poly_mono_ratios_protein_coding, protein_coding, by.x="target_id", by.y="gene_caps")

#### so even for the coding protein genes for CELF2 dataframe, each protein coding gene may have multiple transcripts
# need to select the longest transcript for each gene
# there are 12351 unique protein coding genes in the CELF2 dataframe, 
#but 72361 rows when merged with transcripts info
#indicating multiple transcripts per gene



### new dataframe for storing the selected transcripts for each gene
new <- data.frame(colnames(merge_celf2_genes))
new <-rbind(new, selected_keep)
new_t <- t(new)
colnames(new_t)<- colnames(merge_celf2_genes)
new_t <- data.frame(new_t)


unique_protein_coding_genes <- (unique(merge_celf2_genes$target_id))

for (gene in unique_protein_coding_genes){
  
  #select rows containing specific gene 
  selected <- merge_celf2_genes[merge_celf2_genes$target_id == gene,]
  
  # Sort the selected genes dataframe by transcripts length 
  selected_ordered <- selected[order(selected$Transcript.length..including.UTRs.and.CDS., decreasing = TRUE),]
  
  #keep just the longest transcript for the gene
  selected_keep <- selected_ordered[1,]
  
  #add this row for the longest transcript into new dataframe
  new_t <- rbind(new_t, selected_keep)
  new_t <<- new_t
  
}

new_filtered <- na.omit(new_t)
new2 <- new_filtered[2:12375,]

write.csv(new2, "results/celf2_wtpolymono_ratios_kipolymono_ratios_longest_transcript.csv")

new2$Transcript.length..including.UTRs.and.CDS. <-as.numeric(new2$Transcript.length..including.UTRs.and.CDS.)

#log2 transform the transcript length
new2$log2_transcript_length <- log2(new2$Transcript.length..including.UTRs.and.CDS.)

new2$WT_poly_to_mono_avg_ratio <- as.numeric(new2$WT_poly_to_mono_avg_ratio)

new2$KI_poly_to_mono_avg_ratio<- as.numeric(new2$KI_poly_to_mono_avg_ratio)


#plot this transcript length with the CELF2-KI-Polysome seq data
ggplot(new2, aes(x=WT_poly_to_mono_avg_ratio, y=KI_poly_to_mono_avg_ratio, col=log2_transcript_length, alpha= 0.5)) +
  geom_point() +
  labs(x=("average WTpoly_WTmono ratio"),
       y=("average KIpoly_KImono ratio"),
       title=("CELF2-KI_Polysome-seq")) +
  scale_colour_gradient(low = "blue", high = "red", na.value = NA)

####################overlay CSDE vs IgG binding b value on the CELF2 WTpoly vs WT total (x-axis) and WTpoly vs WTmono (y-axis) plot 

wald_test_csde1_vs_IgG$target_id <- toupper(wald_test_csde1_vs_IgG$target_id) 

#filter for protein coding genes 
csde1_vs_IgG_protein_coding <- wald_test_csde1_vs_IgG[wald_test_csde1_vs_IgG$target_id %in% protein_coding$gene_caps, ]

###### merge CSDE vs IgG b value to to merged CELF2 WTpoly vs WT total and WTpoly vs WTmono dataframe
merged_celf2_wt_and_csde1_vs_IgG <- merge(merged_WTpolyWTtotal_and_WTpolyWTmono, csde1_vs_IgG_protein_coding, by.x = "target_id", by.y = "target_id", all = FALSE)

legend_title <- "b value CSDE1 vs IgG"

p <- ggplot(merged_celf2_wt_and_csde1_vs_IgG, aes(x=b.x, y=b.y, col=b, alpha=0.1)) +
  geom_point() +
  labs(x=("wald test b value WT poly vs WT total"),
       y=("wald test b value WT poly vs WT mono"),
      
       title=("CELF2-KI_Polysome-seq")) +
  scale_colour_gradient(legend_title, low = "grey", high = "darkblue", na.value = NA)

plot(p)

################## for this plot try putting the b value CSDE1 vs IgG into categories

# b value > 2 : enriched in CSDE1
# b value < -2: enriched in IgG
# b value between -2 and 2

merged_celf2_wt_and_csde1_vs_IgG$enrichment <- "NA"

library(dplyr)
library(stringr)

merged_celf2_wt_and_csde1_vs_IgG_new <-merged_celf2_wt_and_csde1_vs_IgG %>% 
  mutate(enrichment = case_when(
    b > 2 ~ "enriched in CSDE1",
    b < -2 ~ "enriched in IgG",
    (b > -2 & b < 2) ~ "not enriched"))

#re-plot scatterplot with new CSDE1 enrichment categories, filtered for protein coding genes
p <- ggplot(merged_celf2_wt_and_csde1_vs_IgG_new %>%
              arrange(desc(enrichment)),                 ## this alters the order in which the points are plotted, so that CSDE1 plots are on top
            aes(x=b.x, y=b.y, col=enrichment, alpha=0.1)) +
  scale_color_manual(values = c( "red", "darkblue","grey"))+
  geom_point(size = 3) +
  labs(x=("wald test b value WT poly vs WT total"),
       y=("wald test b value WT poly vs WT mono"),
       
       title=("CELF2-KI_Polysome-seq - filtered for protein coding genes")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  #scale_colour_gradient(legend_title, low = "grey", high = "darkblue", na.value = NA)

plot(p)

####### filtering out genes which are not enriched and keeping only genes which are enriched in either CSDE1 or IgG

merged_celf2_wt_and_csde1_vs_IgG_no_NA <- merged_celf2_wt_and_csde1_vs_IgG_new[merged_celf2_wt_and_csde1_vs_IgG_new$enrichment == "enriched in CSDE1" | merged_celf2_wt_and_csde1_vs_IgG_new$enrichment == "enriched in IgG",]
merged_celf2_wt_and_csde1_vs_IgG_no_NA$b.x <- as.numeric(merged_celf2_wt_and_csde1_vs_IgG_no_NA$b.x)
merged_celf2_wt_and_csde1_vs_IgG_no_NA$b.y <- as.numeric(merged_celf2_wt_and_csde1_vs_IgG_no_NA$b.y)

p <- ggplot(merged_celf2_wt_and_csde1_vs_IgG_no_NA, aes(x=b.x, y=b.y, col=enrichment, alpha=0.1)) +
  geom_point() +
  labs(x=("wald test b value WT poly vs WT total"),
       y=("wald test b value WT poly vs WT mono"),
       title=("CELF2-KI_Polysome-seq - filtered for protein coding genes")) +
       geom_point(size = 5)

#scale_colour_gradient(legend_title, low = "grey", high = "darkblue", na.value = NA)

plot(p)

#filter for CSDE1 enrichment with b value > 4 , and for IgG enrichment with b value , -4
merged_celf2_wt_and_csde1_vs_IgG_cutoff_bvalue_4 <-merged_celf2_wt_and_csde1_vs_IgG %>% 
  mutate(enrichment = case_when(
    b > 4 ~ "enriched in CSDE1",
    b < -4 ~ "enriched in IgG",
    (b > -4 & b < 4) ~ "not enriched"))

p <- ggplot(merged_celf2_wt_and_csde1_vs_IgG_cutoff_bvalue_4, aes(x=b.x, y=b.y, col=enrichment, alpha=0.1)) +
  geom_point() +
  labs(x=("wald test b value WT poly vs WT total"),
       y=("wald test b value WT poly vs WT mono"),
       
       title=("CELF2-KI_Polysome-seq - filtered for protein coding genes")) #+
#scale_colour_gradient(legend_title, low = "grey", high = "darkblue", na.value = NA)

plot(p)

### show only enriched genes at b value > 4 or b value < -4
merged_celf2_wt_and_csde1_vs_IgG_cutoff_bvalue_4_enriched_only <- merged_celf2_wt_and_csde1_vs_IgG_cutoff_bvalue_4[merged_celf2_wt_and_csde1_vs_IgG_cutoff_bvalue_4$enrichment %in% c("enriched in CSDE1", "enriched in IgG"),]

#filter for CSDE1 enrichment with b value > 4 , and for IgG enrichment with b value , -4
p <- ggplot(merged_celf2_wt_and_csde1_vs_IgG_cutoff_bvalue_4_enriched_only, aes(x=b.x, y=b.y, col=enrichment, alpha=0.1)) +
  geom_point() +
  labs(x=("wald test b value WT poly vs WT total"),
       y=("wald test b value WT poly vs WT mono"),
       
       title=("CELF2-KI_Polysome-seq - filtered for protein coding genes")) #+
plot(p)

merged_celf2_wt_and_csde1_vs_IgG_no_NA<- na.omit(merged_celf2_wt_and_csde1_vs_IgG_no_NA)
merged_celf2_wt_and_csde1_vs_IgG_no_NA <- as.matrix(merged_celf2_wt_and_csde1_vs_IgG_no_NA)

#sort by enrichment status 
df2 <- emp_df[order(df$price),]

merged_celf2_wt_and_csde1_vs_IgG_no_NA <- as.data.frame(merged_celf2_wt_and_csde1_vs_IgG_no_NA)
df2 <- merged_celf2_wt_and_csde1_vs_IgG_no_NA[order(merged_celf2_wt_and_csde1_vs_IgG_no_NA$enrichment),]

df2 <- as.matrix(df2)

##### save the enriched table for CSDE1 vs IgG cutoff value of b <-2 and b > 2
write.csv(df2, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/merged_celf2_wtpoly_WTtotal_vs_WTpoly_vs_WTmono_and_csde1_vs_IgG_cutoff_bvalue_2.csv")



############## plot CELF2 WT poly vs WT total (x-axis) vs CELF2 WT poly vs WT mono (y-axis) with overlay of wald test b value of CSDE1 vs IgG 
#### DO NOT filter for protein coding genes 

#read in the CELF2 WTpoly vs WT total wald test table values:
setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/")

celf2_WTpoly_WTmono <- read.csv("wald_test_celf2_WTpoly_WTmono_NOT_FILTERED_FOR_PROTEIN_CODING.csv")
celf2_WTpoly_WTtotal <- read.csv("wald_test_celf2_WTpoly_WTtotal_NOT_FILTERED_FOR_PROTEIN_CODING.csv")

setwd("C:/Users/queenie.tsang/Desktop/CSDE1/ribodetector/results") 
wald_test_csde1_vs_IgG <- read.csv("wald_test_csde1_vs_IgG.csv")

#filter for genes which have greater than b value of 2 in wald test CSDE1 vs IgG
csde1_up_bvalue_2 <- wald_test_csde1_vs_IgG[wald_test_csde1_vs_IgG$b > 2, ]
IgG_up_bvalue_2 <- wald_test_csde1_vs_IgG[wald_test_csde1_vs_IgG$b < -2,]

#first merge celf2_WTpoly_WTmono and celf2_WTpoly_WTtotal dataframes
celf2_WTpoly_vs_WTtotal_versus_WTpoly_vs_WTmono_not_filterd_by_protein_coding <- merge(celf2_WTpoly_WTtotal, celf2_WTpoly_WTmono, by.x= "target_id", by.y = "target_id", all = FALSE )

#merge CSDE1_vs_IGG with the combined WTpoly_vs_WTtotal_merged_with_WTpoly_vs_WTmono
WTpoly_vs_WTtotal_versus_WTpoly_vs_WTmono_csde1_IgG_no_protein_coding_filtering <-merge(celf2_WTpoly_vs_WTtotal_versus_WTpoly_vs_WTmono_not_filterd_by_protein_coding, wald_test_csde1_vs_IgG, by.x= "target_id", by.y = "target_id", all = FALSE)

#filter for protein coding genes 
WTpoly_vs_WTtotal_versus_WTpoly_vs_WTmono_csde1_IgG_no_protein_coding_filtering$target_id <- toupper(WTpoly_vs_WTtotal_versus_WTpoly_vs_WTmono_csde1_IgG_no_protein_coding_filtering$target_id)
merged_celf2_WTpoly_vs_WTtotal_versus_WTpoly_vs_WTmono_csde1_IgG_protein_coding<- WTpoly_vs_WTtotal_versus_WTpoly_vs_WTmono_csde1_IgG_no_protein_coding_filtering[WTpoly_vs_WTtotal_versus_WTpoly_vs_WTmono_csde1_IgG_no_protein_coding_filtering$target_id %in% protein_coding$gene_caps,]

#save the merged CELF2 WTpoly-WTtotal vs WTpoly-WTmono with CSDE1 vs IgG b values filtered for protein coding genes 
write.csv(merged_celf2_WTpoly_vs_WTtotal_versus_WTpoly_vs_WTmono_csde1_IgG_protein_coding, "C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/merged_celf2_WTpoly_vs_WTtotal_versus_WTpoly_vs_WTmono_csde1_IgG_protein_coding.csv")

legend_title <- "b value CSDE1 vs IgG"

p <- ggplot(WTpoly_vs_WTtotal_versus_WTpoly_vs_WTmono_csde1_IgG_no_protein_coding_filtering, aes(x=b.x, y=b.y, col=b, alpha=0.1)) +
  geom_point() +
  labs(x=("wald test b value WT poly vs WT total"),
       y=("wald test b value WT poly vs WT mono"),
       
       title=("CELF2-KI_Polysome-seq - no filtering for protein coding genes")) +
  scale_colour_gradient(legend_title, low = "grey", high = "darkblue", na.value = NA)

plot(p)

##### try doing overlay of the CSDE1 genes which are enriched (b value > 2)


merged_celf2_wt_and_csde1_vs_IgG_NOT_filtered_for_protein_coding <-WTpoly_vs_WTtotal_versus_WTpoly_vs_WTmono_csde1_IgG_no_protein_coding_filtering %>% 
  mutate(enrichment = case_when(
    b > 2 ~ "enriched in CSDE1",
    b < 2 ~ "not enriched"))

#adjust the b value cutoff value considered for enrichment in CSDE1 to 1.5
merged_celf2_csde1_protein_coding <-merged_celf2_WTpoly_vs_WTtotal_versus_WTpoly_vs_WTmono_csde1_IgG_protein_coding
merged_celf2_csde1_protein_coding <- merged_celf2_csde1_protein_coding %>%
  mutate(enrichment = case_when(
    b > 1.5 ~ "enriched in CSDE1",
    b < 1.5 ~ "not enriched"))


#save this table to use later for changes with plot aesthetics etc.
setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/")

write.csv(merged_celf2_wt_and_csde1_vs_IgG_NOT_filtered_for_protein_coding, "merged_celf2_wt_and_csde1_vs_IgG_NOT_filtered_for_protein_coding.csv")


celf2_csde1_scatterplot(merged_table = merged_celf2_wt_and_csde1_vs_IgG_NOT_filtered_for_protein_coding)
celf2_csde1_scatterplot3(merged_table = merged_celf2_wt_and_csde1_vs_IgG_NOT_filtered_for_protein_coding)

celf2_csde1_scatterplot(merged_table = merged_celf2_wt_and_csde1_vs_IgG_NOT_filtered_for_protein_coding) + celf2_csde1_scatterplot2(merged_table = merged_celf2_wt_and_csde1_vs_IgG_NOT_filtered_for_protein_coding)

#not filtered for protein coding genes
celf2_csde1_scatterplot3(merged_table = merged_celf2_wt_and_csde1_vs_IgG_NOT_filtered_for_protein_coding) + 
  geom_point(data = subset(merged_celf2_wt_and_csde1_vs_IgG_NOT_filtered_for_protein_coding, b > 2),
             stat = 'identity',
       #aes(x=b.x, y=b.y, alpha=0.1, fill=b, size = 5)) +
       aes(x=b.x, y=b.y, alpha=0.1, colour=b, size = 5)) +
  #scale_fill_gradient(low="cyan4", high = "coral1")+
  scale_colour_gradient(low= "cyan4", high="coral", legend_title = )

  #scale_fill_gradient2(low="cyan4", mid = "white", high = "coral1")+
  #geom_point() +
  labs(x=("wald test b value WT poly vs WT total"),
       y=("wald test b value WT poly vs WT mono"),
       
       title=("CELF2-KI_Polysome-seq - filtered for protein coding genes")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

####################### do the same thing but now DO FILTER FOR PROTEIN CODING GENES

merged_celf2_wt_and_csde1_vs_IgG_NOT_filtered_for_protein_coding$target_id <-toupper(merged_celf2_wt_and_csde1_vs_IgG_NOT_filtered_for_protein_coding$target_id)
merged_celf2_wt_and_csde1_vs_IgG_protein_coding <-merged_celf2_wt_and_csde1_vs_IgG_NOT_filtered_for_protein_coding[merged_celf2_wt_and_csde1_vs_IgG_NOT_filtered_for_protein_coding$target_id %in% protein_coding$gene_caps,]
  
#re-plot this with the table filtered for PROTEIN CODING GENES:
legend_title<-("b value CSDE1 vs IgG")

#sorting the merged wald test b value table by b value will allow the genes with higher B value (more upregulated in CSDE1)
#to be plotted last when rendering so that these higher CSDE1 targets appear on top of the other genes (easier to visualize the upregulated genes)
merged_celf2_wt_and_csde1_vs_IgG_protein_coding_sort_b <- merged_celf2_wt_and_csde1_vs_IgG_protein_coding[order(merged_celf2_wt_and_csde1_vs_IgG_protein_coding$b),]
merged_celf2_csde1_protein_coding_sort_b <- merged_celf2_csde1_protein_coding[order(merged_celf2_csde1_protein_coding$b),]
  
  
write.csv(merged_celf2_wt_and_csde1_vs_IgG_protein_coding_sort_b, "merged_celf2_wt_and_csde1_vs_IgG_protein_coding_sort_b.csv")

celf2_csde1_scatterplot3(merged_table = merged_celf2_wt_and_csde1_vs_IgG_protein_coding_sort_b)+
  geom_point(data = subset(merged_celf2_wt_and_csde1_vs_IgG_protein_coding_sort_b, b > 2 & qval <= 0.01),
             stat = 'identity',
             shape=16,
             aes(x=b.x, y=b.y, alpha=0.1, colour=b, size = 5)) + #shape16 should not have outline around the points 
  scale_colour_gradient(low= "cyan4", high="coral", legend_title)+
  labs(x=("wald test b value WT poly vs WT total"),
     y=("wald test b value WT poly vs WT mono"),
     title=("CELF2-KI_Polysome-seq - filtered for protein coding genes; b value cutoff > 2, qvalue <= 0.01")) +
    
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
##### try filtering by different b value cutoffs 

#b value cutoff of 1.5 
celf2_csde1_scatterplot3(merged_table = merged_celf2_csde1_protein_coding_sort_b)+
  geom_point(data = subset(merged_celf2_csde1_protein_coding_sort_b, b > 1.5 & qval <= 0.01),
             stat = 'identity',
             shape=16,
             aes(x=b.x, y=b.y, alpha=0.1, colour=b, size = 5)) + #shape16 should not have outline around the points 
  scale_colour_gradient(low= "cyan4", high="coral", legend_title)+
  labs(x=("wald test b value WT poly vs WT total"),
       y=("wald test b value WT poly vs WT mono"),
       title=("CELF2-KI_Polysome-seq - filtered for protein coding genes; b value cutoff > 1.5, qvalue <= 0.01")) +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(alpha= none, size= none)  ### remove the transparency and point size info from the legend 


#b value cutoff of 0.58 
celf2_csde1_scatterplot3(merged_table = merged_celf2_wt_and_csde1_vs_IgG_protein_coding_sort_b)+
  geom_point(data = subset(merged_celf2_wt_and_csde1_vs_IgG_protein_coding_sort_b, b > 0.58 & qval <= 0.01),
             stat = 'identity',
             shape=16,
             aes(x=b.x, y=b.y, alpha=0.1, colour=b, size = 5)) + #shape16 should not have outline around the points 
  scale_colour_gradient(low= "cyan4", high="coral", legend_title) +
  labs(x=("wald test b value WT poly vs WT total"),
     y=("wald test b value WT poly vs WT mono"),
     title=("CELF2-KI_Polysome-seq - protein coding genes; b value cutoff > 0.58, qvalue <= 0.01")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#b value cutoff of 1 
celf2_csde1_scatterplot3(merged_table = merged_celf2_wt_and_csde1_vs_IgG_protein_coding_sort_b)+
  geom_point(data = subset(merged_celf2_wt_and_csde1_vs_IgG_protein_coding_sort_b, b > 1 & qval <= 0.01),
             stat = 'identity',
             shape=16,
             aes(x=b.x, y=b.y, alpha=0.1, colour=b, size = 5)) + #shape16 should not have outline around the points 
  scale_colour_gradient(low= "cyan4", high="coral", legend_title) +
  labs(x=("wald test b value WT poly vs WT total"),
       y=("wald test b value WT poly vs WT mono"),
       title=("CELF2-KI_Polysome-seq - protein coding genes; b value cutoff > 1, qvalue <= 0.01")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

###### try different q value cutoffs: 0.05 bvalue cutoff of 1:



################## December 13 2023
# Double checking WT vs KI calculations with volcano plots

### read in the wald test results for WT total vs KI total 
setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/results/ribodetector/")

wald_test_celf2_KI_WTtotal_vs_KItotal <- read.csv("wald_test_celf2_KI_WTtotal_vs_KItotal_protein_coding.csv")

#plot this on volcano plot 
print(volcano_plot_sleuth_wald_test(top_table=wald_test_celf2_KI_WTtotal_vs_KItotal, title="CELF2_WT_total_vs_KI_total_protein_coding_", fold_change_cutoff=1.33, pvalue_cutoff=0.05))

##### wald test for WT polysome vs KI polysome 

metadata_WTpoly_KIpoly <- metadata[metadata$condition %in% c("WT-poly", "Het-poly"),]

#volcano plot for WT polysome vs KI polysome
print(volcano_plot_sleuth_wald_test(top_table = wald_test_WTpoly_KIpoly_protein_coding, title = "CELF2_WTpoly_vs_KIpoly_protein_coding", fold_change_cutoff = 1.33, pvalue_cutoff = 0.05))

#volcano plot for WT monosome vs KI monosome
print(volcano_plot_sleuth_wald_test(top_table= sleuth_wald_test_WTmono_KImono_protein_coding, title = "CELF2_WTmono_vs_KImono_protein_coding", fold_change_cutoff = 1.33, pvalue_cutoff = 0.05))






##### Example of how to subset only some of the data for colouring, and colouring using a gradient based on column values:
####### example https://stackoverflow.com/questions/47398391/apply-ggplot-color-scale-gradient-to-portion-of-data
tr_sim <- data.frame(site_id = seq(1,100,1), estimated_impact = 
                       rnorm(100,18,200), impact_group = rep(c(1,2),each = 50))

rng_full <- range(tr_sim$estimated_impact)
#produce plot showing the full range of impacts across all sites and then 
#over the subsetted sites

impact_plot_full <- ggplot(data = tr_sim, aes(x = factor(site_id, levels = 
                                                           site_id[order(estimated_impact)]), y = estimated_impact)) +
  geom_bar(stat = "identity",width = 1, fill = "grey80") 

impact_plot_full 

impact_plot_full + 
  geom_bar(stat = "identity", width = 1, position = "stack", aes(y  = 
                                                                   estimated_impact[impact_group == 1])) +
  scale_fill_gradient2(low="firebrick", mid="yellow", high = "green4") +
  labs(y = "Estimated Impact ($/week)", x = "Total number of sites with estimate 
  is 100", title = "Sites with the greatest impact after coverage loss") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(breaks = 
                       round(seq(rng_full[1],rng_full[2],by=100),digits=0)) 

impact_plot_full + 
  geom_bar(data = subset(tr_sim, estimated_impact < 0), 
           stat = "identity",
           aes(y = estimated_impact, fill = estimated_impact)) + 
  scale_fill_gradient2(low = "firebrick", mid = "yellow", high = 
                         "green4") +
  theme(axis.text.x = element_blank()) +
  xlab("site_id")
