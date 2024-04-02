# Derived from Dr.Gordon Paul R-script History
# Modified by Drayden Kopp for Dr. Guang Yang's RIP-Seq Data
# Date of acquisition of Paul's Script: November 19, 2019
# Modified by Reza Aghanoori June 2022
#modified by Queenie Tsang September 2023

library(tidyr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)

#-------------------------------- Data Fetching -----------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rhdf5", force = TRUE)
install.packages("devtools")
remotes::install_github("pachterlab/sleuth#260")
library(sleuth)
#install.packages("ggplot2")
library(ggplot2)


#read in functions 
setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/R")
source("Cefl2_polysome_seq_functions.R")


# folder containing all *.kallisto folders used with a metadata table
#setwd("/Users/rezaaghanoori/Desktop/June2022, Csde1 RIP Seq data/CSDE1 Kallisto") 
setwd("C:/Users/queenie.tsang/Desktop/June2022, Csde1 RIP Seq data/1.CSDE1 Kallisto")

setwd("C:/Users/queenie.tsang/Desktop/CSDE1/ribodetector")


#load in the kallisto files with fastq sample files filtered to remove rRNA  using Ribodetector

sample_id <- dir(pattern = ".kallisto")
kal_dirs <- sample_id
kal_dirs

#import data from metadata table
s2c <- read.csv(file.path("sampleinfo.csv"), header = TRUE)

#ribodetector version of metadata table
s2c <- read.csv(file.path("sampleinfo_ribodetector.csv"), header=TRUE)

#assigning covariates 
s2c <- dplyr::select(s2c, sample = Sample, Condition = State, Csde = Csde, Capture = Capture)


#link the metadata with kallisto paths : 
s2c <- dplyr::mutate(s2c, path = kal_dirs)

#confirm the metadata file matches the correct file location
print(s2c)

# -------------------------------Defining the Annotation-------------------------------------

t2g <- read.table("refseq2gene_mouse", colClasses=rep("character", 2), header=TRUE)
t2g <- dplyr::distinct(t2g, target_id, .keep_all= TRUE)

#-------------------------------------filtering function---------------------------------------

#Possibly a redundant filtering method. Refer to ?sleuth_prep
#filter out any features that do not have at least 5 estimated counts in at least 47
design_filter <- function(design, row, min_reads=5, min_prop = 0.47){
  sum(apply(design, 2, function(x){
    y <- as.factor(x);
    return(max(tapply(row, y, function(f){sum(f >= min_reads)})/
                 tapply(row, y, length)) == 1 
           || basic_filter(row, min_reads, min_prop)
    )
  })) > 0}


#---------Generate count table for Csde and Capture-----------------------------------------

#desiginated the matrix to be used for sleuth_prep, verify the conditions
design_matrix <- model.matrix(~Capture+Csde, s2c)
design_matrix

# creating the prep covariates

library(sleuth)

so <- sleuth_prep(s2c, ~Capture+Csde, aggregation_column="gene", gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, filter_fun=function(x){design_filter(design_matrix, x)})

#sept 21 2023 edit:  transform counts using log2 instead of the default natural log, using the parameter: transform_fun_counts=function (x) log2(x+0.5) 
so <- sleuth_prep(s2c, ~Capture+Csde, aggregation_column="gene", gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5), filter_fun=function(x){design_filter(design_matrix, x)})

so <- sleuth_fit(so)
#########Oct 13 2023 analyze CSDE1 vs IgG and CSDE1 vs Capture separately
#subset s2c for CSDE1 and IgG
meta_csde_igg <-s2c[s2c$Condition %in% c("IgG", "Csde"),]
design_matrix_csde_igg <- design_matrix[1:6,]

#subset metadata table for just CSDE1 and Input samples
meta_csde_capture <- s2c[s2c$condition %in% c("Csde", "Capture"),]



##### first relevel the factors so that IgG is the base level 
meta_csde_igg$condition <- as.factor(meta_csde_igg$condition)
meta_csde_igg$condition <- relevel(meta_csde_igg$condition, ref= "IgG")

so_csde_igg <- sleuth_prep(meta_csde_igg, aggregation_column = "gene", gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts= function (x) log2(x+0.5), filter_fun=function(x){design_filter(design_matrix_csde_igg, x)})
so_csde_igg <- sleuth_fit(so_csde_igg, ~Condition, 'full')

#### for CSDE vs Input samples 
meta_csde_capture$Condition <- as.factor(meta_csde_capture$Condition)
meta_csde_capture$Condition <- relevel(meta_csde_capture$Condition, ref = "Capture")


#create the sleuth object
so_csde_capture <- sleuth_prep(meta_csde_capture, aggregation_column = "gene", gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping = t2g, transform_fun_counts= function(x) log2(x+0.5), filter_fun=function(x){design_filter(design_matrix_csde_input, x)})
so_csde_capture <- sleuth_fit(so_csde_capture, ~Condition, 'full')

#Plot PCA 
sleuth_live(so)
sleuth_live(so_csde_igg)

#look at the genes contributing the most to the variances in PC1 and PC2
plot_loadings(so, pc_input=1)
plot_loadings(so, pc_input=2)

#wald test for Csde1 versus IgG differential expression analysis 
csde_vs_igg_wt <-sleuth_wt(so_csde_igg, which_beta = "ConditionCsde")

#output results for csde1 vs IgG wald test
sleuth_wald_test_csde_vs_IgG <- sleuth_results(csde_vs_igg_wt, test = "ConditionCsde", show_all = TRUE)

#remove the genes with NA for all samples 
new_sleuth_wald_test_csde_vs_IgG <- na.omit(sleuth_wald_test_csde_vs_IgG)
write.csv(new_sleuth_wald_test_csde_vs_IgG, "results/wald_test_csde1_vs_IgG.csv")

#wald test for Csde1 versus Capture differential expression analysis
csde1_vs_input_wt <- sleuth_wt(so_csde_capture, which_beta = "ConditionCsde")
sleuth_wald_test_csde_vs_input <- sleuth_results(csde1_vs_input_wt, test = "ConditionCsde", show_all = TRUE)

#remove all the genes which have NA for all the samples
new_sleuth_wald_test_csde_vs_input <- na.omit(sleuth_wald_test_csde_vs_input)
write.csv(new_sleuth_wald_test_csde_vs_input, "results/wald_test_csde1_vs_input.csv")

##### try filtering with different threshold values for qvalue and b (log2 fold change)

#### filter with q value =< 0.01 and b value (log2 fold change) > 1.3
wald_test_csde_vs_IgG_qvalue_0.01 <- new_sleuth_wald_test_csde_vs_IgG[new_sleuth_wald_test_csde_vs_IgG$qval <= 0.01, ]
wald_test_csde_vs_IgG_bvalue_1.3 <- wald_test_csde_vs_IgG_qvalue_0.01[wald_test_csde_vs_IgG_qvalue_0.01$b >1.3, ]
#1775 genes

wald_test_csde_vs_Input_qvalue_0.01 <- new_sleuth_wald_test_csde_vs_input[new_sleuth_wald_test_csde_vs_input$qval <= 0.01, ]
wald_test_csde_vs_Input_bvalue_1.3 <- wald_test_csde_vs_Input_qvalue_0.01[wald_test_csde_vs_Input_qvalue_0.01$b > 1.3,]
#2355 genes

#### see which genes overlap between the CSDE vs IgG and CSDE vs Input genes when filtered by qvalue <= 0.01 and b value (log 2 fold change) > 1.3
csde1_enriched_qvalue_0.01_bvalue_1.3 <- intersect(wald_test_csde_vs_IgG_bvalue_1.3$target_id, wald_test_csde_vs_Input_bvalue_1.3$target_id) 
write.csv(csde1_enriched_qvalue_0.01_bvalue_1.3, "results/csde1_overlap_qvalue_0.01_bvalue_1.3.csv")
#1034 genes 


############ try filtering by qvalue <= 0.01 and b value (log 2 )
csde_vs_IgG_bvalue_2 <- wald_test_csde_vs_IgG_qvalue_0.01[wald_test_csde_vs_IgG_qvalue_0.01$b > 2,]
csde_vs_input_bvalue_2 <- wald_test_csde_vs_Input_qvalue_0.01[wald_test_csde_vs_Input_qvalue_0.01$b > 2,]

#see which genes overlap between CSDE vs IgG and CSDE1 vs Input genes when filtered by qvalue <=0.01 and 
csde_enriched_qvalue_0.01_bvalue_2 <- intersect(csde_vs_IgG_bvalue_2$target_id, csde_vs_input_bvalue_2$target_id)

#write.csv
write.csv(csde_enriched_qvalue_0.01_bvalue_2, "results/csde1_overlap_qvalue_0.01_bvalue_2.csv")




############ grouped analysis ie. Csde vs IgG/Input may not be suitable for the purposes of this current paper


#wald test for Csde versus IgG/Input
so <- sleuth_wt(so, "Csde")
wald_Csde <- sleuth_results(so, "Csde")


#Wald test for Capture vs Csde/Input
so <- sleuth_wt(so, "Capture")
wald_Capture <- sleuth_results(so, "Capture")


wald_Capture

#Genes which with FDR<.05 and ln(FC) > 0 for Capture vs Csde/input & Csde vs Capture/Input
##Sept 21 2023 EDIT: Genes which with FDR<.05 and log 2 (FC) > 0 for Capture vs Csde/input & Csde vs Capture/Input
#adjust log 2 fold change cutoff to 1.3 (which is 2.5 X fold change)
Capture_ids <-  wald_Capture$target_id[which(wald_Capture$qval<.05 & wald_Capture$b > 1.3)]
Csde_ids <- wald_Csde$target_id[which(wald_Csde$qval < 0.05 & wald_Csde$b > 1.3)]

##### try different values for log 2 fold change cutoff values ie. fold change cut off of 1.5 ( which is equal to 0.58 for log 2 fold change)
Capture_ids_log2fc_0.58 <-  wald_Capture$target_id[which(wald_Capture$qval<.05 & wald_Capture$b > 0.58)]
Csde_ids_log2fc_0.58 <- wald_Csde$target_id[which(wald_Csde$qval < 0.05 & wald_Csde$b > 0.58)]

#match the names of the rows for both conditions
rownames(wald_Capture) <- wald_Capture$target_id
rownames(wald_Csde) <- wald_Csde$target_id

#filter genes which are expressed higher in Capture than Csde
Csde_ids_b_net_pos <- Csde_ids[which(-1*wald_Capture[Csde_ids,]$b < wald_Csde[Csde_ids,]$b)]  #subset the CSDE1 genes where the log2FC of Capture genes is less than the CSDE1 log 2 FC 
Csde_ids_b_net_pos1 <- Csde_ids[which(wald_Capture[Csde_ids,]$b < wald_Csde[Csde_ids,]$b)]    #subset the CSDE1 genes where the log2FC of Capture (input) genes is less than the CSDE1 log2FC
Csde_ids_b_net_pos <- intersect(Csde_ids_b_net_pos, Csde_ids_b_net_pos1)

Csde_ids_b_net_pos_log2fc_0.58 <- Csde_ids[which(-1*wald_Capture[Csde_ids_log2fc_0.58,]$b < wald_Csde[Csde_ids_log2fc_0.58,]$b)]  #subset the CSDE1 genes where the log2FC of Capture genes is less than the CSDE1 log 2 FC 
Csde_ids_b_net_pos1_log2fc_0.58 <- Csde_ids[which(wald_Capture[Csde_ids_log2fc_0.58,]$b < wald_Csde[Csde_ids_log2fc_0.58,]$b)]    #subset the CSDE1 genes where the log2FC of Capture (input) genes is less than the CSDE1 log2FC
Csde_ids_b_net_pos_log2fc_0.58 <- intersect(Csde_ids_b_net_pos_log2fc_0.58, Csde_ids_b_net_pos1_log2fc_0.58)


length(Csde_ids_b_net_pos)
#5243
# for log2 fold change cutoff of 1.3, this value is 2180

#for log 2 fold change of 0.58: 
length(Csde_ids_b_net_pos_log2fc_0.58)
#2618

#Comparison of Capture without Csde to IgG
so <- sleuth_fit(so, ~Capture, "no_Csde")
so <- sleuth_lrt(so, "no_Csde", "full")

lrt_Capture <- sleuth_results(so, "no_Csde:full", test_type="lrt")
table(lrt_Capture$qval < .05)
lrt_Csde <- sleuth_results(so, "no_Csde:full", test_type="lrt")
lrt_Csde_ids <- lrt_Csde$target_id[which(lrt_Csde$qval < .05)]

length(intersect(Csde_ids_b_net_pos, lrt_Csde_ids))
#5243
#2180

length(intersect(Csde_ids_b_net_pos_log2fc_0.58, lrt_Csde_ids))
#2617

write.table(wald_Capture, file = "results/wald.txt", sep="\t", quote=FALSE)

write.table(wald_Capture, file = "results/wald_test_Capture.txt", sep="\t", quote=FALSE)

wald_Csde

#Generate Table of results,
write.table(data.frame(wald_Csde[Csde_ids_b_net_pos,], wald_Capture[Csde_ids_b_net_pos, c("qval","b")]), sep="\t", quote=FALSE, "results/Correct_Csde_Input_pos_log2fc_cutoff_1.3.txt")
write.table(data.frame(wald_Csde[Csde_ids_b_net_pos_log2fc_0.58,], wald_Capture[Csde_ids_b_net_pos_log2fc_0.58, c("qval","b")]), sep="\t", quote=FALSE, "results/Correct_Csde_Input_pos_log2fc_cutoff_0.58.txt")


#table for wald capture
write.table(wald_Capture, "results/Corrected_Wald_capture.txt", sep = "\t", quote = FALSE)
#log2 transformed table for wald capture 
write.table(wald_Capture, "results/Corrected_Wald_capture_log2.txt", sep = "\t", quote = FALSE)


#table for wald Csde
write.table(wald_Csde, "Corrected_Wald_Csde.txt", sep = "\t", quote = FALSE)
#table for wald CSDE log2
write.table(wald_Csde, "Corrected_Wald_Csde_log2.txt", sep = "\t", quote = FALSE)

sleuth_live(so)

#table of all counts

write.table(sleuth_to_matrix(so, "obs_raw", "scaled_reads_per_base"), "Corrected_all_counts.txt", sep = "\t", quote = FALSE)

write.table(sleuth_to_matrix(so, "obs_raw", "tpm"), "Corrected_all_tpm.txt", sep = "\t", quote = FALSE)

length(sig_ids)
sleuth_live(so)


############# Plot scatter plot of Input versus RIP (CSDE1) Sept 21 2023


so_matrix <- sleuth_to_matrix(obj= so, which_df = "obs_norm", which_units = "scaled_reads_per_base")
so_matrix <- as.data.frame(so_matrix)

write.csv(so_matrix, "CSDE1_norm_counts_scaled_reads_per_base.csv")

so_matrix <- read.csv("CSDE1_norm_counts_scaled_reads_per_base.csv")
so_matrix <- as.data.frame(so_matrix)
row.names(so_matrix) <- so_matrix[,1]

##########for Input vs RIP (CSDE1)

#subset for CSDE1 and Input columns
library(dplyr)
csde_values  <- so_matrix[, c("Csde1",  "Csde2",  "Csde3")]
input_values <- so_matrix [,c("Input1", "Input2", "Input3")]
IGG_values <- so_matrix[,c("IgG1", "IgG2", "IgG3")]

#find the average value of the scaled_reads_per_base for the 3 replicates of CSDE1 samples for each of the genes
csde_values$average_csde <- rowMeans(csde_values)
input_values$average_input <- rowMeans(input_values)
IGG_values$average_IGG <- rowMeans(IGG_values)

#merge the average scaled reads per base values for csde samples and input samples into a single df
merged_csde_input <- data.frame(csde_values$average_csde, input_values$average_input)

#merge the average scaled reads per base values for CSDE samples and IGG samples into a single df
merged_csde_IGG <- data.frame(csde_values$average_csde, IGG_values$average_IGG)

row.names(merged_csde_input) <-row.names(csde_values)
row.names(merged_csde_IGG)<- row.names(csde_values)


#create a scatter plot for Input average scaled_reads_per_base versus CSDE1 scaled_reads_per_base 
library(ggplot2)

View(merged_csde_input)
write.csv(merged_csde_input, "merged_csde_input_average_scaled_reads_per_base.csv")
write.csv(merged_csde_IGG, "merged_csde_IGG_average_scaled_reads_per_base.csv")

merged_csde_IGG <- read.csv("merged_csde_IGG_average_scaled_reads_per_base.csv")

ggplot(merged_csde_input, aes(x=input_values.average_input, y=csde_values.average_csde)) + 
  geom_point()


ggplot(merged_csde_IGG, aes(x=IGG_values.average_input, y=csde_values.average_csde)) + 
  geom_point()

dev.off()
dev.new()

####### remove outliers that make the plot hard to visualize the majority of the genes

input_sorted<- input_values[order(input_values$average_input, decreasing=TRUE),]
input_outliers<- row.names(head(input_sorted)) 
#"Rn7s1"  "Rn7s2"  "RN7SK"  "Rpph1"  "Rn45s"  "MALAT1"
input_values_filtered <- input_values[!(row.names(input_values) %in% input_outliers),] 

#see the what the top expressed genes are in the IGG samples 
IGG_sorted <- IGG_values[order(IGG_values$average_IGG, decreasing = TRUE),]
IGG_outliers<- row.names(head(IGG_sorted))
#"Rn45s"  "Rn7s1"  "Rn7s2"  "Rn28s1" "RN7SK"  "SAMD5"


csde_sorted <- csde_values[order(csde_values$average_csde, decreasing=TRUE),]
head(csde_sorted)
#"Rn45s"        "Rn28s1"       "NCAN"         "RMRP"         "LRP1"         "LOC101056014"
csde_outliers <-c("Rn45s","Rn28s1")

outliers <- c(csde_outliers, input_outliers)

#merge the average scaled reads per base values for csde samples and input samples into a single df
merged_csde_input <- data.frame(csde_values$average_csde, input_values$average_input)
row.names(merged_csde_input) <-row.names(csde_values)
csde_input_filtered <- merged_csde_input[!(row.names(merged_csde_input) %in% outliers),]

#remove the outliers noncoding RNA from the merged df of the average scaled reads per base values for CSDE and IGG samples 
merged_csde_IGG_filtered <- merged_csde_IGG[!(row.names(merged_csde_IGG) %in% IGG_outliers), ]


#log2 transform the scaled_reads_per_base 
csde_IGG_filtered_log2 <- log2(merged_csde_IGG_filtered+0.5)

ggplot(csde_IGG_filtered_log2, aes(x=IGG_values.average_IGG, y=csde_values.average_csde)) + 
  geom_point()


ggplot(csde_input_filtered, aes(x=input_values.average_input, y=csde_values.average_csde)) + 
  geom_point()


#try log 2 transforming the scaled_reads_per_base_values and replot scatterplot
csde_input_filtered_log2 <-log2(csde_input_filtered+0.5)

ggplot(csde_input_filtered_log2, aes(x=input_values.average_input, y=csde_values.average_csde)) + 
  geom_point()



######### October 6 2023 PCA plot with ribodetector filtered samples
#
so_matrix <- sleuth_to_matrix(obj= so, which_df = "obs_norm", which_units = "scaled_reads_per_base")
so_matrix <- as.data.frame(so_matrix)

#remove the genes which have 0 variance between the samples
so_matrix <- so_matrix[ which(apply(so_matrix, 1, var) != 0),]
samples_log2 <- log2(so_matrix+0.5)

str(csde1_pca)

#PCA
csde1_pca <-prcomp(t(samples_log2), center = TRUE, scale = TRUE)$x
summary(csde1_pca)

#plot PC1 by PC2 
#Oct 16 2023 - make modifications to PCA plot based on suggested changes from meeting
pc_data = data.frame(x= csde1_pca[,1], y = csde1_pca[,2], sample = row.names(csde1_pca))

pc_data$condition <- c("Csde1", "Csde1", "Csde1", "IgG", "IgG", "IgG", "Input", "Input", "Input")
  
p <- ggplot (pc_data, aes(x=x, y=y, color = condition, label = sample))+
  geom_point(size=15, alpha = 0.7) +    #transparency adjusted with "alpha" parameter
  scale_color_manual("Sample", values = c("blue", "darkgreen", "red") )+
  theme_bw() +
  xlab("PC1") +
  ylab("PC2") + 
  theme(text = element_text(size = 16), 
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        aspect.ratio=1)  #make the plot square

p + geom_text()


############# Perform Gene Set Enrichment Analysis 

#first extract the matrix of normalized values from Sleuth
so_matrix <- sleuth_to_matrix(obj= so, which_df = "obs_norm", which_units = "scaled_reads_per_base")
so_matrix <- as.data.frame(so_matrix)

#remove the genes which have 0 variance between the samples
so_matrix <- so_matrix[ which(apply(so_matrix, 1, var) != 0),]
samples_log2 <- log2(so_matrix+0.5)

#write the log 2 matrix of all samples into a csv file
write.csv(samples_log2, "results/csde1_RIP_all_samples_ribodetector_log2.csv")

#filter for just the genes which are high in CSDE1 mRNAs
csde1_genes <- samples_log2[Csde_ids_b_net_pos,]

row.names(csde1_genes)

########### volcano plot for CSDE1 vs IgG  December 20 2023

setwd("C:/Users/queenie.tsang/Desktop/CSDE1/ribodetector/results")

csde1_vs_IgG_wald_test <- read.csv("wald_test_csde1_vs_IgG.csv")

#read in protein coding genes list from ENSEMBL
protein_coding <- read.csv("protein_coding_genes_biomart_ensembl.txt")
protein_coding$gene_caps <- toupper(protein_coding$Gene.name)

# filter the CSDE1 vs IgG for protein coding genes 
csde1_vs_IgG_wald_test$target_id <- toupper(csde1_vs_IgG_wald_test$target_id)

csde1_vs_IgG_protein_coding <- csde1_vs_IgG_wald_test[csde1_vs_IgG_wald_test$target_id %in% protein_coding$gene_caps,]

volcano_plot_sleuth_wald_test(top_table = csde1_vs_IgG_protein_coding, title = "CSDE1 vs IgG wald test qvalue cutoff 0.0001 log2fc cutoff ", fold_change_cutoff= 1.5, pvalue_cutoff= 0.001)

volcano_plot2(top_table = csde1_vs_IgG_protein_coding, title= "CSDE1 vs IgG ")

volcano_plot2 (csde1_vs_IgG_protein_coding, title= "CSDE1 vs IgG pvalue_cutoff = 0.0001, log2fc_cutoff = 1.5", pvalue_cutoff = 0.0001, log2fc_cutoff = 1.5)

library(ggrepel)
library(ggplot2)

#for the plot just label select genes of interest, and italicize the labels
lab_italics <- paste0("italic('", csde1_vs_IgG_protein_coding$target_id, "')") 

selectLab_italics = paste0(
  "italic('",
  c('Tle4', 'Tbr1', 'Pou3f3', 'Satb2', 'Fezf2', 'Fezf1', 'Bcl11b'), "')")

#for the top table, convert gene names to Capital first letter and lowercase the rest of the gene name to indicate it's a mouse gene (not human)
csde1_vs_IgG_protein_coding$target_id[csde1_vs_IgG_protein_coding$target_id == 'BCL11B' ] <- 'Bcl11b'
csde1_vs_IgG_protein_coding$target_id[csde1_vs_IgG_protein_coding$target_id == 'TLE4' ] <- 'Tle4'
csde1_vs_IgG_protein_coding$target_id[csde1_vs_IgG_protein_coding$target_id == 'TBR1' ] <- 'Tbr1'
csde1_vs_IgG_protein_coding$target_id[csde1_vs_IgG_protein_coding$target_id == 'POU3F3' ] <- 'Pou3f3'
csde1_vs_IgG_protein_coding$target_id[csde1_vs_IgG_protein_coding$target_id == 'SATB2' ] <- 'Satb2'
csde1_vs_IgG_protein_coding$target_id[csde1_vs_IgG_protein_coding$target_id == 'FEZF2' ] <- 'Fezf2'
csde1_vs_IgG_protein_coding$target_id[csde1_vs_IgG_protein_coding$target_id == 'FEZF1' ] <- 'Fezf1'

#selectLabel = c('Rpl7', 'Rpl34', 'Rpl13a', 'Gapdh', 'Actb', 'Camk2a', 'Rbfox1', 'Bcl11b', 'Tubb3', 'Tuba1a')
selectLabel =  c('Tle4', 'Tbr1', 'Pou3f3', 'Satb2', 'Fezf2', 'Fezf1', 'Bcl11b')

#volcano_plot_sleuth(top_table=WT_mono_vs_poly_wald_test, title="WT-mono_versus_WT-poly_Wald_Test", fold_change_cutoff=1.32, pvalue_cutoff=0.05)
volcano_plot_sleuth(top_table=csde1_vs_IgG_protein_coding, title="CSDE1 vs IgG protein coding", fold_change_cutoff=1.5, pvalue_cutoff=0.0001)

