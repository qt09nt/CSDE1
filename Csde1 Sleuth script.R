# Derived from Dr.Gordon Paul R-script History
# Modified by Drayden Kopp for Dr. Guang Yang's RIP-Seq Data
# Date of acquisition of Paul's Script: November 19, 2019
# Modified by Reza Aghanoori June 2022

#-------------------------------- Data Fetching -----------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rhdf5", force = TRUE)
install.packages("devtools")
remotes::install_github("pachterlab/sleuth#260")
library("sleuth")
#install.packages("ggplot2")
library(ggplot2)

# folder containing all *.kallisto folders used with a metadata table
#setwd("/Users/rezaaghanoori/Desktop/June2022, Csde1 RIP Seq data/CSDE1 Kallisto") 
setwd("C:/Users/queenie.tsang/Desktop/June2022, Csde1 RIP Seq data/1.CSDE1 Kallisto")

sample_id <- dir(pattern = ".kallisto")
kal_dirs <- sample_id
kal_dirs

#import data from metadata table
s2c <- read.csv(file.path("sampleinfo.csv"), header = TRUE)
#assigning covariates 
s2c <- dplyr::select(s2c, sample = Sample, Condition = State, Csde = Csde, Capture = Capture)
#link the metadata with kallisto paths
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

#Plot PCA 
sleuth_live(so)

#wald test for Csde versus IgG/Input
so <- sleuth_wt(so, "Csde")
wald_Csde <- sleuth_results(so, "Csde")


#Wald test for Capture vs Csde/Input
so <- sleuth_wt(so, "Capture")
wald_Capture <- sleuth_results(so, "Capture")


wald_Capture

#Genes which with FDR<.05 and ln(FC) > 0 for Capture vs Csde/input & Csde vs Capture/Input
##Sept 21 2023 EDIT: Genes which with FDR<.05 and log 2 (FC) > 0 for Capture vs Csde/input & Csde vs Capture/Input
Capture_ids <-  wald_Capture$target_id[which(wald_Capture$qval<.05 & wald_Capture$b > 0)]
Csde_ids <- wald_Csde$target_id[which(wald_Csde$qval < 0.05 & wald_Csde$b > 0)]

#match the names of the rows for both conditions
rownames(wald_Capture) <- wald_Capture$target_id
rownames(wald_Csde) <- wald_Csde$target_id

#filter genes which are expressed higher in Capture than Csde
Csde_ids_b_net_pos <- Csde_ids[which(-1*wald_Capture[Csde_ids,]$b < wald_Csde[Csde_ids,]$b)]  #subset the CSDE1 genes where the log2FC of Capture genes is less than the CSDE1 log 2 FC 
Csde_ids_b_net_pos1 <- Csde_ids[which(wald_Capture[Csde_ids,]$b < wald_Csde[Csde_ids,]$b)]  #subset the CSDE1 genes where the the log2FC of Capture (input) genes is less than the CSDE1 log2FC

Csde_ids_b_net_pos <- intersect(Csde_ids_b_net_pos, Csde_ids_b_net_pos1)

length(Csde_ids_b_net_pos)
#5243

#Comparison of Capture without Csde to IgG
so <- sleuth_fit(so, ~Capture, "no_Csde")
so <- sleuth_lrt(so, "no_Csde", "full")
so <- sleuth_lrt(so, "no_Csde", "full")
lrt_Capture <- sleuth_results(so, "no_Csde:full", test_type="lrt")
table(lrt_Capture$qval < .05)
lrt_Csde <- sleuth_results(so, "no_Csde:full", test_type="lrt")
lrt_Csde_ids <- lrt_Csde$target_id[which(lrt_Csde$qval < .05)]
length(intersect(Csde_ids_b_net_pos, lrt_Csde_ids))
#5243
write.table(wald_Capture, file = "wald.txt", sep="\t", quote=FALSE)

wald_Csde

#Generate Table of results,
write.table(data.frame(wald_Csde[Csde_ids_b_net_pos,], wald_Capture[Csde_ids_b_net_pos, c("qval","b")]), sep="\t", quote=FALSE, "Correct_Csde_Input_pos.txt")

#table for wald capture
write.table(wald_Capture, "Corrected_Wald_capture.txt", sep = "\t", quote = FALSE)
#log2 transformed table for wald capture 
write.table(wald_Capture, "Corrected_Wald_capture_log2.txt", sep = "\t", quote = FALSE)


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


########## look into software to remove rRNA genes  
#https://github.com/hzi-bifo/RiboDetector
# https://www.biostars.org/p/159959/
