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
krausher_WTpoly_WTtotal_meta <- krausher_meta[krausher_meta$condition %in% c("P_P0_WT", "T_P0_WT"),]
krausher_WTmono_WTtotal_meta <- krausher_meta[krausher_meta$condition %in% c("M_P0_WT", "T_P0_WT"),]

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

##### WT poly vs WT total 
krausher_WTpoly_WTtotal_meta$condition <- as.factor(krausher_WTpoly_WTtotal_meta$condition)
krausher_WTpoly_WTtotal_meta$condition <- relevel(krausher_WTpoly_WTtotal_meta$condition, ref = "T_P0_WT")

###### WT mono vs WT total
krausher_WTmono_WTtotal_meta$condition <- as.factor(krausher_WTmono_WTtotal_meta$condition)
krausher_WTmono_WTtotal_meta$condition <- relevel(krausher_WTmono_WTtotal_meta$condition, ref= "T_P0_WT")

#double check the condition factor levels
krausher_mono_poly_meta$condition
krausher_WTmono_WTtotal_meta$condition

#####wild type poly vs mono - create sleuth object
so_WTpoly_vs_WTmono <- sleuth_prep(krausher_mono_poly_meta, extra_bootstrap_summary = TRUE, aggregation_column="gene", gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))
so_WTpoly_vs_WTtotal <- sleuth_prep(krausher_WTpoly_WTtotal_meta, extra_bootstrap_summary = TRUE, aggregation_column="gene", gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))
so_WTmono_vs_WTtotal <- sleuth_prep(krausher_WTmono_WTtotal_meta, extra_bootstrap_summary = TRUE, aggregation_column="gene", gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))

#fit the full model
so_WTpoly_vs_WTmono <- sleuth_fit(so_WTpoly_vs_WTmono, ~condition, 'full')
so_WTpoly_vs_WTtotal <- sleuth_fit(so_WTpoly_vs_WTtotal, ~condition, 'full')
so_WTmono_vs_WTtotal <- sleuth_fit(so_WTmono_vs_WTtotal, ~condition, 'full')

models(so_WTpoly_vs_WTmono)
models(so_WTpoly_vs_WTtotal)
models(so_WTmono_vs_WTtotal)

#Wald Test for differential expression of WT mono versus WT total 
wald_test_WTpoly_vs_WTmono <- sleuth_wt(so_WTpoly_vs_WTmono, which_beta = 'conditionP_P0_WT')
wald_test_WTpoly_vs_WTtotal <- sleuth_wt(so_WTpoly_vs_WTtotal, which_beta= 'conditionP_P0_WT')
wald_test_WTmono_vs_WTtotal <- sleuth_wt(so_WTmono_vs_WTtotal, which_beta= 'conditionM_P0_WT')

#output wald test results
sleuth_wald_test_WTpoly_vs_WTmono <- sleuth_results(wald_test_WTpoly_vs_WTmono, test = 'conditionP_P0_WT', show_all=TRUE)
sleuth_wald_test_WTpoly_vs_WTtotal <- sleuth_results(wald_test_WTpoly_vs_WTtotal, test = 'conditionP_P0_WT', show_all=TRUE)
sleuth_wald_test_WTmono_vs_WTtotal <- sleuth_results(wald_test_WTmono_vs_WTtotal, test = 'conditionM_P0_WT', show_all=TRUE)

#remove genes which are just NA for all samples
sleuth_wald_test_WTpoly_vs_WTmono_filtered <- na.omit(sleuth_wald_test_WTpoly_vs_WTmono)
sleuth_wald_test_WTpoly_vs_WTtotal_filtered <- na.omit(sleuth_wald_test_WTpoly_vs_WTtotal)
sleuth_wald_test_WTmono_vs_WTtotal_filtered <- na.omit(sleuth_wald_test_WTmono_vs_WTtotal)

#filter for protein coding genes only
sleuth_wald_test_WTpoly_vs_WTmono_filtered$target_id <- toupper(sleuth_wald_test_WTpoly_vs_WTmono_filtered$target_id)
sleuth_wald_test_WTpoly_vs_WTmono_protein_coding <- sleuth_wald_test_WTpoly_vs_WTmono_filtered[sleuth_wald_test_WTpoly_vs_WTmono_filtered$target_id %in% protein_coding$gene_caps,]

sleuth_wald_test_WTpoly_vs_WTtotal_filtered$target_id <- toupper(sleuth_wald_test_WTpoly_vs_WTtotal_filtered$target_id)
sleuth_wald_test_WTpoly_vs_WTtotal_protein_coding <- sleuth_wald_test_WTpoly_vs_WTtotal_filtered[sleuth_wald_test_WTpoly_vs_WTtotal_filtered$target_id %in% protein_coding$gene_caps,]

sleuth_wald_test_WTmono_vs_WTtotal_filtered$target_id <- toupper(sleuth_wald_test_WTmono_vs_WTtotal_filtered$target_id)
sleuth_wald_test_WTmono_vs_WTtotal_protein_coding <- sleuth_wald_test_WTmono_vs_WTtotal_filtered[sleuth_wald_test_WTmono_vs_WTtotal_filtered$target_id %in% protein_coding$gene_caps,]


setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Kraushar et al 2014 dataset/results/")

write.csv(sleuth_wald_test_WTpoly_vs_WTmono_protein_coding, "C:/Users/queenie.tsang/Desktop/CSDE1/Kraushar et al 2014 dataset/results/sleuth_wald_test_WTpoly_vs_WTmono_protein_coding.csv")
write.csv(sleuth_wald_test_WTpoly_vs_WTtotal_protein_coding, "sleuth_wald_test_WTpoly_vs_WTtotal_protein_coding.csv")
write.csv(sleuth_wald_test_WTmono_vs_WTtotal_protein_coding, "sleuth_wald_test_WTmono_vs_WTtotal_protein_coding.csv")

########################## 

#################################### get the normalized Kallisto Expression Matrix

krausher_meta_WTmono_WTpoly <- krausher_meta[krausher_meta$condition %in% c("M_P0_WT", "P_P0_WT"),]
krausher_meta_WTmono_WTpoly$condition <- as.factor(krausher_meta_WTmono_WTpoly$condition)
krausher_meta_WTmono_WTpoly$condition <- relevel(krausher_meta_WTmono_WTpoly$condition, ref="M_P0_WT")

krausher_meta_WTmono_WTpoly$condition

so_WTpoly_vs_WTmono <- sleuth_prep(krausher_meta_WTmono_WTpoly, extra_bootstrap_summary = TRUE, aggregation_column="gene", gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))

#get the kallisto table 
WTpoly_vs_WTmono_normalized_expr <- kallisto_table(so_WTpoly_vs_WTmono, use_filtered = FALSE, normalized = TRUE)
WTpoly_vs_WTtotal_normalized_expr <- kallisto_table(so_WTpoly_vs_WTtotal, use_filtered = FALSE, normalized = TRUE)
WTmono_vs_WTtotal_norm_expr <- kallisto_table(so_WTmono_vs_WTtotal, use_filtered = FALSE, normalized = TRUE)

#re-format this table to have samples as columns
WTpoly_vs_WTmono_normalized_expr_select <- WTpoly_vs_WTmono_normalized_expr[,c("sample", "target_id", "scaled_reads_per_base")]
WTpoly_vs_WTmono_expr_reformatted <-spread(WTpoly_vs_WTmono_normalized_expr_select, sample, scaled_reads_per_base)

WTpoly_vs_WTtotal_normalized_expr_selected <- WTpoly_vs_WTtotal_normalized_expr[,c("sample", "target_id", "scaled_reads_per_base")]
WTpoly_vs_WTtotal_reformatted <- spread(WTpoly_vs_WTtotal_normalized_expr_selected, sample, scaled_reads_per_base)

WTmono_vs_WTtotal_norm_expr <-WTmono_vs_WTtotal_norm_expr[,c("sample", "target_id", "scaled_reads_per_base")]
WTmono_vs_WTtotal_norm_expr_reformatted <- spread(WTmono_vs_WTtotal_norm_expr,  sample, scaled_reads_per_base)

#filter the Kallisto table expression matrix for protein coding genes only
WTpoly_vs_WTmono_expr_reformatted$target_id <- toupper(WTpoly_vs_WTmono_expr_reformatted$target_id)
WTpoly_vs_WTmono_expr_protein_coding <- WTpoly_vs_WTmono_expr_reformatted[WTpoly_vs_WTmono_expr_reformatted$target_id %in% protein_coding$gene_caps,]

WTpoly_vs_WTtotal_reformatted$target_id<- toupper(WTpoly_vs_WTtotal_reformatted$target_id)
WTpoly_vs_WTtotal_normalized_expr_protein_coding <- WTpoly_vs_WTtotal_reformatted[WTpoly_vs_WTtotal_reformatted$target_id %in% protein_coding$gene_caps, ] 

WTmono_vs_WTtotal_norm_expr_reformatted$target_id <- toupper(WTmono_vs_WTtotal_norm_expr_reformatted$target_id)
WTmono_vs_WTtotal_expr_protein_coding <- WTmono_vs_WTtotal_norm_expr_reformatted[WTmono_vs_WTtotal_norm_expr_reformatted$target_id %in% protein_coding$gene_caps,]

#remove rows where expression is 0 across all samples
ans = WTpoly_vs_WTmono_expr_protein_coding[rowSums(WTpoly_vs_WTmono_expr_protein_coding[,2:7])>0,] 
ans2 <- ans

ans = WTpoly_vs_WTtotal_normalized_expr_protein_coding[rowSums(WTpoly_vs_WTtotal_normalized_expr_protein_coding[,2:7])>0,]
ans3 <- ans

ans = WTmono_vs_WTtotal_expr_protein_coding[rowSums(WTmono_vs_WTtotal_expr_protein_coding[,2:7])>0, ]
ans4 <- ans

#log2 transform the expression matrix 
ans2[,2:7] <- log2(ans2[,2:7] + 0.5)
WTpoly_vs_WTmono_expr_protein_coding_log2 <- ans2

#log2 transform expression matrix WTpoly vs WT total
ans3[,2:7]<- log2(ans3[,2:7] + 0.5)
WTpoly_vs_WTtotal_expr_matrix_protein_coding_log2 <- ans3

ans4[,2:7]<-log2(ans4[,2:7]+0.5)
WTmono_vs_WTtotal_expr_matrix_protein_coding_log2 <-ans4


#write the expression matrix to the results directory:
setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Kraushar et al 2014 dataset/results/")

write.csv(WTpoly_vs_WTmono_expr_protein_coding_log2, "Krausher_WTpoly_vs_WTmono_expr_protein_coding_log2.csv")
write.csv(WTpoly_vs_WTtotal_expr_matrix_protein_coding_log2,"WTpoly_vs_WTtotal_expr_matrix_protein_coding_log2.csv")
write.csv(WTmono_vs_WTtotal_expr_matrix_protein_coding_log2, "WTmono_vs_WTtotal_expr_matrix_protein_coding_log2.csv")



##################### try to overlay CDS information over the mono/total (x-axis) and poly/total (y-axis) for Krausher et al  (2014) dataset 

##### read in the Krausher et al (2014) wild type mono/total table and poly/total tables from Kallisto-Sleuth
setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Kraushar et al 2014 dataset/results/")

krausher_wtmono_Wttotal <- read.csv("sleuth_wald_test_WTmono_vs_WTtotal_protein_coding.csv", header=TRUE)
krausher_wtpoly_wttotal <- read.csv("sleuth_wald_test_WTpoly_vs_WTtotal_protein_coding.csv", header= TRUE)

#read in ENSEMBL Biomart Mouse genes (GRCm39) protein coding genes information
genes <- read.csv("mouse_ensembl_biomart_mouse_genes_GRCm39.txt", header=TRUE)

#first plot the mono/total numbers on the x-axis and the poly/total numbers on the y-axis from 
#the Krausher (2014) dataset
#merge into single dataframe before plotting

merged_krausher <- merge(krausher_wtmono_Wttotal, krausher_wtpoly_wttotal, by.x = "target_id", by.y = "target_id")

library(ggplot2)
library(RColorBrewer)

ggplot(merged_krausher, aes(x=b.x, y=b.y)) +
  geom_point() +
  theme_bw() +
  ylab("Wald test b value Krausher et al (2014) WT-poly vs WT-total") +
  xlab("Wald test b value Krausher et al (2014) WT-mono vs WT-total")

#convert gene names to capitals to match those of the merged_krausher table
genes$Gene.name <- toupper(genes$Gene.name)

merge_krausher_genes <- merge(merged_krausher, genes, by.x="target_id", by.y="Gene.name")

ggplot(merge_krausher_genes, aes(x=b.x, y=b.y)) +
  geom_point(aes(colour = Transcript.length..including.UTRs.and.CDS.,  alpha = 0.5 )) +  #use alpha parameter to adjust transparency of the dots
  theme_bw() +
  #scale_fill_gradient(low="black", high="green")+
  scale_colour_gradient2(low = "#3288bd", mid = "#fee08b", high = "#d53e4f", midpoint = median(merge_krausher_genes$`Transcript.length..including.UTRs.and.CDS.`, rm.na = TRUE))+
  ylab("Wald test b value Krausher et al (2014) WT-poly vs WT-total") +
  xlab("Wald test b value Krausher et al (2014) WT-mono vs WT-total")

set.seed(1)
df <- data.frame(
  x = runif(100),
  y = runif(100),
  z1 = rnorm(100),
  z2 = abs(rnorm(100))
)

df_na <- data.frame(
  value = seq(1, 20),
  x = runif(20),
  y = runif(20),
  z1 = c(rep(NA, 10), rnorm(10))
)

# Default colour scale colours from light blue to dark blue
ggplot(df, aes(x, y)) +
  geom_point(aes(colour = z2))


# For diverging colour scales use gradient2
ggplot(df, aes(x, y)) +
  geom_point(aes(colour = z1)) +
  scale_colour_gradient2()


# Use your own colour scale with gradientn
ggplot(df, aes(x, y)) +
  geom_point(aes(colour = z1)) +
  scale_colour_gradientn(colours = terrain.colors(10))


# Equivalent fill scales do the same job for the fill aesthetic
ggplot(faithfuld, aes(waiting, eruptions)) +
  geom_raster(aes(fill = density)) +
  scale_fill_gradientn(colours = terrain.colors(10))


# Adjust colour choices with low and high
ggplot(df, aes(x, y)) +
  geom_point(aes(colour = z2)) +
  scale_colour_gradient(low = "white", high = "black")

# Avoid red-green colour contrasts because ~10% of men have difficulty
# seeing them

# Use `na.value = NA` to hide missing values but keep the original axis range
ggplot(df_na, aes(x = value, y)) +
  geom_bar(aes(fill = z1), stat = "identity") +
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)


ggplot(df_na, aes(x, y)) +
  geom_point(aes(colour = z1)) +
  scale_colour_gradient(low = "yellow", high = "red", na.value = NA)


##############
### try to plot ORF length (x-axis) and b value from wald test between WT polysome vs WT monosome

krausher_wtpoly_wtmono_wald_test <- read.csv("C:/Users/queenie.tsang/Desktop/CSDE1/Kraushar et al 2014 dataset/results/sleuth_wald_test_WTpoly_vs_WTmono_protein_coding.csv")

genes$Gene.name <- toupper(genes$Gene.name)

#merge the Genes table and the krausher_wtpoly_wtmono_wald_test df
krausher_WTpoly_WTmono_waldtest_genes <- merge(krausher_wtpoly_wtmono_wald_test, genes, by.x = "target_id", by.y = "Gene.name")

#convert the Transcript length from numeric to factor type
krausher_WTpoly_WTmono_waldtest_genes$Transcript.length..including.UTRs.and.CDS. <- as.factor(krausher_WTpoly_WTmono_waldtest_genes$Transcript.length..including.UTRs.and.CDS.)


#plot transcript length on the x-axis and the wald test bvalue for krausher (2014) WTpoly vs WTmono on y-axis
ggplot(krausher_WTpoly_WTmono_waldtest_genes, aes(x=Transcript.length..including.UTRs.and.CDS., y= b)) +
  geom_point()+
  labs(x= "Transcript.length..including.UTRs.and.CDS.", 
       y = "Wald Test B value Wild Type Polysome vs Wild Type Monosome",
       title= "Krausher et al (2014)")
