#### Re-do CD1 polysome seq and CSDE1 targets overlay scatterplot
#### results from kallisto alignment with ENSEMBL GRCm38.p6, tximport and DESEQ2
## June 7 2024
# Queenie Tsang

##use only set 2 and set 3 of total, poly, mono from CD1 dataset
setwd("C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/")

#functions for processing data:
source("new_Celf2_polysome_seq_functions_feb_2024.R")

library(dplyr)
library(ggplot2)
library(readxl)
library(ggrepel)
library(stringr)
library(tidyverse)

#read in the table for CD1 Poly vs Total and Poly vs Mono table, filtered for protein coding genes
CD1_protein_coding_genes<-read.csv("results_kallisto_ENSEMBL_transcriptomes_v96/CD1_protein_coding_genes.csv")

#read in the table for CSDE1 RIP table 
CSDE1_protein_coding_genes <- read.csv("results_kallisto_ENSEMBL_transcriptomes_v96/CSDE1_protein_coding_genes.csv")

### for the CSDE1 table filter for CSDE1 targets; 2.5 fold change for CSDE1 RIP vs IgG RIP = log2FC 1.32
CSDE1_targets_from_CSDE1_vs_IgG <- CSDE1_protein_coding_genes[CSDE1_protein_coding_genes$CSDE1_vs_IgG_log2FoldChange >= 1.32, ]

#merge the CD1 table with the CSDE1 vs IgG table
merged_CD1_CSDE1_vs_IgG<- merge(CD1_protein_coding_genes, CSDE1_targets_from_CSDE1_vs_IgG, by.x="gene_ID", by.y="gene_ID", all = TRUE)

### save merged table of CD1 poly vs total and poly vs mono; CSDE1 vs IgG targets
write.csv(merged_CD1_CSDE1_vs_IgG, "C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/results_kallisto_ENSEMBL_transcriptomes_v96/merged_CD1_CSDE1_vs_IgG.csv")


#this first plots all genes in light grey in a scatterplot x-axis WTpoly-WTtotal, and y-axis WTpoly-WTmono
cd1_csde1_scatterplot4<- function(merged_table, title){
  p <- ggplot(data = merged_table,               ## this alters the order in which the points are plotted, so that CSDE1 plots are on top
              aes(x=CD1_poly_vs_mono_log2FoldChange, y=CD1_poly_vs_total_log2FoldChange, alpha=0.5, label = GENE_NAME.x)) +
    
    #label select genes that are layer marker genes
    #geom_label(data=subset(merged_table, target_id %in% c("TBR1", "SATB2", "TLE4", "POU3F3", "FEZF2", "FEZF1")))+
    
    #geom_point(stat="identity", size = 5, colour="lightgrey", show.legend = FALSE) +
    geom_point(stat="identity", size = 1, colour="darkgrey", show.legend = FALSE) +
    
    labs(x=("Log 2 Fold Change WT poly vs WT total"),
         y=("Log 2 Fold Change WT poly vs WT mono"),
         
         title=(title)) +
    theme_classic() # blank white background
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #      panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  plot(p)
}

#first plot all the CD1 genes using the log2FC value from DESEQ2
cd1_csde1_scatterplot4(merged_table = merged_CD1_CSDE1_vs_IgG, title = "CD1 E17 polysome")

#plot the CSDE1 target genes on top of this CD1 polysome seq E17 plot
colnames(merged_CD1_CSDE1_vs_IgG)

#order the rows of the  merged table by the CSDE1 vs IgG b values:
merged_CD1_CSDE1_vs_IgG_sort_log2FC<- merged_CD1_CSDE1_vs_IgG[order(merged_CD1_CSDE1_vs_IgG$`CSDE1_vs_IgG_log2FoldChange`),]

legend_title<-("Log2FC value CSDE1 vs IgG")

colnames(merged_CD1_CSDE1_vs_IgG_sort_log2FC)

cd1_csde1_scatterplot4(merged_table=merged_CD1_CSDE1_vs_IgG_sort_log2FC, title="CD1 Polysomeseq E17")+
  #geom_point(data = subset(merged_CD1_CSDE1_vs_IgG, Csde1.target.only.igg.b_iGG > 1.32 & Csde1.target.only.igg.qval_iGG <= 0.01),
  geom_point(data = subset(merged_CD1_CSDE1_vs_IgG_sort_log2FC, CSDE1_vs_IgG_log2FoldChange > 1.32 & CSDE1_vs_IgG_padj <= 0.001),
             
             stat = 'identity',
             shape=16,
             aes(x=CD1_poly_vs_total_log2FoldChange, y=CD1_poly_vs_mono_log2FoldChange, alpha=0.1, colour=CSDE1_vs_IgG_log2FoldChange, size = 0.5)) + #shape16 should not have outline around the points 
  
  #scale_colour_gradient(low= "cyan4", high="coral", legend_title)+
  
  scale_colour_gradient2(low="black", mid="purple", high="red", legend_title)+
  
  #black border outline only around marker genes
  #  geom_point(data=subset(merged_CD1_CSDE1_vs_IgG_sort_b, target_id %in% c("TBR1", "SATB2", "TLE4", "POU3F3", "FEZF2", "BCL11B")), aes(colour = Csde1.target.only.igg.b_iGG, size = 5))+
  
  #label select genes that are layer marker genes
  #geom_label(data=subset(merged_celf2_csde1_protein_coding_sort_b, target_id %in% c("TBR1", "SATB2", "TLE4", "POU3F3", "FEZF2", "FEZF1")))+
  
  # geom_label_repel(data=subset(merged_CD1_CSDE1_vs_IgG_sort_b, target_id %in% c("TBR1", "SATB2", "TLE4", "POU3F3", "FEZF2", "FEZF1", "BCL11B")), 
  #                  aes(label = target_id), 
  #                  segment.color= "black", colour = "black",
  #                  arrow = arrow(length = unit(0.25, 'cm'), type = 'closed'),
  #                  box.padding =  0.5,
  #                  max.overlaps = Inf)+
  # 
  labs(x=("LOG 2 Fold Change WT poly vs WT total"),
       y=("LOG 2 Fold Change WT poly vs WT mono"),
       title=("CD1_Polysomeseq E17 - filtered for protein coding genes; Log 2 Fold Change value cutoff > 1.32, qval <= 0.001")) +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
