#####  Get the ENSEMBL CDS length for  CSDE1 aligned with  Kallisto and ENSEMBL transcriptomes v96

setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/")

#install the bioconductor BiomaRt R package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt", force = TRUE)

BiocManager::install(version = "3.18")
devtools::install_version("dbplyr", version = "2.3.4")

library(biomaRt)
library(GenomicFeatures)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

listEnsembl()

#https://useast.ensembl.org/info/data/biomart/biomart_r_package.html

ensembl = useEnsembl(biomart="ensembl")

(listDatasets(ensembl))

list(listDatasets(ensembl))

GRCm38=useEnsembl(biomart="ensembl", GRCm=38)

grch37 = useEnsembl(biomart="ensembl",GRCh=37)
listDatasets(grch37)[31:35,]

library(biomaRt)

ensembl=useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
head(listFilters(ensembl))
#             name              description
# 1 chromosome_name Chromosome/scaffold name
# 2           start                    Start
# 3             end                      End
# 4      band_start               Band Start
# 5        band_end                 Band End

listDatasets(grcm38)
# 6          strand                   Strand

#To USE mouse Genome GRCm38 for biomart R
#https://www.biostars.org/p/9550266/

#https://useast.ensembl.org/info/website/archives/assembly.html

#ENSEMBL transcriptomes v96 corresponds to mouse GRCm38.p6, which was use from version 92 to version 102 
ensembl96 <- useEnsembl(biomart = 'genes', dataset="mmusculus_gene_ensembl", version = 96)  ### 
# Error: Specified Ensembl version is not available.
# Use listEnsemblArchives() to view available versions.

listEnsemblArchives()
ensembl102 <- useEnsembl(biomart = 'genes', dataset="mmusculus_gene_ensembl", version = 102)

#The "listFilters" function will give you the list of available filters for a given mart and species:
listFilters(ensembl102)

#filters: The "listAttributes" function will give you the list of the available attributes for a given mart and species:
listAttributes(ensembl102)

# ensembl_gene_id
# ensembl_transcript_id
# 5_utr_start
# 5_utr_end
# cds_length

#The "getBM" function allow you to build a BioMart query using a list of mart filters and attributes.

ENSEMBL_GRCm38.p6_genes <- getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', '5_utr_start', '5_utr_end', '3_utr_start', '3_utr_end', 'cds_length'),
                                 mart=ensembl102)

colnames(ENSEMBL_GRCm38.p6_genes)
ENSEMBL_GRCm38.p6_genes$`5_utr_length` <- ENSEMBL_GRCm38.p6_genes$`5_utr_end` - ENSEMBL_GRCm38.p6_genes$`5_utr_start`

ENSEMBL_GRCm38.p6_genes$`3_utr_length` <- ENSEMBL_GRCm38.p6_genes$`3_utr_end` - ENSEMBL_GRCm38.p6_genes$`3_utr_start`
  
### save the results: 
write.csv(ENSEMBL_GRCm38.p6_genes, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/results/ENSEMBL_GRCm38.p6_genes.csv")

ENSEMBL_GRCm38.p6_genes <- read.csv("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/results/ENSEMBL_GRCm38.p6_genes.csv")

#### filter ENSEMBL table for protein coding genes:
## read in the list of protein coding genes downloaded from ENSEMBL Biomart tool
protein_coding_genes <- read.csv("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/protein_coding_genes_biomart_ensembl.txt")

ENSEMBL_GRCm38.p6_protein_coding_genes <- ENSEMBL_GRCm38.p6_genes[ENSEMBL_GRCm38.p6_genes$external_gene_name %in% protein_coding_genes$Gene.name, ]

### read in the CSDE1 protein coding genes DESEQ2 results:
CSDE1_protein_coding_genes <- read.csv("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/results/CSDE1_protein_coding_genes.csv")

#Total CDS length is the sum of all constitutive and alternately spliced exons; average CDS length is the average CDS length of each transcript

#### for the gene version of this table: keep the row (transcript) with the longest CDS length
#https://support.bioconductor.org/p/9140688/

##### Extract the longest CDS for each gene 
cds_length_per_gene <- splitAsList(ENSEMBL_GRCm38.p6_protein_coding_genes$cds_length, ENSEMBL_GRCm38.p6_protein_coding_genes$ensembl_gene_id)

# IntegerList of length 6
# [["ENSMUSG00000000001"]] 1065 1065 1065 1065
# [["ENSMUSG00000000003"]] 525 525 525 525 414 414 414 414
# [["ENSMUSG00000000028"]] 1701 1701 1701 1701 1701 1563 1563 1563 1563 410 410 <NA>
# [["ENSMUSG00000000037"]] 1101 1101 2478 2478 2478 2484 2484 2484 2538 ... 2910 <NA> 392 392 2706 2706 2706 2706
# [["ENSMUSG00000000049"]] 1038 1038 1038 154 154 154 99 99 99 464 464 464
# [["ENSMUSG00000000056"]] 1389 1389 1389 <NA> <NA>

transcript_id_per_gene <- splitAsList(ENSEMBL_GRCm38.p6_protein_coding_genes$ensembl_transcript_id, ENSEMBL_GRCm38.p6_protein_coding_genes$ensembl_gene_id)

transcript_id_per_gene

which_max <- which.max(cds_length_per_gene)
length_of_longest_cds <- unlist(cds_length_per_gene[as.list(which_max)])

#length of longest cds is a named integer vector that maps each gene id to the length of the longest CDS for that gene:
head(length_of_longest_cds)
# ENSMUSG00000000001 ENSMUSG00000000003 ENSMUSG00000000028 ENSMUSG00000000037 ENSMUSG00000000049 ENSMUSG00000000056 
# 1065                525               1701               2910               1038               1389 

#'transcript id of longest cds' is a named integer vector that maps each gene id to the transcript id of the longest CDS for that gene
transcript_id_of_longest_cds <- unlist(transcript_id_per_gene[as.list(which_max)])

head(transcript_id_of_longest_cds)
# ENSMUSG00000000001   ENSMUSG00000000003   ENSMUSG00000000028   ENSMUSG00000000037   ENSMUSG00000000049 
# "ENSMUST00000000001" "ENSMUST00000000003" "ENSMUST00000000028" "ENSMUST00000077375" "ENSMUST00000000049" 
# ENSMUSG00000000056 
# "ENSMUST00000103015"

### 'length of longest cds' and 'transcript_id_of_longest length are parallel vectors with identical names:
identical(names(length_of_longest_cds), names(transcript_id_of_longest_cds))  #TRUE

#summarize the results in a 3-column dataframe with 1 row per gene:
longest_cds <- data.frame(gene_id=names(length_of_longest_cds), length_of_longest_cds=length_of_longest_cds, transcript_id_of_longest_cds=transcript_id_of_longest_cds)

#let's remove rows for which the longest CDS has length 0 (non-coding genes):
longest_cds <- subset(longest_cds, length_of_longest_cds !=0 )

write.csv(longest_cds, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/results/ENSEMBL_GRCm38.p6_longest_cds_for_protein_coding_genes.csv")


########## for the ENSEMBL GRCm38.p6 table keep just the rows where transcript ID is for the longest CDS
ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds <- ENSEMBL_GRCm38.p6_protein_coding_genes[ENSEMBL_GRCm38.p6_protein_coding_genes$ensembl_transcript_id %in% longest_cds$transcript_id_of_longest_cds,]

#for ENSEMBL GRCm38.p6 table, keep just the gene_id, transcript_id, gene name, cds_length, 5' UTR length and 3' UTR length columns
colnames(ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds)
ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_keep <- ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds[,c("ensembl_gene_id", "ensembl_transcript_id","external_gene_name", "cds_length", "X5_utr_length", "X3_utr_length")]

library(tidyr)
library(dplyr)

#### tidy the table so that rows that belong to the same transcript id gather into 1 row; if there's multiple values, then have them comma separated
ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_keep2 = ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_keep %>% 
  group_by(ensembl_transcript_id) %>%
  summarise_all(~ toString(na.omit(.)))

#remove the duplicate values for ensembl gene id, and external gene name 
ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_keep2$ensembl_gene_id <- gsub("^(.*?),.*", "\\1", ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_keep2$ensembl_gene_id)

ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_keep2$external_gene_name <-  gsub("^(.*?),.*", "\\1", ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_keep2$external_gene_name)
ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_keep2$cds_length <-  gsub("^(.*?),.*", "\\1", ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_keep2$cds_length)  

### save the results
write.csv(ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_keep2, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/results/ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_5UTR_3UTR_length.csv")


#### merge the CSDE1 protein coding genes table with the ENSEMBL GRCm38.p6 protein coding table
ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_keep2$external_gene_name <- toupper(ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_keep2$external_gene_name)

CSDE1_protein_coding_genes_ENSEMBL_merged<-merge(CSDE1_protein_coding_genes, ENSEMBL_GRCm38.p6_protein_coding_genes_longest_cds_keep2, by.x="GENE_NAME", by.y="external_gene_name")

colnames(CSDE1_protein_coding_genes_ENSEMBL_merged)

## keep columns 
CSDE1_protein_coding_genes_ENSEMBL_merged_keep <- CSDE1_protein_coding_genes_ENSEMBL_merged[,c("GENE_NAME", "gene_ID", "CSDE1_vs_IgG_log2FoldChange", "CSDE1_vs_IgG_pvalue", 
                                                                                                    "CSDE1_vs_IgG_padj",  "CSDE1_vs_Capture_log2FoldChange", "CSDE1_vs_Capture_pvalue",
                                                                                                    "CSDE1_vs_Capture_padj", "ensembl_transcript_id", "cds_length", 
                                                                                                    "X5_utr_length", "X3_utr_length")]
#count the number of NA rows in the data frame
colSums(is.na(CSDE1_protein_coding_genes_ENSEMBL_merged_keep))
#there are 28 genes which are missing CSDE1_vs_IgG_pvalue, CSDE1_vs_IgG_padj, CSDE1_vs_Capture_pvalue, and  CSDE1_vs_Capture_padj


## remove the genes that have NA for the CSDE1 vs IGG, and CSDE1 vs Capture columns
CSDE1_protein_coding_genes_ENSEMBL_merged_keep2 <- na.omit(CSDE1_protein_coding_genes_ENSEMBL_merged_keep)
dim(CSDE1_protein_coding_genes_ENSEMBL_merged_keep2)
#15104    12

dim(CSDE1_protein_coding_genes_ENSEMBL_merged_keep)
#15132    12

## save this table for CSDE1 protein coding genes 
write.csv(CSDE1_protein_coding_genes_ENSEMBL_merged_keep2, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/results/CSDE1_protein_coding_genes_ENSEMBL_longest_cds_transcript_5UTR_3UTR.csv")
