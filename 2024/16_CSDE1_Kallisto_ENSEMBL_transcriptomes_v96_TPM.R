#####  Get the ENSEMBL CDS length for  CSDE1 aligned with  Kallisto and ENSEMBL transcriptomes v96

setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/")

#install the bioconductor BiomaRt R package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt", force = TRUE)

BiocManager::install(version = "3.18")
devtools::install_version("dbplyr", version = "2.3.4")

library(biomaRt)

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
