############# Kallisto files for CSDE1 RIP aligned with Kallisto and reference transcriptome from ENSEMBL transcriptiomes v96

######### try importing the Kallisto files using tximportData package
#https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto_with_TSV_files

library(tximportData)
library(tximport)
library(readr)
library(tibble)
library(DESeq2)

setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/")

samples <- read.table("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/CSDE1_kallisto_ENSEMBL_transcriptomes_v96_metadata.txt", header = TRUE, fill=TRUE)

### import the transcripts IDs and Gene IDs 
tx2gene <- read.table("tx2gene_ENSEMBL_transcriptomes_v96.txt", header=TRUE)
tx2gene <- as_tibble(tx2gene)

dir <- ("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/kallisto_ENSEMBL_transcriptomes_v96/")


### import kallisto with TSV files
samples$files <- paste0(dir, samples$kallisto, "/abundance.tsv")

write.csv(samples, "CSDE1_kallisto_with_ENSEMBL_transcriptome_v96_metadata.csv")

files<-file.path(samples$files)
names(files) <-paste0(samples$sample)


txi.kallisto.tsv <- tximport(files, type="kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
head(txi.kallisto.tsv$counts)



####  https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#transcript-abundance-files-and-tximport-tximeta
#https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#Downstream_DGE_in_Bioconductor

#There are two suggested ways of importing estimates for use with differential gene expression (DGE) methods. 
#The first method, which we show below for edgeR and for DESeq2, is to use the gene-level estimated counts from the 
#quantification tools, and additionally to use the transcript-level abundance estimates to calculate a gene-level offset
#that corrects for changes to the average transcript length across samples. The code examples below accomplish these steps
#for you, keeping track of appropriate matrices and calculating these offsets. For edgeR you need to assign a matrix 
#to y$offset, but the function DESeqDataSetFromTximport takes care of creation of the offset for you. Let’s call this 
#method “original counts and offset”.

#https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#transcript-abundance-files-and-tximport-tximeta
#Note that the tximport-to-DESeq2 approach uses estimated gene counts from the transcript abundance quantifiers, but not normalized counts.

#https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#tximport

sampleTable <- data.frame(condition = samples$state)
rownames(sampleTable) <- colnames(txi.kallisto.tsv$counts)

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, sampleTable, ~condition)

### double check that the order of the samples from sampleTable matches the samples order for txi.kallisto.tsv.CSDE1.IgG
sampleTable
# 
#         condition
# Input1    Capture
# Input2    Capture
# Input3    Capture
# IgG1          IgG
# IgG2          IgG
# IgG3          IgG
# Csde1_1      Csde
# Csde1_2      Csde
# Csde1_3      Csde


head(txi.kallisto.tsv$counts, 2)

#                       Input1 Input2 Input3 IgG1 IgG2 IgG3 Csde1_1 Csde1_2 Csde1_3
# ENSMUSG00000000001.4    2609   2605   2434  411  489  356     886    2199     586
# ENSMUSG00000000003.15      0      0      0    0    0    0       0       0       0

#### Prefiltering 

# While it is not necessary to pre-filter low count genes before running the DESeq2 functions, 
# there are two reasons which make pre-filtering useful: by removing rows in which there are very few reads, 
# we reduce the memory size of the dds data object, and we increase the speed of count modeling within DESeq2. 
# It can also improve visualizations, as features with no information for differential expression are not plotted 
# in dispersion plots or MA-plots.
# 
# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples.
# The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation for the minimal number of samples is to
# specify the smallest group size, e.g. here there are 3 treated samples. If there are not discrete groups, one can 
# use the minimal number of samples where non-zero counts would be considered interesting. One can also omit this step
# entirely and just rely on the independent filtering procedures available in results(), either IHW or genefilter. See 
# independent filtering section.

## For this analysis we will use 3 for the smallest group size, since each condition has 3 samples
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]





######## use Contrasts for pairwise comparisons

# A contrast is a linear combination of estimated log2 fold changes, which can be used to test if differences between groups are equal to zero. 
#The simplest use case for contrasts is an experimental design containing a factor with three levels, say A, B and C. Contrasts enable the user 
#to generate results for all 3 possible differences: log2 fold change of B vs A, of C vs A, and of C vs B. The contrast argument of results function 
#is used to extract test results of log2 fold changes of interest, for example:
#   
#   results(dds, contrast=c("condition","C","B"))
# 
# Log2 fold changes can also be added and subtracted by providing a list to the contrast argument which has two elements: the names of the log2 fold 
#changes to add, and the names of the log2 fold changes to subtract. The names used in the list should come from resultsNames(dds). Alternatively,
#a numeric vector of the length of resultsNames(dds) can be provided, for manually specifying the linear combination of terms. A tutorial describing 
#the use of numeric contrasts for DESeq2 explains a general approach to comparing across groups of samples. Demonstrations of the use of contrasts for 
#various designs can be found in the examples section of the help page ?results. The mathematical formula that is used to generate the contrasts can 
#be found below.


#### on the factor levels 
#you can either explicitly tell results which comparison to make using the contrast argument (this will be shown later), 
#or you can explicitly set the factors levels. In order to see the change of reference levels reflected in the results names, 
#you need to either run DESeq or nbinomWaldTest/nbinomLRT after the re-leveling operation. Setting the factor levels can
#be done in two ways, either using factor:

### specify the reference level, for the first comparison, set IgG as the reference: 
dds$condition <- relevel(dds$condition, ref = "IgG")

dds$condition 
# [1] Capture Capture Capture IgG     IgG     IgG     Csde    Csde    Csde   
# Levels: IgG Capture Csde

dds <- DESeq(dds)

# Differential expression analysis
# 
# The standard differential expression analysis steps are wrapped into a single function, DESeq. The estimation steps 
#performed by this function are described below, in the manual page for ?DESeq and in the Methods section of the DESeq2 
#publication (Love, Huber, and Anders 2014).
# 
# Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values
#and adjusted p values. With no additional arguments to results, the log2 fold change and Wald test p value will be for
#the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable 
#over the reference level (see previous note on factor levels). However, the order of the variables of the design do not
#matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of 
#results.
# 
# Details about the comparison are printed to the console, directly above the results table. The text, condition treated 
#vs untreated, tells you that the estimates are of the logarithmic fold change log2(treated/untreated).

CSDE1_vs_IgG <-results(dds, contrast=c("condition","Csde","IgG"))
CSDE1_vs_IgG

write.csv(CSDE1_vs_IgG, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/results/CSDE1_vs_IgG_DESEQ2.csv")


#### for the next comparison of CSDE1 vs Input
### specify the reference level, for the first comparison, set IgG as the reference: 
dds$condition <- relevel(dds$condition, ref = "Capture")

dds$condition 
# Capture Capture Capture IgG     IgG     IgG     Csde    Csde    Csde   
# Levels: Capture IgG Csde

dds <- DESeq(dds)
CSDE1_vs_Capture <- results(dds, contrast = c("condition", "Csde", "Capture"))


CSDE1_vs_Capture

# log2 fold change (MLE): condition Csde vs Capture 
# Wald test p-value: condition Csde vs Capture 
# DataFrame with 18393 rows and 6 columns
#                         baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#                       <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# ENSMUSG00000000001.4  1055.7763      -0.312959  0.238074  -1.31454 1.88663e-01 2.95959e-01
# ENSMUSG00000000028.15  166.8900      -1.648850  0.376437  -4.38015 1.18597e-05 6.47252e-05
# ENSMUSG00000000037.16  123.9243      -1.950859  0.474880  -4.10811 3.98913e-05 1.96213e-04
# ENSMUSG00000000056.7  1634.8435       0.350754  0.194823   1.80037 7.18022e-02 1.37498e-01
# ENSMUSG00000000058.6    87.3251      -1.026984  0.486719  -2.11002 3.48570e-02 7.61437e-02
# ...                         ...            ...       ...       ...         ...         ...
# ENSMUSG00000118332.1   602.0909       1.358585  0.292623  4.642783 3.43748e-06 2.07045e-05
# ENSMUSG00000118345.1    13.4581      -0.133917  0.903656 -0.148194 8.82189e-01 9.19285e-01
# ENSMUSG00000118346.1   133.3456      -0.861128  0.569831 -1.511200 1.30738e-01 2.20855e-01
# ENSMUSG00000118353.1    68.3191      -0.574606  0.907118 -0.633442 5.26445e-01 6.38923e-01
# ENSMUSG00000118382.1    11.3976      -1.713286  1.169552 -1.464908 1.42946e-01 2.37659e-01

#### save the results:
write.csv(CSDE1_vs_Capture, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/results/CSDE1_vs_Capture.csv")
