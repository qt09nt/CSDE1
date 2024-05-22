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


### import the kallisto files using tximport
#https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
# The tximport package has a single function for importing transcript-level estimates. The type argument is used
# to specify what software was used for estimation. A simple list with matrices, "abundance", "counts", and "length",
# is returned, where the transcript level information is summarized to the gene-level. Typically, abundance is provided by the 
# quantification tools as TPM (transcripts-per-million), while the counts are estimated counts (possibly fractional), and the "length"
# matrix contains the effective gene lengths. The "length" matrix can be used to generate an offset matrix for downstream gene-level differential 
# analysis of count matrices, as shown below.

#the argument "ignoreAfterBar" is for facilitating the summarization to gene level
txi.kallisto.tsv <- tximport(files, type="kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
head(txi.kallisto.tsv$counts)

write.csv(txi.kallisto.tsv$counts, "results/CSDE1_RIP_Kallisto_estimated_counts_tximport.csv")


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
# specify the smallest group size, e.g. here there are 3 treated samples in the example.
#For the CSDE1 RIP dataset, each condition group also has 3 samples.

#If there are not discrete groups, one can 
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


#### Note on the factor levels 
#you can either explicitly tell results which comparison to make using the contrast argument (this will be shown later), 
#or you can explicitly set the factors levels. In order to see the change of reference levels reflected in the results names, 
#you need to either run DESeq or nbinomWaldTest/nbinomLRT after the re-leveling operation. Setting the factor levels can
#be done in two ways, either using factor:

### specify the reference level, for the first comparison, set IgG as the reference: 
dds$condition <- relevel(dds$condition, ref = "IgG")

dds$condition 
# [1] Capture Capture Capture IgG     IgG     IgG     Csde    Csde    Csde   
# Levels: IgG Capture Csde

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

## from the ?results help page:
# contrast	
# this argument specifies what comparison to extract from the object to build a results table. one of either:
#   
#   a character vector with exactly three elements: the name of a factor in the design formula, 
#the name of the numerator level for the fold change, and the name of the denominator level for the fold change (simplest case)
# 
# a list of 2 character vectors: the names of the fold changes for the numerator, and the names of the fold changes for the 
#denominator. these names should be elements of resultsNames(object). if the list is length 1, a second element is added which 
#is the empty character vector, character(). (more general case, can be to combine interaction terms and main effects)
# 
# a numeric contrast vector with one element for each element in resultsNames(object) (most general case)
# 
# If specified, the name argument is ignored.


dds <- DESeq(dds)

CSDE1_vs_IgG <-results(dds, contrast=c("condition","Csde","IgG"))
CSDE1_vs_IgG

# log2 fold change (MLE): condition Csde vs IgG 
# Wald test p-value: condition Csde vs IgG 
# DataFrame with 18393 rows and 6 columns
#                       baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#                         <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# ENSMUSG00000000001.4  1055.7763       0.076160  0.240978  0.316046 7.51968e-01 8.26691e-01
# ENSMUSG00000000028.15  166.8900      -1.689625  0.383374 -4.407246 1.04693e-05 5.60035e-05
# ENSMUSG00000000037.16  123.9243      -0.697021  0.490326 -1.421548 1.55157e-01 2.49111e-01
# ENSMUSG00000000056.7  1634.8435       0.185037  0.196678  0.940814 3.46800e-01 4.65087e-01
# ENSMUSG00000000058.6    87.3251      -0.910550  0.497998 -1.828421 6.74864e-02 1.27359e-01
# ...                         ...            ...       ...       ...         ...         ...
# ENSMUSG00000118332.1   602.0909       0.440730  0.295026  1.493870    0.135210  0.22336342
# ENSMUSG00000118345.1    13.4581      -0.676406  0.935306 -0.723192    0.469562  0.58545958
# ENSMUSG00000118346.1   133.3456      -0.803971  0.576016 -1.395745    0.162791  0.25894063
# ENSMUSG00000118353.1    68.3191       2.991575  0.985474  3.035673    0.002400  0.00740811
# ENSMUSG00000118382.1    11.3976      -1.133974  1.213126 -0.934754    0.349915  0.46816062

write.csv(CSDE1_vs_IgG, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/results/CSDE1_vs_IgG_DESEQ2.csv")

#get the dataframe for CSDE1 vs IgG comparison
CSDE1_vs_IgG_df <- data.frame(CSDE1_vs_IgG)
CSDE1_vs_IgG_df$gene_ID <- row.names(CSDE1_vs_IgG_df)

#keep just the Gene ID and Gene Name columns from the tx2gene file
tx2gene_keep <-tx2gene[,c("GENE_ID", "GENE_NAME")]
#remove duplicate rows
tx2gene_keep <- unique(tx2gene_keep)

## add the gene name /symbol from tx2 gene file
CSDE1_vs_IgG_merged <- merge(CSDE1_vs_IgG_df, tx2gene_keep, by.x="gene_ID", by.y="GENE_ID")

#keep the CSDE1_vs_IgG_merged columns 
colnames(CSDE1_vs_IgG_merged)

CSDE1_vs_IGG_keep <- CSDE1_vs_IgG_merged[,c("gene_ID", "log2FoldChange", "pvalue", "padj", "GENE_NAME")]
colnames(CSDE1_vs_IGG_keep)<-c("gene_ID", "CSDE1_vs_IgG_log2FoldChange", "CSDE1_vs_IgG_pvalue", "CSDE1_vs_IgG_padj", "GENE_NAME")

#### for the next comparison of CSDE1 vs Input
### specify the reference level, for the first comparison, set IgG as the reference: 
dds$condition <- relevel(dds$condition, ref = "Capture")

dds$condition 
# Capture Capture Capture IgG     IgG     IgG     Csde    Csde    Csde   
# Levels: Capture IgG Csde

dds <- DESeq(dds)
CSDE1_vs_Capture <- results(dds, contrast = c("condition", "Csde", "Capture"))

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


#####
CSDE1_vs_Capture


CSDE1_vs_Capture_df <- data.frame(CSDE1_vs_Capture)
CSDE1_vs_Capture_df$gene_ID <- row.names(CSDE1_vs_Capture_df)

#merge the text2gene file to add the gene names to this dataframe
CSDE1_Capture_df_merged <- merge(CSDE1_vs_Capture_df, tx2gene_keep, by.x = "gene_ID", by.y="GENE_ID")
colnames(CSDE1_Capture_df_merged)

#keep just the columns for Log2FC
CSDE1_Capture_df_merged_keep <- CSDE1_Capture_df_merged[,c("gene_ID", "log2FoldChange", "pvalue", "padj", "GENE_NAME")]

#rename the column names for CSDE1 vs Capture dataframe:
colnames(CSDE1_Capture_df_merged_keep) <- c("gene_ID", "CSDE1_vs_Capture_log2FoldChange", "CSDE1_vs_Capture_pvalue", "CSDE1_vs_Capture_padj", "GENE_NAME")

#### merge the CSDE1 vs IgG and the CSDE1 vs Capture dataframes together
CSDE1_dataframe <- merge(CSDE1_vs_IGG_keep, CSDE1_Capture_df_merged_keep, by.x="gene_ID", by.y="gene_ID")
colnames(CSDE1_dataframe)

#keep just the GENE ID, "CSDE1_vs_IgG_log2FoldChange", "CSDE1_vs_IgG_pvalue", "CSDE1_vs_IgG_padj", "GENE_NAME.x", "CSDE1_vs_Capture_log2FoldChange","CSDE1_vs_Capture_pvalue", "CSDE1_vs_Capture_padj", "GENE_NAME.y"                    
CSDE1_dataframe_keep <- CSDE1_dataframe[,c("gene_ID", "CSDE1_vs_IgG_log2FoldChange", "CSDE1_vs_IgG_pvalue", "CSDE1_vs_IgG_padj", "CSDE1_vs_Capture_log2FoldChange", "CSDE1_vs_Capture_pvalue", "CSDE1_vs_Capture_padj", "GENE_NAME.y")]

colnames(CSDE1_dataframe_keep)[8] <- "GENE_NAME"

write.csv(CSDE1_dataframe_keep, "results/CSDE1_dataframe_keep.csv")


#### filter for protein coding genes only
setwd("C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/")

protein_coding_genes <- read.csv("protein_coding_genes_biomart_ensembl.txt")

CSDE1_dataframe_keep$GENE_NAME <- toupper(CSDE1_dataframe_keep$GENE_NAME)

#convert the protein coding genes to capital letters to better use as a filter for the CSDE1 dataframe
protein_coding_genes$protein_caps <- toupper(protein_coding_genes$Gene.name)

#####filter for just the protein coding genes in the CSDE1 dataframe
CSDE1_protein_coding_genes <- CSDE1_dataframe_keep[CSDE1_dataframe_keep$GENE_NAME %in% protein_coding_genes$protein_caps, ]

###### save the results 
write.csv(CSDE1_protein_coding_genes, "C:/Users/queenie.tsang/Desktop/CSDE1/Csde1RNA_IP/kallisto_ENSEMBL_transciptomes_v96/results/CSDE1_protein_coding_genes.csv")


### Log fold change shrinkage for visualization and ranking
# Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. To shrink the LFC, we pass the dds object to the function
# lfcShrink. Below we specify to use the apeglm method for effect size shrinkage (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.
# 
# We provide the dds object and the name or number of the coefficient we want to shrink, where the number refers to the order of the coefficient as
# it appears in resultsNames(dds).
resultsNames(dds)
# "Intercept"                 "condition_IgG_vs_Capture"  "condition_Csde_vs_Capture"

resLFC_CSDE1_vs_Capture <- lfcShrink(dds, coef="condition_Csde_vs_Capture", type="apeglm")
resLFC_CSDE1_vs_Capture

# log2 fold change (MAP): condition Csde vs Capture 
# Wald test p-value: condition Csde vs Capture 
# DataFrame with 18393 rows and 5 columns
#                         baseMean log2FoldChange     lfcSE      pvalue        padj
#                         <numeric>      <numeric> <numeric>   <numeric>   <numeric>
# ENSMUSG00000000001.4  1055.7763      -0.281810  0.227920 1.88663e-01 2.95959e-01
# ENSMUSG00000000028.15  166.8900      -1.515799  0.382639 1.18597e-05 6.47252e-05
# ENSMUSG00000000037.16  123.9243      -1.752110  0.491741 3.98913e-05 1.96213e-04
# ENSMUSG00000000056.7  1634.8435       0.329275  0.189578 7.18022e-02 1.37498e-01
# ENSMUSG00000000058.6    87.3251      -0.790380  0.473125 3.48570e-02 7.61437e-02
# ...                         ...            ...       ...         ...         ...
# ENSMUSG00000118332.1   602.0909      1.2762428  0.294984 3.43748e-06 2.07045e-05
# ENSMUSG00000118345.1    13.4581     -0.0497702  0.553895 8.82189e-01 9.19024e-01
# ENSMUSG00000118346.1   133.3456     -0.5739599  0.512250 1.30738e-01 2.20855e-01
# ENSMUSG00000118353.1    68.3191     -0.2173907  0.580614 5.26445e-01 6.38923e-01
# ENSMUSG00000118382.1    11.3976     -0.5240065  0.803744 1.42946e-01 2.37659e-01

### p values and adjusted p-values
#order the results by the smallest p value:
resOrdered_CSDE1_vs_Capture <- CSDE1_vs_Capture[order(CSDE1_vs_Capture$pvalue),]

summary(CSDE1_vs_Capture)

#How many adjusted p-values were less than 0.05?
sum(CSDE1_vs_Capture$padj < 0.05, na.rm=TRUE)
#7729

### exploring results:

##MA plot
# In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts 
# for all the samples in the DESeqDataSet. Points will be colored blue if the adjusted p value is less than 0.1. Points which fall
# out of the window are plotted as open triangles pointing either up or down.
plotMA(CSDE1_vs_Capture, ylim=c(-2,2))

#It is more useful to visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from 
#low count genes without requiring arbitrary filtering thresholds.
plotMA(resLFC_CSDE1_vs_Capture, ylim=c(-2,2))
