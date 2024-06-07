#### CD1 E17.5 Polysome seq 
#### Kallisto pseudoalignment done with ENSEMBL transcriptomes v96 which is GRCm38.p6

### Differential Expression analysis with tximport and DESEQ2

######### try importing the Kallisto files using tximportData package
#https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto_with_TSV_files

library(tximportData)
library(tximport)
library(readr)
library(tibble)
library(DESeq2)

setwd("C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/kallisto_ENSEMBL_transcriptomes_v96/")

samples <- read.table("CD1_E17.5_polysome_kallisto_ENSEMBL_transcriptomes_v96_metadata.txt", header=TRUE, fill=TRUE)

### for the  CD1 samples keep only set 2 and set 3 since set contained CD1-E17-1-Poly-S2 which was an outlier on the PCA plots
samples <- samples[4:9,]

#import the transcripts IDs and Gene IDs
tx2gene <- read.table("tx2gene_ENSEMBL_transcriptomes_v96.txt", header=TRUE)
tx2gene <- as_tibble(tx2gene)

dir <- ("C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/kallisto_ENSEMBL_transcriptomes_v96/")

### import kallisto files with TSV files 
samples$files <- paste0(dir, samples$kallisto, "/abundance.tsv")

write.csv(samples, "CD1_E17.5_polysomeseq_kallisto_with_ENSEMBL_transcriptome_v96_metadata.csv")

files <- file.path(samples$files)
names(files) <- paste0(samples$sample)

### import the kallisto files using tximport
#https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
# The tximport package has a single function for importing transcript-level estimates. The type argument is used
# to specify what software was used for estimation. A simple list with matrices, "abundance", "counts", and "length",
# is returned, where the transcript level information is summarized to the gene-level. Typically, abundance is provided by the 
# quantification tools as TPM (transcripts-per-million), while the counts are estimated counts (possibly fractional), and the "length"
# matrix contains the effective gene lengths. The "length" matrix can be used to generate an offset matrix for downstream gene-level differential 
# analysis of count matrices, as shown below.
#import the kallisto files using tximport

txi.kallisto.tsv <- tximport(files, type="kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
head(txi.kallisto.tsv$counts)


write.csv(txi.kallisto.tsv$counts, "C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/results_kallisto_ENSEMBL_transcriptomes_v96/CD1_E17.5_polysomeseq_kallisto_estimated_counts_tximport.csv")

head(txi.kallisto.tsv$abundance)
write.csv(txi.kallisto.tsv$abundance, "C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/results_kallisto_ENSEMBL_transcriptomes_v96/CD1_E17.5_polysomeseq_kallisto_abundance_tximport.csv")

####### save the abundance summarized to gene level, rearranged to have the different conditions together:
CD1_kallisto_abundance_rearranged <- txi.kallisto.tsv$abundance[,c("CD1-E17-2-Poly_S5", "CD1-E17-3-Poly_S8", "CD1-E17-2-Mono_S6", "CD1-E17-3-Mono_S9", "CD1-E17-2-Total_S4", "CD1-E17-3-Total_S7")]
colSums(CD1_kallisto_abundance_rearranged)
# CD1-E17-2-Poly_S5  CD1-E17-3-Poly_S8  CD1-E17-2-Mono_S6  CD1-E17-3-Mono_S9 CD1-E17-2-Total_S4 CD1-E17-3-Total_S7 
# 993752.7           993525.5           995388.8           994603.6           992567.0           991899.2 

write.csv(CD1_kallisto_abundance_rearranged, "C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/results_kallisto_ENSEMBL_transcriptomes_v96/CD1_set_2_set_3_kallisto_abundance.csv")

#https://compgenomr.github.io/book/gene-expression-analysis-using-high-throughput-sequencing-technologies.html
# "TPM also controls for both the library size and the gene lengths, however, with the TPM method, the read counts are first normalized 
# by the gene length (per kilobase), and then gene-length normalized values are divided by the sum of the gene-length normalized values 
# and multiplied by 10^6. Thus, the sum of normalized values for TPM will always be equal to 10^6 for each library, while the sum of
# RPKM/FPKM values do not sum to 10^6. Therefore, it is easier to interpret TPM values than RPKM/FPKM values."


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

sampleTable
#                     condition
# CD1-E17-2-Total_S4     total
# CD1-E17-2-Poly_S5       poly
# CD1-E17-2-Mono_S6       mono
# CD1-E17-3-Total_S7     total
# CD1-E17-3-Poly_S8       poly
# CD1-E17-3-Mono_S9       mono

head(txi.kallisto.tsv$counts, 2)
#                     #   CD1-E17-2-Total_S4 CD1-E17-2-Poly_S5 CD1-E17-2-Mono_S6 CD1-E17-3-Total_S7 CD1-E17-3-Poly_S8
# ENSMUSG00000000001.4                5827              7860              3109               5737              9040
# ENSMUSG00000000003.15                  0                 0                 0                  0                 0
#                         CD1-E17-3-Mono_S9
# ENSMUSG00000000001.4               3346
# ENSMUSG00000000003.15                 0                    
# #### Prefiltering 

# While it is not necessary to pre-filter low count genes before running the DESeq2 functions, 
# there are two reasons which make pre-filtering useful: by removing rows in which there are very few reads, 
# we reduce the memory size of the dds data object, and we increase the speed of count modeling within DESeq2. 
# It can also improve visualizations, as features with no information for differential expression are not plotted 
# in dispersion plots or MA-plots.
# 
# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples.
# The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation for the minimal number of samples is to
# specify the smallest group size, e.g. here there are 3 treated samples in the example.
#For the CD1 polysome dataset, each condition group also has 3 samples.

#If there are not discrete groups, one can 
# use the minimal number of samples where non-zero counts would be considered interesting. One can also omit this step
# entirely and just rely on the independent filtering procedures available in results(), either IHW or genefilter. See 
# independent filtering section.

## For this analysis we will use 2 for the smallest group size, since each condition has 2 samples
smallestGroupSize <- 2
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
dds$condition <- relevel(dds$condition, ref = "total")

dds$condition 
# total poly  mono  total poly  mono
# Levels: total mono poly

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

CD1_poly_vs_total <-results(dds, contrast=c("condition","poly","total"))

CD1_poly_vs_total


# log2 fold change (MLE): condition poly vs total 
# Wald test p-value: condition poly vs total 
# DataFrame with 17852 rows and 6 columns
#                       baseMean log2FoldChange     lfcSE      stat      pvalue       padj
#                         <numeric>      <numeric> <numeric> <numeric>   <numeric>  <numeric>
# ENSMUSG00000000001.4   5690.977      0.3534846 0.0913009  3.871645 0.000108104 0.00090957
# ENSMUSG00000000028.15   278.729      0.2266814 0.2417069  0.937836 0.348328732 0.54082707
# ENSMUSG00000000037.16   280.510     -0.6581650 0.3015885 -2.182328 0.029085336 0.09085152
# ENSMUSG00000000056.7   5959.063      0.0724718 0.0911992  0.794654 0.426814664 0.61473056
# ENSMUSG00000000058.6    309.724     -0.0400676 0.2416521 -0.165807 0.868308858 0.92931457
# ...                         ...            ...       ...       ...         ...        ...
# ENSMUSG00000118318.1   159.9397     -1.0818002  0.311116 -3.477157 0.000506762  0.0033277
# ENSMUSG00000118332.1  2944.9711     -0.0530707  0.126393 -0.419887 0.674567913  0.8063105
# ENSMUSG00000118345.1    12.0240     -0.1908204  0.993185 -0.192130 0.847640582         NA
# ENSMUSG00000118346.1   282.1591     -0.5989453  0.232622 -2.574762 0.010030906  0.0393930
# ENSMUSG00000118382.1    10.7981     -2.1163668  1.284906 -1.647099 0.099537686         NA

write.csv(CD1_poly_vs_total, "C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/results_kallisto_ENSEMBL_transcriptomes_v96/CD1_poly_vs_total_DESEQ2.csv")

CD1_poly_vs_total_df <-data.frame(CD1_poly_vs_total)

CD1_poly_vs_total_df$gene_ID <- row.names(CD1_poly_vs_total_df)

### keep just Gene ID and Gene Name columns from tx2gene file
tx2gene_keep <-tx2gene[,c("GENE_ID", "GENE_NAME")]

#remove duplicate rows
tx2gene_keep <- unique(tx2gene_keep)

## add the gene name /symbol from tx2 gene file
CD1_poly_vs_total_merged <- merge(CD1_poly_vs_total_df, tx2gene_keep, by.x="gene_ID", by.y="GENE_ID")

colnames(CD1_poly_vs_total_merged)

CD1_poly_vs_total_keep <- CD1_poly_vs_total_merged[,c("gene_ID", "log2FoldChange", "pvalue", "padj", "GENE_NAME")]
colnames(CD1_poly_vs_total_keep) <- c("gene_ID", "CD1_poly_vs_total_log2FoldChange", "CD1_poly_vs_total_pvalue", "CD1_poly_vs_total_padj", "GENE_NAME")

write.csv(CD1_poly_vs_total_keep, "C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/results_kallisto_ENSEMBL_transcriptomes_v96/CD1_poly_vs_total.csv")



########### for the comparison of CD1 Poly vs Mono
####### specify the reference level, set Mono as the reference
dds$condition <- relevel(dds$condition, ref="mono")

dds$condition
# total poly  mono  total poly  mono 
# Levels: mono total poly

dds <- DESeq(dds)

CD1_poly_vs_mono <- results(dds, contrast = c("condition", "poly", "mono"))
CD1_poly_vs_mono

# log2 fold change (MLE): condition poly vs mono 
# Wald test p-value: condition poly vs mono 
# DataFrame with 17852 rows and 6 columns
#                       baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#                       <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# ENSMUSG00000000001.4   5690.977       1.152113 0.0920830 12.511683 6.44439e-36 6.93461e-35
# ENSMUSG00000000028.15   278.729       0.569766 0.2451209  2.324427 2.01026e-02 3.12824e-02
# ENSMUSG00000000037.16   280.510      -1.721835 0.2974736 -5.788193 7.11475e-09 2.22986e-08
# ENSMUSG00000000056.7   5959.063       0.657419 0.0916854  7.170377 7.47917e-13 3.01055e-12
# ENSMUSG00000000058.6    309.724       0.100726 0.2435523  0.413571 6.79189e-01 7.29447e-01
# ...                         ...            ...       ...       ...         ...         ...
# ENSMUSG00000118318.1   159.9397      -2.539421  0.301482  -8.42312 3.66603e-17 1.85190e-16
# ENSMUSG00000118332.1  2944.9711      -0.476085  0.126152  -3.77391 1.60709e-04 3.40976e-04
# ENSMUSG00000118345.1    12.0240      -1.818639  0.941968  -1.93068 5.35227e-02 7.64513e-02
# ENSMUSG00000118346.1   282.1591      -2.124463  0.224528  -9.46192 3.02327e-21 1.83889e-20
# ENSMUSG00000118382.1    10.7981      -3.389668  1.258535  -2.69334 7.07391e-03 1.19768e-02

write.csv(CD1_poly_vs_mono, "C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/results_kallisto_ENSEMBL_transcriptomes_v96/CD1_poly_vs_mono.csv")

CD1_poly_vs_mono

# log2 fold change (MLE): condition poly vs mono 
# Wald test p-value: condition poly vs mono 
# DataFrame with 17852 rows and 6 columns
#                       baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#                       <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   ENSMUSG00000000001.4   5690.977       1.152113 0.0920830 12.511683 6.44439e-36 6.93461e-35
# ENSMUSG00000000028.15   278.729       0.569766 0.2451209  2.324427 2.01026e-02 3.12824e-02
# ENSMUSG00000000037.16   280.510      -1.721835 0.2974736 -5.788193 7.11475e-09 2.22986e-08
# ENSMUSG00000000056.7   5959.063       0.657419 0.0916854  7.170377 7.47917e-13 3.01055e-12
# ENSMUSG00000000058.6    309.724       0.100726 0.2435523  0.413571 6.79189e-01 7.29447e-01
# ...                         ...            ...       ...       ...         ...         ...
# ENSMUSG00000118318.1   159.9397      -2.539421  0.301482  -8.42312 3.66603e-17 1.85190e-16
# ENSMUSG00000118332.1  2944.9711      -0.476085  0.126152  -3.77391 1.60709e-04 3.40976e-04
# ENSMUSG00000118345.1    12.0240      -1.818639  0.941968  -1.93068 5.35227e-02 7.64513e-02
# ENSMUSG00000118346.1   282.1591      -2.124463  0.224528  -9.46192 3.02327e-21 1.83889e-20
# ENSMUSG00000118382.1    10.7981      -3.389668  1.258535  -2.69334 7.07391e-03 1.19768e-02

CD1_poly_vs_mono_df <- data.frame(CD1_poly_vs_mono)
CD1_poly_vs_mono_df$gene_ID <- row.names(CD1_poly_vs_mono_df)

#merge the text2gene file to add the gene names to this dataframe
CD1_poly_vs_mono_merged <- merge(CD1_poly_vs_mono_df, tx2gene_keep, by.x = "gene_ID", by.y="GENE_ID")
colnames(CD1_poly_vs_mono_merged)

##keep just the columns for Log2FC 
CD1_poly_vs_mono_keep <- CD1_poly_vs_mono_merged[,c("gene_ID", "log2FoldChange", "pvalue", "padj", "GENE_NAME")]

#rename the column names for CD1 poly vs mono 
colnames(CD1_poly_vs_mono_keep) <- c("gene_ID", "CD1_poly_vs_mono_log2FoldChange", "CD1_poly_vs_mono_pvalue", "CD1_poly_vs_mono_padj", "GENE_NAME")

#merge the CD1 poly vs total and CD1 poly vs mono tables together
CD1_dataframe <- merge(CD1_poly_vs_mono_keep, CD1_poly_vs_total_keep, by.x="gene_ID", by.y="gene_ID")

colnames(CD1_dataframe)

###keep the relevant columns
CD1_df_keep <- CD1_dataframe[,c("gene_ID", "CD1_poly_vs_mono_log2FoldChange", "CD1_poly_vs_mono_pvalue",      
  "CD1_poly_vs_mono_padj", "CD1_poly_vs_total_log2FoldChange", "CD1_poly_vs_total_pvalue", "CD1_poly_vs_total_padj", "GENE_NAME.y")]

#rename the last Gene names column to gene name 
colnames(CD1_df_keep)[8]<-"GENE_NAME"

write.csv(CD1_df_keep, "C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/results_kallisto_ENSEMBL_transcriptomes_v96/CD1_DESEQ2_all_genes.csv")

### filter for protein coding genes only
setwd("C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/kallisto_ENSEMBL_transcriptomes_v96/")

protein_coding_genes <- read.csv("protein_coding_genes_biomart_ensembl.txt")

CD1_df_keep$GENE_NAME <- toupper(CD1_df_keep$GENE_NAME)

#convert the protein coding genes to capital letters to better filter for CD1 
protein_coding_genes$protein_caps <- toupper(protein_coding_genes$Gene.name)

#filter for just the protein coding genes in the CD1 dataframe:
CD1_protein_coding_genes <-CD1_df_keep[CD1_df_keep$GENE_NAME %in% protein_coding_genes$protein_caps,]

## save the results
write.csv(CD1_protein_coding_genes, "C:/Users/queenie.tsang/Desktop/wt_CD1_E17_cortex_polysome_feb_15_2024/results_kallisto_ENSEMBL_transcriptomes_v96/CD1_protein_coding_genes.csv")


