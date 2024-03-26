# Create sample data ===================================================
set.seed(43)
data <- matrix(rnorm(500), 50, 10)
colnames(data) <- paste0("Sample_", 1:10)
rownames(data) <- paste0("Gene_", 1:50)
head(data)

# create a data frame for column annotation
ann_df <- data.frame(Group = rep(c("Disease", "Control"), c(5, 5)),
                     Lymphocyte_count = rnorm(10))
row.names(ann_df) <- colnames(data)
head(ann_df)
gene_functions_df <- data.frame(gene_functions = rep(c('Oxidative_phosphorylation', 
                                                       'Cell_cycle',
                                                       'Immune_regulation',
                                                       'Signal_transduction',
                                                       'Transcription'), rep(10, 5)))
row.names(gene_functions_df) <- rownames(data)
ann_colors <- list(
  gene_functions = c("Oxidative_phosphorylation" = "#F46D43",
                     "Cell_cycle" = "#708238",
                     "Immune_regulation" = "#9E0142",
                     "Signal_transduction" = "beige", 
                     "Transcription" = "violet"), 
  Group = c("Disease" = "darkgreen",
            "Control" = "blueviolet"),
  Lymphocyte_count = brewer.pal(5, 'PuBu')
)
