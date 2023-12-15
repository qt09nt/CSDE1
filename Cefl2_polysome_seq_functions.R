#functions for analysis of Celf2_KI_polysome seq 

#volcano plot where p. value cut off is 0.1
volcano_plot = function(top_table, title){
  v = EnhancedVolcano(top_table, x="log2FoldChange", y= "padj", lab=row.names(top_table), title = title,
                      pCutoff = 0.05, FCcutoff = 1.333, 
                      ylim = c(0, max(-log10(top_table$padj), na.rm=TRUE)),
                      col = c("grey30", "#A4D49C", "#4E80A5", "#F68D91"),
                      colAlpha = 1.0)
  #pdf(paste("../figures/", title, "volcano.pdf"), width = 8, height = 8)
  pdf(paste("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/figures/", title, "volcano.pdf"), width = 8, height = 8)
  print(v)
  dev.off()
  return(v)
}

#this function is for plotting a volcano plot using the top table generated from Sleuth Wald test
#the input top table has the columns: 
#"target_id"       "pval"            "qval"            "b"               "se_b"            "mean_obs"        "var_obs"        
# "tech_var"        "sigma_sq"        "smooth_sigma_sq" "final_sigma_sq" 
volcano_plot_sleuth = function(top_table, title, fold_change_cutoff, pvalue_cutoff){
  v = EnhancedVolcano(#top_table, x="b", y= "qval",  lab=top_table$target_id, title = title,
                      #top_table, x="b", y= "qval",  lab=lab_italics, title = title,
                      
                      top_table, x="avg_log2fc", y= "qval",  lab=lab_italics, title = title,
                      
                      pCutoff = pvalue_cutoff, FCcutoff = fold_change_cutoff, 
                     # ylim = c(0, max(-log10(top_table$qval), na.rm=TRUE)),
                      #col = c("grey30", "#A4D49C", "#4E80A5", "#F68D91"),
                     ylab= '-log10(q value)',
                     #selectLab = selectLabel,
                     selectLab = selectLab_italics,
                     boxedLabels = FALSE,
                     labSize = 5.0,           #size of label text
                     pointSize = 2.0,
                     drawConnectors = TRUE,
                     widthConnectors = 2.0,
                     colConnectors = 'black',
                     directionConnectors = 'both',
                     arrowheads = TRUE,
                     parseLabels = TRUE,
                     colAlpha = 0.30)  #this adjusts the transparency of the points 
  #pdf(paste("../figures/", title, "volcano.pdf"), width = 8, height = 8)
  pdf(paste("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/figures/", title, "likelihood_ratio_test_norm_fun_counts.pdf"), width = 8, height = 8)
  print(v)
  dev.off()
  return(v)
}


#this function is for plotting a volcano plot using the top table generated from Sleuth Wald test
#the input top table has the columns: 
#"target_id"       "pval"            "qval"            "b"               "se_b"            "mean_obs"        "var_obs"        
# "tech_var"        "sigma_sq"        "smooth_sigma_sq" "final_sigma_sq" 
volcano_plot_sleuth_wald_test = function(top_table, title, fold_change_cutoff, pvalue_cutoff){
  v = EnhancedVolcano(#top_table, x="b", y= "qval",  lab=top_table$target_id, title = title,
    #top_table, x="b", y= "qval",  lab=lab_italics, title = title,
    
    top_table, x="b", y= "qval", lab=top_table$target_id, title = title,
    
    pCutoff = pvalue_cutoff, FCcutoff = fold_change_cutoff, 
    # ylim = c(0, max(-log10(top_table$qval), na.rm=TRUE)),
    #col = c("grey30", "#A4D49C", "#4E80A5", "#F68D91"),
    ylab= '-log10(q value)',
    #selectLab = selectLabel,
    #selectLab = selectLab_italics,
    boxedLabels = FALSE,
    #labSize = 5.0,           #size of label text
    pointSize = 2.0,
    drawConnectors = TRUE,
    widthConnectors = 2.0,
    colConnectors = 'black',
    directionConnectors = 'both',
    arrowheads = TRUE,
    parseLabels = TRUE,
    colAlpha = 0.30)  #this adjusts the transparency of the points 
  #pdf(paste("../figures/", title, "volcano.pdf"), width = 8, height = 8)
  pdf(paste("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/figures/", title, "sleuth_wald_test_volcano.pdf"), width = 8, height = 8)
  print(v)
  dev.off()
  return(v)
}


#this function is for calculating the log2 fold change for pairs of monosome and polysome samples ie.WT-polysome-1 vs WT-monosome-1 fold change
#the input log2_df is a dataframe which is extracted with Sleuth package kallisto_table() where the samples have been normalized, containing
#columns for gene name (target_id) and the other columns are sample names which contain the log 2 transformed normalized scaled_counts_per_base.
#condition_A, condition_B indicate the 2 condition samples you want to compare
calculate_fold_change<-function(log2_df, condition_A, condition_B){
  log2_df$new <- (log2_df[condition_A] - log2_df[condition_B])
  
  #rename last column with conditions in fold change calculation
  colnames(log2_df)[ncol(log2_df)]<-paste(condition_A, "_vs_", condition_B, "_log2FC", sep="")
  log2_df <<-log2_df
  return(log2_df)
}


# Helper functions for main_analysis
#remove duplicated genes from table:
get_gene_names = function(data, gene_colname = "Gene names", protein_colname = "target_id"){
  # go through the gene names and pick the first one if there are multiple
  # if there are no gene names, take the protein name
  # expects data to be perseus output with columns 'Gene names' and 'Protein names'
  
  multi_names = strsplit(data[[gene_colname]], ";")
  protein_multi_names = strsplit(data[[protein_colname]], ";")
  gene_names = rep(NA, nrow(data))
  for(i in 1:length(multi_names)){
    gene_names[i] = multi_names[[i]][1]
    if(is.na(gene_names[i])){
      gene_names[i] = protein_multi_names[[i]][1]
    }
  }
  return(make.names(gene_names, unique = T))
}


# gets the top differentially expressed proteins using the limma package
get_top_table = function(data, design){
  fit = lmFit(data, design) 
  fit = eBayes(fit) 
  top_table =  topTable(fit, coef="difference", adjust="fdr", number = nrow(data))
  return(top_table)
}

# helper for getting differentially expressed proteins at specific timepoint
get_tp_top_table = function(data, groups, cell_type, timepoint){
  cell_type_data = data[, groups$cell_type == cell_type]
  tps = groups$tp[groups$cell_type == cell_type]
  
  difference = rep(0, ncol(cell_type_data))
  difference[tps == timepoint] = 1
  design = cbind(control=1, difference=difference)
  
  top_table =  get_top_table(cell_type_data, design)
  return(top_table)
}

volcano_plot = function(top_table, title){
  v = EnhancedVolcano(top_table, x="logFC", y= "adj.P.Val", lab=row.names(top_table), title = title,
                      pCutoff = 0.1, FCcutoff = 1.333, 
                      ylim = c(0, max(-log10(top_table$adj.P.Val), na.rm=TRUE)),
                      col = c("grey30", "#A4D49C", "#4E80A5", "#F68D91"),
                      colAlpha = 1.0)
  pdf(paste("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/figures/", title, "volcano.pdf"), width = 8, height = 8)
  print(v)
  dev.off()
  return(v)
}

get_pscore_table = function(top_table){
  direction_mult = sign(top_table$logFC)
  signed_pscore = -log10(top_table$adj.P.Val) * direction_mult
  pscore_table = cbind.data.frame(row.names(top_table), signed_pscore)
  pscore_table = pscore_table[order(pscore_table$signed_pscore, decreasing = T),]
  return(pscore_table)
}

get_boxplot = function(data, groups, title, adj_mult = 1,
                       order_groups = c(), comparisons = list(),
                       show_pair_pval = F, font_size = 15, line_thickness = 1.0,
                       x_lab = NULL, y_lab = "Expression", jitter = T){
  
  df = cbind.data.frame(data = as.numeric(data), groups = groups)
  
  if(length(order_groups) > 0){
    df$groups = factor(df$groups, levels = order_groups)
  }
  
  if(jitter){
    add = "jitter"
  }else{
    add = NULL
  }
  
  g =  ggboxplot(df, x = "groups", y = "data", color = "groups", 
                 add = add, legend = "none", size = line_thickness,
                 palette = c("#9DB4DD", "#F69679", "#F7AFC3",  "#A4D49C", "#FDC98F",
                             "#AB92C5", "#4E80A5", "#BD7691", "#FFF79F")) +
    rotate_x_text(angle = 45) 
  
  if(!is.null(x_lab)){
    g = g + xlab(x_lab)
  }
  if(!is.null(y_lab)){
    g = g + ylab(y_lab)
  }
  
  if(length(unique(df$groups)) == 2){
    df$groups = as.character(df$groups)
    groups1 = df$data[df$groups == unique(df$groups)[1]]
    groups2 = df$data[df$groups == unique(df$groups)[2]]
    test = t.test(groups1, groups2)
    
    if(adj_mult == 1){
      tab = cbind.data.frame(test$p.value)
      colnames(tab) = c("p.value")
      p_tab <- tableGrob(tab, rows = NULL)
    }else{
      tab = rbind.data.frame(test$p.value, (test$p.value)*adj_mult)
      rownames(tab) = c("p", "adj.p")
      p_tab <- tableGrob(tab, cols = NULL)
    }
    
    g = g + geom_signif(comparisons = list(unique(df$groups)), 
                        test = "t.test",
                        map_signif_level=c("***"=(0.001/adj_mult), "**"=(0.01/adj_mult), "*"=(0.05/adj_mult)), 
                        step_increase = c(0.1, 0.1, 0.1),
    )
    
    if(show_pair_pval){
      g = ggdraw() + draw_plot(g, width = 0.7) + draw_plot(p_tab, x = 0.7, width = 0.3)
    }
  }else if(length(comparisons > 0)){
    g = g + geom_signif(comparisons = comparisons, 
                        test = "t.test",
                        map_signif_level=c("***"=(0.001/adj_mult), "**"=(0.01/adj_mult), "*"=(0.05/adj_mult)), 
                        step_increase = c(0.1, 0.1, 0.1), +
                          geom_hline(yintercept = mean(df$data), linetype = 2) # Add horizontal line at base mean
                        
    )
  }else{
    if(show_pair_pval){
      label_type = "p.format"
    }else{
      label_type = "p.signif"
    }
    
    g = g + stat_compare_means(method = "anova", label.y = (max(df$data) + sd(df$data)),
                               label.x.npc = 0.1)+ # Add global annova p-value
      stat_compare_means(label = label_type, method = "t.test",  # alternatively: label = "p.format"
                         ref.group = ".all.", hide.ns = T, label.y.npc = "top")+# Pairwise comparison against all
      geom_hline(yintercept = mean(df$data), linetype = 2) # Add horizontal line at base mean
  }
  
  # apply font size
  g = g + theme(text = element_text(size = font_size))
  g = g + ggtitle(title)
  
  return(g)
}

venn_diagram = function(list_of_sets, names, title, output=T){
  myCol <- brewer.pal(max(3, length(names)), "Pastel2")[1:length(names)]
  
  v = venn.diagram(
    x = list_of_sets,
    category.names = names,
    filename = paste("figures/", title, ".tiff"),
    output=output,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans"
  )
  return(v)
}

select_columns = function(data, colname1, colname2){
  first_col = which(colnames(data) == colname1)
  last_col = which(colnames(data) == colname2)
  return(data[,first_col : last_col])
}

### re-write this as a function to input different combinations of principle components
#the input prcomp_df is the dataframe containing the PCA plot values, which is the output of the prcomp function
#first_pc is the first principle component you want to look at (numeric), and second_pc is the second principle
#component you want to look at 
pca_plot_combo <- function(prcomp_df, first_pc, second_pc, first_pc_string, second_pc_string){
  
  celf2_data_pc = data.frame(x=prcomp_df[,first_pc], y=prcomp_df[,second_pc], condition=metadata$condition, sample=row.names(prcomp_df))                   
  
  #setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/figures")
  
  #png(file=paste0("Celf2_PCA_", first_pc, "_", second_pc, ".png"))
  
  p <- (ggplot(celf2_data_pc, aes(x=x, y=y, color=condition, label = sample)) +
          geom_point(size=6) +
          scale_color_manual("Condition",
                             values = c("#A4D49C", "#4E80A5", "darkblue", "red", "gold", "purple")) +
          theme_bw() +
          # alpha = 0.1 +
          xlab(paste0("PC", first_pc_string)) +
          ylab(paste0("PC", second_pc_string)) +
          theme(text = element_text(size = 14)))
  
  
  
  p + geom_text(check_overlap = TRUE, colour = "black", vjust = "inward", hjust = "inward", nudge_x = 0.05) + theme(aspect.ratio=1) 
  
  #dev.off()
  
}

#function for plotting scatterplot of CELF2 WTpoly-vs_WTtotal(x-axis) and CELF2 WTpoly-vs_WTmono(y-axis);
#overlay colour of dots with CSDE1 vs IgG 
#input is the merged sleuth table wald test results of WTpoly-vs_WTtotal, CELF2 WTpoly-vs_WTmono, and CSDE1 vs IgG
celf2_csde1_scatterplot<- function(merged_table){
    p <- ggplot(merged_table %>%
                arrange(desc(enrichment)),                 ## this alters the order in which the points are plotted, so that CSDE1 plots are on top
              aes(x=b.x, y=b.y, col=enrichment, alpha=0.1)) +
    #scale_color_manual(values = c( "red", "darkblue","grey"))+
    scale_color_manual(values = c( "coral1", "lightgrey"))+

    geom_point(size = 5) +
    labs(x=("wald test b value WT poly vs WT total"),
         y=("wald test b value WT poly vs WT mono"),
         
         title=("CELF2-KI_Polysome-seq - not filtered for protein coding genes")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  plot(p)
}

### this version of scatterplot only colours the enriched CSDE1 gene targets in a gradient
celf2_csde1_scatterplot2<- function(merged_table){
   p <- geom_point(data = subset(merged_table, b > 2), 
                aes(x=b.x, y=b.y, col=enrichment, alpha=0.1, fill=b, size = 5)) +
      scale_fill_gradient2(low="cyan4", mid = "white", high = "coral1")+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))
   
   plot(p)
}

#this first plots all genes in light grey in a scatterplot x-axis WTpoly-WTtotal, and y-axis WTpoly-WTmono
celf2_csde1_scatterplot3<- function(merged_table){
  p <- ggplot(data = merged_table,               ## this alters the order in which the points are plotted, so that CSDE1 plots are on top
              aes(x=b.x, y=b.y, alpha=0.1)) +
    
    geom_point(stat="identity", size = 5, colour="lightgrey", show.legend = FALSE) +
    labs(x=("wald test b value WT poly vs WT total"),
         y=("wald test b value WT poly vs WT mono"),
         
         title=("CELF2-KI_Polysome-seq - filtered for protein coding genes")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  plot(p)
}

