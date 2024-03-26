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
                      
                      #top_table, x="avg_log2fc", y= "qval",  lab=lab_italics, title = title,
                      top_table, x="b", y= "qval",  lab=lab_italics, title = title,
    
                      pCutoff = pvalue_cutoff, FCcutoff = fold_change_cutoff, 
                     # ylim = c(0, max(-log10(top_table$qval), na.rm=TRUE)),
                      #col = c("grey30", "#A4D49C", "#4E80A5", "#F68D91"),
                     ylab= '-log10(q value)',
                     col = c("lightgrey", "cyan3", "#A4D49C",  "coral2"),
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
                     colAlpha = 0.50)  #this adjusts the transparency of the points 
  #pdf(paste("../figures/", title, "volcano.pdf"), width = 8, height = 8)
  #pdf(paste("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/figures/", title, "wald_test_volcano.pdf"), width = 8, height = 8)
  #pdf(paste("C:/Users/queenie.tsang/Desktop/CSDE1/ribodetector/figures", title, "wald_test_volcano.pdf"), width = 8, height = 8)
  
  #print(v)
  plot(v)
  #dev.off()
  #return(v)
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
    xlab='b value',
    #selectLab = selectLabel,
    #selectLab = selectLab_italics,
    boxedLabels = FALSE,
    #labSize = 5.0,           #size of label text
    
    #legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", expression(p - value ~ and ~ log[2] ~ FC)),
    legendLabels = c("NS", expression(Log[] ~ FC), "p-value", expression(p - value ~ and ~ log[] ~ FC)),
    
    pointSize = 2.0,
    drawConnectors = TRUE,
    widthConnectors = 2.0,
    colConnectors = 'black',
    directionConnectors = 'both',
    arrowheads = TRUE,
    parseLabels = TRUE,
    colAlpha = 0.30)  #this adjusts the transparency of the points 
  #pdf(paste("../figures/", title, "volcano.pdf"), width = 8, height = 8)
  #pdf(paste("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/figures/", title, "sleuth_wald_test_volcano.pdf"), width = 8, height = 8)
  
  #pdf(paste("C:/Users/queenie.tsang/Desktop/CSDE1/ribodetector/figures", title, "wald_test_volcano.pdf"), width = 8, height = 8)
  
  print(v)
  dev.off()
  return(v)
}

volcano_plot_sleuth_wald_test_xaxis_log2 = function(top_table, title, fold_change_cutoff, pvalue_cutoff){
  v = EnhancedVolcano(#top_table, x="b", y= "qval",  lab=top_table$target_id, title = title,
    #top_table, x="b", y= "qval",  lab=lab_italics, title = title,
    
    top_table, x="b", y= "qval", lab=top_table$target_id, title = title,
    
    pCutoff = pvalue_cutoff, FCcutoff = fold_change_cutoff, 
    # ylim = c(0, max(-log10(top_table$qval), na.rm=TRUE)),
    #col = c("grey30", "#A4D49C", "#4E80A5", "#F68D91"),
    ylab= '-log10(q value)',
    xlab='log2 b value',
    #selectLab = selectLabel,
    #selectLab = selectLab_italics,
    boxedLabels = FALSE,
    #labSize = 5.0,           #size of label text
    
    #legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", expression(p - value ~ and ~ log[2] ~ FC)),
    legendLabels = c("NS", expression(Log[2] ~ FC), "q-value", expression(q - value ~ and ~ log[2] ~ FC)),
    
    pointSize = 2.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    directionConnectors = 'both',
    arrowheads = TRUE,
    parseLabels = TRUE,
    colAlpha = 0.30)  #this adjusts the transparency of the points 
  #pdf(paste("../figures/", title, "volcano.pdf"), width = 8, height = 8)
  #pdf(paste("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/figures/", title, "sleuth_wald_test_volcano.pdf"), width = 8, height = 8)
  
  #pdf(paste("C:/Users/queenie.tsang/Desktop/CSDE1/ribodetector/figures", title, "wald_test_volcano.pdf"), width = 8, height = 8)
  
  print(v)
  dev.off()
  return(v)
}

volcano_plot_sleuth_wald_test_xaxis_log2_labels = function(top_table, title, fold_change_cutoff, pvalue_cutoff){
  v = EnhancedVolcano(#top_table, x="b", y= "qval",  lab=top_table$target_id, title = title,
    #top_table, x="b", y= "qval",  lab=lab_italics, title = title,
    
    top_table, x="b", y= "qval", lab=top_table$target_id, title = title,
    
    pCutoff = pvalue_cutoff, 
    
    FCcutoff = fold_change_cutoff, 
    # ylim = c(0, max(-log10(top_table$qval), na.rm=TRUE)),
    #col = c("grey30", "#A4D49C", "#4E80A5", "#F68D91"),
    ylab= '-log10(q value)',
    xlab='log2 b value',
    #selectLab = selectLabel,
    #selectLab = selectLab_italics,
    boxedLabels = FALSE,
    #labSize = 5.0,           #size of label text
    
    #legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", expression(p - value ~ and ~ log[2] ~ FC)),
    legendLabels = c("NS", expression(Log[2] ~ FC), "q-value", expression(q - value ~ and ~ log[2] ~ FC)),

    pointSize = 2.0,
   # drawConnectors = TRUE,
    widthConnectors = 2.0,
    colConnectors = 'black',
    #directionConnectors = 'both',
    #arrowheads = TRUE,
    parseLabels = TRUE,
    colAlpha = 0.30)  #this adjusts the transparency of the points 
  #pdf(paste("../figures/", title, "volcano.pdf"), width = 8, height = 8)
  #pdf(paste("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/figures/", title, "sleuth_wald_test_volcano.pdf"), width = 8, height = 8)
  
  #pdf(paste("C:/Users/queenie.tsang/Desktop/CSDE1/ribodetector/figures", title, "wald_test_volcano.pdf"), width = 8, height = 8)
  
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

volcano_plot = function(top_table, title, pvalue_cutoff, fc_cutoff){
  v = EnhancedVolcano(top_table, x="logFC", y= "adj.P.Val", lab=row.names(top_table), title = title,
                      pCutoff = pvalue_cutoff, FCcutoff = fc_cutoff, 
                      ylim = c(0, max(-log10(top_table$adj.P.Val), na.rm=TRUE)),
                      col = c("grey30", "#A4D49C", "#4E80A5", "#F68D91"),
                      colAlpha = 1.0)
  pdf(paste("C:/Users/queenie.tsang/Desktop/CELF2-KI-Polysomeseq-P5/Celf2-KI-Polysome-seq/figures/2024/", title, "volcano.pdf"), width = 8, height = 8)
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
  celf2_data_pc <<- celf2_data_pc
  #setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/figures")
  
  #png(file=paste0("Celf2_PCA_", first_pc, "_", second_pc, ".png"))
  
  p <- (ggplot(celf2_data_pc, aes(x=x, y=y, color=condition, label = sample)) +
          geom_point(size=6) +
          scale_color_manual("Condition",
                             values = c("#A4D49C", "#4E80A5", "darkblue", "red", "gold", "purple", "darkgreen", "pink")) +
          theme_bw() +
          # alpha = 0.1 +
          xlab(paste0("PC", first_pc_string)) +
          ylab(paste0("PC", second_pc_string)) +
          theme(text = element_text(size = 14)))
  
  
  
  p + geom_text(check_overlap = TRUE, colour = "black", vjust = "inward", hjust = "inward", nudge_x = 0.05) + theme(aspect.ratio=1) 
  
  #dev.off()
  
}

pca_plot_combo2 <- function(prcomp_df, s2c, first_pc, second_pc, first_pc_string, second_pc_string){
  
  celf2_data_pc = data.frame(x=prcomp_df[,first_pc], y=prcomp_df[,second_pc], sample=row.names(prcomp_df))    
  s2c <- s2c[order(s2c$condition),]
  celf2_data_pc$condition <- s2c$condition
  celf2_data_pc <<- celf2_data_pc
  #setwd("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/figures")
  
  #png(file=paste0("Celf2_PCA_", first_pc, "_", second_pc, ".png"))
  
  p <- (ggplot(celf2_data_pc, aes(x=x, y=y, color=condition, label = sample)) +
          geom_point(size=6) +
          scale_color_manual("Condition",
                             values = c("#A4D49C", "#4E80A5", "darkblue", "red", "gold", "purple", "darkgreen", "pink")) +
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
              aes(x=b.x, y=b.y, alpha=0.1, label = target_id)) +
    
    #label select genes that are layer marker genes
    geom_label(data=subset(merged_table, target_id %in% c("TBR1", "SATB2", "TLE4", "POU3F3", "FEZF2", "FEZF1")))+
  
    
    #geom_point(stat="identity", size = 5, colour="lightgrey", show.legend = FALSE) +
    geom_point(stat="identity", size = 2, colour="cyan4", show.legend = FALSE) +
    
    labs(x=("wald test b value WT poly vs WT total"),
         y=("wald test b value WT poly vs WT mono"),
         
         title=("CELF2-KI_Polysome-seq - filtered for protein coding genes")) #+
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #      panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  plot(p)
}

#this function is for retrieving the list of genes within a specified box on the scatterplot
#of CELF2 WTpoly-WTtotal vs WTpoly-vs-WTmono and overlay of CSDE1 vs IgG b value
#input min.x is the minimum x value; max.x is maximum x value
#input min.y is the minimum y value; and max.y is the maximum y value
get_box_genes <- function(merged_df, max.x, min.x, max.y, min.y ){
  box_x <- merged_df[merged_df$b.x <= max.x & merged_df$b.x >= min.x,]
  box_y <- box_x[box_x$b.y <= max.y & box_x$b.y >= min.y,]
  return(box_y)
}

volcano_plot2 = function(top_table, title, pvalue_cutoff, log2fc_cutoff){
  v = EnhancedVolcano(top_table, x="b", y= "qval", lab=(top_table$target_id), title = title,
                      pCutoff =pvalue_cutoff, FCcutoff = log2fc_cutoff, 
                      ylab= '-log10(q value)',
                      ylim = c(0, max(-log10(top_table$qval), na.rm=TRUE)),
                      col = c("lightgrey", "cyan3", "lightgrey",  "lightcoral"),
                      colAlpha = 1.0)
  #pdf(paste("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/figures/", title, "volcano.pdf"), width = 8, height = 8)
 # pdf(paste("C:/Users/queenie.tsang/Desktop/CSDE1/ribodetector/figures", title, "volcano.pdf"), width = 8, height = 8)
  print(v)
  dev.off()
  return(v)
}

volcano_plot3 = function(top_table, title, pvalue_cutoff, log2fc_cutoff){
  v = EnhancedVolcano(top_table, x="logFC", y= "adj.P.Val", lab=(top_table$target_id), title = title,
                      pCutoff =pvalue_cutoff, FCcutoff = log2fc_cutoff, 
                      ylab= '-log10(q value)',
                      ylim = c(0, max(-log10(top_table$`adj.P.Val`), na.rm=TRUE)),
                      col = c("lightgrey", "cyan3", "lightgrey",  "lightcoral"),
                      colAlpha = 1.0)
  #pdf(paste("C:/Users/queenie.tsang/Desktop/CELF2/Celf2-KI-Polysome-seq/figures/", title, "volcano.pdf"), width = 8, height = 8)
  # pdf(paste("C:/Users/queenie.tsang/Desktop/CSDE1/ribodetector/figures", title, "volcano.pdf"), width = 8, height = 8)
  print(v)
  dev.off()
  return(v)
}

volcano_plot4 = function(top_table, title, pvalue_cutoff, fc_cutoff){
  v = EnhancedVolcano(top_table, x="logFC", y= "adj.P.Val", lab=row.names(top_table), title = title,
                      pCutoff = pvalue_cutoff, FCcutoff = fc_cutoff, 
                      ylab= '-log10(q value)',
                      legendLabels=c('Not sig.','Log (base 2) FC','q-value',
                                     'q-value & Log (base 2) FC'),
                      ylim = c(0, max(-log10(top_table$adj.P.Val), na.rm=TRUE)),
                      col = c("grey30", "#A4D49C", "#4E80A5", "#F68D91"),
                      colAlpha = 1.0)
  pdf(paste("C:/Users/queenie.tsang/Desktop/CELF2-KI-Polysomeseq-P5/Celf2-KI-Polysome-seq/figures/2024/", title, "volcano.pdf"), width = 8, height = 8)
  print(v)
  dev.off()
  return(v)
}

#this function is for getting the sleuth object for 2 sample condition
#the input metadata_table is the subset of the metadata table containing only the rows that are the for the 
#conditions being tested ie. wild type poly (wt_poly) wild type total(wt_total)
#the input "reference" is the condition that you want to set as the reference (control condition)
set_metadata_control_condition <- function(metadata_table, reference){
  metadata_table$condition <- as.factor(metadata_table$condition)
  metadata_table$condition <- relevel(metadata_table$condition, ref= reference)
  return(metadata_table)
}

#this function is to run sleuth_prep (which initializes the sleuth object,and then run sleuth_fit which will fit the full model)
#the input metadata_table is the same as the output of set_metadata_control_condition function which is a metadata df containing 
# only samples for the 2 conditions tested, and that the reference/control condition as been set as the reference, as a factor type
get_sleuth_object <- function(metadata_table){
  so <- sleuth_prep(metadata_table,  extra_bootstrap_summary = TRUE, aggregation_column="gene", gene_mode=TRUE, extra_bootstrap_summary = TRUE, target_mapping=t2g, transform_fun_counts=function (x) log2(x+0.5))
  so <- sleuth_fit(so, ~condition, 'full')
  return(so)
}


#this function is for running the Sleuth wald test and the subsequent processing steps for cleaning up the dataframe
#so that it contains only protein coding genes
#input "so" is the sleuth_object output by the get_sleuth_object function
#input "condition_tested" is the test condition 
sleuth_wald_test_processing <- function(so, condition_tested ){
  wald_test_df <- sleuth_wt(so, which_beta = condition_tested)
  sleuth_wald_test <- sleuth_results(wald_test_df, test = condition_tested, show_all = TRUE)
  
  #remove the NA genes
  sleuth_wald_test_no_NAs <- na.omit(sleuth_wald_test)
  
  #filter for protein coding genes
  sleuth_wald_test_no_NAs$target_id <- toupper(sleuth_wald_test_no_NAs$target_id)
  sleuth_wald_test_protein_coding <- sleuth_wald_test_no_NAs[sleuth_wald_test_no_NAs$target_id %in% protein_coding$gene_caps,] 
  
  return(sleuth_wald_test_protein_coding)
}

# function to plot the NES scores for DCX GSEA and SOX2 GSEA pathways on separate plots
#cell_type is a string specifying cell type ie. "DCX" or "SOX2"
#filtered_df is the filtered gsea results dataframe (filtered by NOM pvalue and FDR q value)
plot_NES <-function(filtered_df, plottitle){
  p <- ggplot(filtered_df, aes(x = condition, y = Pathways, size = NES, color = FDR.q.val )) +
    geom_point(alpha = 0.7) + #scatterplot
    scale_size(range = c(2, 5), name = "NES") +  #change the size of the points and create a name for the size legend
    scale_color_viridis(discrete = FALSE)
  plot(p) +
    ggtitle = plottitle
  #function for saving figures as pdf, png or jpg
  #set directory for saving GSEA results comparing organoid secretome timepoints
  #setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/figures/timepoint comparisons GSEA")
  
  #set directory for GSEA results for comparing organoid DCX vs SOX2 pathways
  setwd("C:/Users/queenie.tsang/Desktop/CELF2-KI-Polysomeseq-P5/Celf2-KI-Polysome-seq/figures/2024/")
  
  
  #save as png file, to save as pdf write "pdf", or jpg, write "jpg"; res = resolution
  pdf(filename=paste0(plottitle, "_gsea_plot.pdf"), width=3000, height=3000, res=300)
  
  #function that makes the plot
  plot(p)
  dev.off()
}  

unique.comparisons <- function(dataframe1,dataframe2){
  #'unique.comparisons
  #'
  #'This function takes two vectors as inputs and returns the elements in
  #'vector1 that are unique to both cases
  #'
  #'@param dataframe1 this is the first vector to be compared to the second vector1
  #'@param dataframe2 this is the second vector to be compared
  #splice out the names of the processes or whatever is being compared.
  #this was originally designed to compare targets in GSEA
  vector1 <- dataframe1[,'Pathways']
  vector2 <- dataframe2[,'Pathways']
  #generates a vector of booleans that show elements that are unique to vector1
  #when compared to each of the other vectors
  # boolean.comparison1 <- (!(vector1 %in% vector2) & !(vector1 %in% vector3)
  #                         & !(vector1 %in% vector4) & !(vector1 %in% vector5))
  
  boolean.comparison1 <- (!(vector1 %in% vector2))
  #gets only the rows that are in unique from the original dataframe of vector1
  unique.dataframe <- dataframe1[boolean.comparison1,]
  #returns the dataframe for future use
  return(unique.dataframe)
}


# ##### #function for filtering the GSEA results to extract only one particular Gene Ontology database ie. GOCC only
# #input GO_db is the name of the Gene Ontology database ie. "GOCC", "GOBP", "GOMF"
#input test_condition is the GSEA tsv results for the test condition (ie. mutant samples)
#input ref_condition is the GSEA tsv results for the reference condition (ie. wild type)
extract_specific_GO_db <-function(test_condition, ref_condition, GO_db, sampletype1, sampletype2){
  
  #add a column for sample type:
  test_condition$condition <- sampletype1
  ref_condition$condition <-sampletype2
  
  #rename "NAME" column to "Pathways"
  colnames(test_condition)[1]<- "Pathways"
  colnames(ref_condition)[1]<- "Pathways"
  
  test_GO <- filter(test_condition, grepl(GO_db, test_condition$Pathways))
  ref_GO <- filter(ref_condition, grepl(GO_db, ref_condition$Pathways))

  #extract the GO databases to the global environment to access 
  test_GO <<- test_GO
  ref_GO <<- ref_GO
}


#function to filter GSEA pathways by NOMp.value
#the input is a dataframe of pathways from GSEA for a particular condition ie. mutant_mono_GOBP
#input parameter "nompvalue" is the nominal p value cutoff by which you want to filter the GSEA results
NOM_pvalue_filter<- function(condition_GO_df, nom_pvalue){
  
  condition_pvalue <- condition_GO_df[condition_GO_df$`NOM.p.val` < nom_pvalue, ]
  return(condition_pvalue)
}
