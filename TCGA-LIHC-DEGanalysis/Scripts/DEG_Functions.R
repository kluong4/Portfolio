library(ggplot2)
library(ggrepel)

nicePCA <- function(pca_object, group_vector){
  ### Extract PCA data ###
  # Convert principal components matrix to data.frame
  pca_data <- as.data.frame(pca_object)
  # Add grouping information to pca_data
  pca_data$Group <- factor(group_vector)
  
  ### PCA plot ###
  pca_plot <- ggplot(pca_data, aes(x = x.PC1, y = x.PC2, color = Group, shape = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = "PCA Plot", color = "Group", shape = "Group",
         x = paste0("PC1 (", round(pca_data$pve[1]*100, 2), "% variance explained)"),
         y = paste0("PC2 (", round(pca_data$pve[2]*100, 2), "% variance explained)")) +
    scale_color_discrete(labels = c("1", "2")) +  # Change legend labels without changing colors
    scale_shape_manual(values = c(16, 17), labels = c("1", "2")) +  # Change legend labels without changing shapes
    coord_fixed() +  # Ensure 1:1 aspect ratio
    theme_bw()  # Use a black-and-white theme
  
  ### Return PCA Plot ###
  return(pca_plot)
}

niceVolcano <- function(df, padj_cutoff, logfc_cutoff = 1.5) {
  # Num of Loci for title
  loci_analyzed <- nrow(df)
  # Remove rows with missing p-values
  df <- df[!is.na(df$padj), ]
  
  # Cutoffs
  df$padj_fdrCutoff <- df$padj <= padj_cutoff 
  
  # Select top significant genes for labeling
  label_genes <- head(df[order(df$padj), ], 10)
  
  # Randomize order of the rows
  df <- df[sample(nrow(df), size = nrow(df)), ]
  
  # Create the volcano plot
  volcano_plot <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), 
                                 color = baseMean, shape = padj_fdrCutoff)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = -log10(padj_cutoff), 
               linetype = "dashed", color = "black") +
    geom_vline(xintercept = -logfc_cutoff, 
               linetype = "dashed", color = "blue") +
    geom_vline(xintercept = logfc_cutoff, 
               linetype = "dashed", color = "red") +
    geom_text_repel(data = label_genes, aes(label = gene_name), 
                    size = 3, max.overlaps = Inf, color = "black") +
    labs(title = paste0("TCGA-LIHC, Differential Gene Expression\n",
                        "Between AJCC Pathological Stage II vs\n",
                        "Stage III, IIIA, IIIB, IIIC\n(",
                        loci_analyzed, " loci analyzed)"), 
         x = "log2 (Fold Change)", 
         y = "-log10(FDR-Adjusted p-values)", 
         color = "Base Mean Count", 
         shape = "Sig. DE") +
    scale_color_viridis_c(trans = "log10") + 
    theme_bw()
  
  # Return the Volcano plot
  return(volcano_plot)
}

resSummary <- function(res_df) {
  ### list of the questions ###
  question <- list()
  question[1] <- paste0("How many genes were included in the DESeq2 analysis ", 
                        "after filtering?")
  question[2] <- paste0("How many genes included in the DESeq2 analysis had ", 
                        "non-missing FDR-adjusted p-values?")
  question[3] <- paste0("Note: using 0.05 as the p-adjusted cutoff for ", 
                        "significance for the following calculations:")
  question[4] <- paste0("How many genes were significantly down-regulated in ", 
                        "group2 relative to group1?")
  question[5] <- paste0("How many genes were significantly up-regulated in ", 
                        "group2 relative to group1?")
  question[6] <- paste0("What were the gene names of the top 10 most ", 
                        "differentially expressed genes ranked by ", 
                        "significance?")
  
  ### Printing the Results ###
  # Total genes included in the DESeq2 analysis
  message(question[1])
  totGenes <- nrow(res_df)
  print(totGenes)
  
  # Total genes with non-missing FDR-adjusted p-values
  message(question[2])
  non_missing_fdr <- sum(!is.na(res_df$padj))
  print(non_missing_fdr)
  
  # Identify significantly down-regulated genes (FDR < 0.05, log2FoldChange < 0)
  message(question[3])
  sig_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)
  message(question[4])
  print(sig_down)
  
  # Identify significantly up-regulated genes (FDR < 0.05, log2FoldChange > 0)
  sig_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
  message(question[5])
  print(sig_up)
  
  # Extract top 10 most significantly, differentially expressed genes 
  message(question[6])
  top_genes <- head(res_df[order(res_df$padj), ], 10)
  print(top_genes$gene_name)
}