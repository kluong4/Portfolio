# =============================================================================
# Script 02: Differential Expression Analysis
# Project:   MBD6 Haploinsufficiency and ASD - Transcriptomic Analysis
# Course:    BMI 5332: Statistical Analysis of Genomic Data
#            UT Health Houston, Spring 2026
# Input:     data/processed/norm_filtered.rds  (from 01_preprocessing.R)
# Output:    data/processed/DEG_results.rds
#            results/tables/Table1_fullDEGresults.csv
#            results/tables/Table2_upregulated_DEGs.csv
#            results/tables/Table3_downregulated_DEGs.csv
#            results/tables/Table4_top10_DEGs.csv
#            results/figures/Fig_VolcanoPlot.png
#            results/figures/Fig_DEG_Heatmap.png
# =============================================================================

# ----- Load Libraries --------------------------------------------------------
library(limma)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)

# ----- Set Working Directory -------------------------------------------------
setwd("~/Desktop/BMI 5332_R/Project")

# ----- Load Preprocessed Data ------------------------------------------------
norm_filtered <- readRDS("Data/Processed/norm_filtered.rds")


# =============================================================================
# STEP 1: Fit Linear Model and Extract DEGs
# =============================================================================
# Design matrix: 2 CTL, 2 MBD6het
design <- model.matrix(~0 + factor(c("CTL", "CTL", "MBD6het", "MBD6het")))
colnames(design) <- c("CTL", "MBD6het")

# Fit linear model
fit <- lmFit(norm_filtered, design)
contrast_matrix <- makeContrasts(MBD6het - CTL, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract all results with BH correction
results <- topTable(fit2, coef = 1, number = Inf, adjust = "BH")

# Keep least p-value for duplicated probes
results_sorted <- results[order(results$adj.P.Val), ]
results_unique <- results_sorted[!duplicated(results_sorted$ProbeName), ]

# Categorize DEGs
results_unique$category <- "Not Significant"
results_unique$category[results_unique$adj.P.Val < 0.05 &
                          results_unique$logFC >  1] <- "Upregulated"
results_unique$category[results_unique$adj.P.Val < 0.05 &
                          results_unique$logFC < -1] <- "Downregulated"
results_unique$category <- factor(results_unique$category,
                                  levels = c("Upregulated",
                                             "Downregulated",
                                             "Not Significant"))

# Subset significant DEGs
sig_DEGs <- results_unique[results_unique$adj.P.Val < 0.05 &
                                  abs(results_unique$logFC) > 1, ]   # 913
upregulated <- sig_DEGs[sig_DEGs$logFC >  1, ]                     # 637
downregulated <- sig_DEGs[sig_DEGs$logFC < -1, ]                     # 276

cat("Total significant DEGs: ", nrow(sig_DEGs), "\n")
cat("Upregulated: ", nrow(upregulated), "\n")
cat("Downregulated: ", nrow(downregulated), "\n")


# =============================================================================
# STEP 2: Volcano Plot
# =============================================================================
sig_DEGs_sorted <- sig_DEGs[order(sig_DEGs$adj.P.Val, -abs(sig_DEGs$logFC)), ]
top10 <- sig_DEGs_sorted[1:10, ]

results_for_volcPlot <- results_unique
results_for_volcPlot$sig <- results_for_volcPlot$adj.P.Val < 0.05 &
  abs(results_for_volcPlot$logFC) > 1
results_for_volcPlot$label <- ifelse(
  results_for_volcPlot$ProbeName %in% top10$ProbeName,
  results_for_volcPlot$GeneSymbol, "")

png("Results/Fig5_VolcanoPlot.png", width = 12, height = 8, units = "in", res = 300)

ggplot(results_for_volcPlot,
       aes(x = logFC, y = -log10(adj.P.Val), color = AveExpr, shape = sig)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_viridis_c(trans = "log10", name = "Ave Expression") +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    name = "Sig. DE",
    labels = c("FALSE" = "FALSE", "TRUE" = "TRUE")) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "blue",  linewidth = 0.6) +
  geom_vline(xintercept =  1, linetype = "dashed", color = "red",   linewidth = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text_repel(aes(label = label),
                  size = 3.5, max.overlaps = Inf,
                  box.padding = 1.0, point.padding = 0.5, force = 15,
                  segment.color = "grey50", segment.size = 0.3,
                  color = "black", fontface = "bold",
                  min.segment.length = 0) +
  labs(title = paste0("MBD6het vs Control, Differential Gene Expression\n",
                      "(", nrow(results_for_volcPlot), " probes analyzed)"),
       x = "log2 (Fold Change)",
       y = "-log10(FDR-Adjusted p-values)") +
  theme_bw() +
  theme(
    plot.title   = element_text(size = 12, face = "bold"),
    axis.title   = element_text(size = 11),
    axis.text    = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9))

dev.off()

# =============================================================================
# STEP 3: DEG Expression Heatmap
# =============================================================================
sig_probes <- norm_filtered$genes$ProbeName %in% sig_DEGs$ProbeName
sig_expr <- norm_filtered$E[sig_probes, ]
sig_genes_info <- norm_filtered$genes[sig_probes, ]

# Match exactly using composite key
sig_genes_info$key <- paste(sig_genes_info$ProbeName,
                            sig_genes_info$Row,
                            sig_genes_info$Col, sep = "_")
sig_DEGs$key <- paste(sig_DEGs$ProbeName,
                            sig_DEGs$Row,
                            sig_DEGs$Col, sep = "_")
exact_match <- sig_genes_info$key %in% sig_DEGs$key
sig_expr <- sig_expr[exact_match, ]
sig_genes_info <- sig_genes_info[exact_match, ]

rownames(sig_expr) <- sig_genes_info$ProbeName
colnames(sig_expr) <- c("CTL_1", "CTL_2", "MBD6het_1", "MBD6het_2")

# Z-score scaling per row
sig_expr_scaled <- t(scale(t(sig_expr)))

top_annotation <- HeatmapAnnotation(
  Samples = c("MBD6het", "MBD6het", "Control", "Control"),
  col     = list(Samples = c(Control = "#A465BA", MBD6het = "#1B9E77")),
  annotation_legend_param = list(
    Samples = list(
      title_gp  = gpar(fontsize = 11, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)))
)

png("Results/Fig6_DEG_Heatmap.png", width = 7, height = 5, units = "in", res = 300)

Heatmap(sig_expr_scaled,
        name             = "z-score",
        col              = colorRamp2(c(-2, 0, 2), c("red", "white", "blue")),
        top_annotation   = top_annotation,
        cluster_rows     = TRUE, cluster_columns = TRUE,
        show_row_names   = FALSE, show_column_names = TRUE,
        column_names_rot = 45,
        column_title     = "DEGs Expression Heatmap (FDR < 0.05, |logFC| > 1)",
        column_title_side = "top",
        column_title_gp  = gpar(fontsize = 13, fontface = "bold"))

dev.off()

# =============================================================================
# STEP 4: Save DEG Tables
# =============================================================================
keep_cols <- c("ProbeName", "GeneSymbol", "logFC", "AveExpr", "t", "adj.P.Val", "category")
full_DEG_clean <- results_unique[, keep_cols]
full_DEG_clean <- full_DEG_clean[order(full_DEG_clean$adj.P.Val), ]
write.csv(full_DEG_clean, "Results/Table1_fullDEGresults.csv", row.names = FALSE)

sig_cols <- c("ProbeName", "GeneSymbol", "logFC", "AveExpr", "t", "adj.P.Val")

upregulated_clean   <- full_DEG_clean[full_DEG_clean$category == "Upregulated", sig_cols]
write.csv(upregulated_clean, "Results/Table2_upregulated_DEGs.csv", row.names = FALSE)

downregulated_clean <- full_DEG_clean[full_DEG_clean$category == "Downregulated", sig_cols]
write.csv(downregulated_clean, "Results/Table3_downregulated_DEGs.csv", row.names = FALSE)

top10_clean <- top10[, keep_cols]
write.csv(top10_clean, "Results/Table4_top10_DEGs.csv", row.names = FALSE)

# =============================================================================
# SAVE: DEG objects for downstream scripts
# =============================================================================
saveRDS(list(
  results_unique = results_unique,
  sig_DEGs = sig_DEGs,
  upregulated = upregulated,
  downregulated = downregulated,
  top10 = top10
), file = "Data/Processed/DEG_results.rds")