# =============================================================================
# Script 04: eQTL-GWAS Integration Analysis
# Project:   MBD6 Haploinsufficiency and ASD - Transcriptomic Analysis
# Course:    BMI 5332: Statistical Analysis of Genomic Data
#            UT Health Houston, Spring 2026
# Input:     data/processed/DEG_results.rds  (from 02_differential_expression.R)
#            GTEx v11 Brain Cortex eGenes + significant eQTL pairs
#            GWAS Catalog ASD associations (EFO_0003756)
# Output:    results/tables/Table6_eQTL_GWAS_overlap.csv
#            results/tables/Table7_overlap_DEG_info.csv
#            results/figures/Fig_eQTL_GWAS_overlap.png
#
# External data required (not included — download manually):
#   GTEx_Analysis_v11_eQTL/Brain_Cortex.v11.eGenes.txt.gz
#   GTEx_Analysis_v11_eQTL/Brain_Cortex.v11.eQTLs.signif_pairs.parquet
#   GTEx_Analysis_v11_eQTL/EFO_0003756_associations_export.tsv
#   Download from: https://gtexportal.org/home/downloads/adult-gtex/qtl
# =============================================================================

# ----- Load Libraries --------------------------------------------------------
library(arrow)
library(dplyr)
library(ggplot2)

# ----- Set Working Directory -------------------------------------------------
setwd("~/Desktop/BMI 5332_R/Project")

# ----- Load DEG Results ------------------------------------------------------
deg_data<- readRDS("Data/Processed/DEG_results.rds")
sig_DEGs<- deg_data$sig_DEGs


# =============================================================================
# STEP 1: Load GTEx eGenes (Brain Cortex, v11)
# =============================================================================
egenes <- read.table(
  "Data/GTEx_Analysis_v11_eQTL/Brain_Cortex.v11.eGenes.txt.gz",
  header = TRUE, sep = "\t")

# Filter for significant eGenes only
egenes_sig <- egenes[egenes$qval <= 0.05, ]
egenes_sig$gene_name <- trimws(egenes_sig$gene_name)

cat("Total eGenes tested: ", nrow(egenes), "\n")
cat("Significant eGenes: ", nrow(egenes_sig), "\n")


# =============================================================================
# STEP 2: Load Significant eQTL Variant-Gene Pairs
# =============================================================================
signif_pairs <- read_parquet(
  "Data/GTEx_Analysis_v11_eQTL/Brain_Cortex.v11.eQTLs.signif_pairs.parquet")

cat("Significant eQTL pairs: ", nrow(signif_pairs), "\n")


# =============================================================================
# STEP 3: Intersect DEGs with Brain Cortex eGenes
# =============================================================================
deg_egenes <- egenes_sig[egenes_sig$gene_name %in% sig_DEGs$GeneSymbol, ]

cat("Total significant eGenes in brain cortex: ", nrow(egenes_sig), "\n")
cat("DEGs that are also eGenes: ", nrow(deg_egenes), "\n")
cat("DEGs not found as eGenes: ", sum(!sig_DEGs$GeneSymbol %in% egenes_sig$gene_name), "\n")
cat("Percentage of DEGs that are eGenes: ", round(nrow(deg_egenes) / nrow(sig_DEGs) * 100, 1), "%\n")

# Get eQTL pairs for DEG eGenes
deg_egene_ids  <- deg_egenes$gene_id
deg_eqtl_pairs <- signif_pairs[signif_pairs$phenotype_id %in% deg_egene_ids, ]

cat("Total eQTL pairs for DEG eGenes: ", nrow(deg_eqtl_pairs), "\n")
cat("Unique eQTL variants: ", length(unique(deg_eqtl_pairs$variant_id)), "\n")
cat("Unique eGenes with pairs: ", length(unique(deg_eqtl_pairs$phenotype_id)), "\n")


# =============================================================================
# STEP 4: Load ASD GWAS Variants
# =============================================================================
gwas_associations <- read.table(
  "Data/GTEx_Analysis_v11_eQTL/EFO_0003756_associations_export.tsv",
  header = TRUE, sep = "\t", quote = "", fill = TRUE)

gwas_associations$traitName <- gsub('"', '', gwas_associations$traitName)
gwas_asd <- gwas_associations[gwas_associations$traitName == "Autism spectrum disorder", ]

# Parse rsIDs and genomic positions
gwas_asd$rsID <- gsub("-.*", "", gwas_asd$riskAllele)
gwas_asd$chr <- paste0("chr", gsub(":.*", "", gwas_asd$locations))
gwas_asd$pos <- as.numeric(gsub(".*:", "", gwas_asd$locations))

cat("Unique ASD GWAS rsIDs: ", length(unique(gwas_asd$rsID)), "\n")
cat("Unique ASD GWAS positions: ", length(unique(gwas_asd$locations)), "\n")


# =============================================================================
# STEP 5: Window-Based Overlap (±500kb)
# =============================================================================
# Extract chr and position from GTEx variant IDs
deg_eqtl_pairs$chr <- gsub("_.*", "", deg_eqtl_pairs$variant_id)
deg_eqtl_pairs$pos <- as.numeric(
  sapply(strsplit(deg_eqtl_pairs$variant_id, "_"), `[`, 2))

cat("eQTL variants with positions extracted: ", nrow(deg_eqtl_pairs), "\n")

# Merge and apply ±500kb window
overlap_window <- merge(
  deg_eqtl_pairs,
  gwas_asd[, c("rsID", "chr", "pos", "mappedGenes", "pValue")],
  by = "chr") %>%
  filter(abs(pos.x - pos.y) <= 500000)

overlap_window$GeneSymbol <- deg_egenes$gene_name[
  match(overlap_window$phenotype_id, deg_egenes$gene_id)]

cat("Window-based matches: ", nrow(overlap_window), "\n")
cat("Unique DEG eGenes in window: ", length(unique(overlap_window$phenotype_id)), "\n")
cat("Unique ASD GWAS variants matched: ", length(unique(overlap_window$rsID)), "\n")

head(overlap_window[, c("GeneSymbol", "variant_id", "rsID",
                        "mappedGenes", "slope", "pValue")], 5)


# =============================================================================
# STEP 6: Summarize Overlap Results
# =============================================================================
overlap_summary <- overlap_window %>%
  group_by(GeneSymbol, rsID, mappedGenes, pValue) %>%
  summarise(
    n_eQTL_variants  = n(),
    mean_eQTL_slope  = round(mean(slope), 4),
    min_pval_nominal = min(pval_nominal),
    .groups          = "drop") %>%
  arrange(pValue)

overlap_genes <- unique(overlap_window$GeneSymbol)
overlap_deg_info <- sig_DEGs[sig_DEGs$GeneSymbol %in% overlap_genes,
                             c("GeneSymbol", "logFC", "adj.P.Val", "category")]
overlap_deg_info <- overlap_deg_info[order(overlap_deg_info$adj.P.Val), ]

cat("\n--- Final Summary ---\n")
cat("Total DEGs: ", nrow(sig_DEGs), "\n")
cat("DEGs that are eGenes: ", nrow(deg_egenes), "\n")
cat("DEGs with eQTLs near ASD GWAS variants: ", length(overlap_genes), "\n")
cat("Percentage of DEGs: ", round(length(overlap_genes) / nrow(sig_DEGs) * 100, 1), "%\n")
cat("Percentage of DEG eGenes: ", round(length(overlap_genes) / nrow(deg_egenes) * 100, 1), "%\n")

cat("\nDEG genes overlapping with ASD GWAS variants:\n")
print(overlap_genes)
print(overlap_deg_info)

# Save tables
write.csv(overlap_summary,  "Results/Table6_eQTL_GWAS_overlap.csv",  row.names = FALSE)
write.csv(overlap_deg_info, "Results/Table7_overlap_DEG_info.csv",   row.names = FALSE)


# =============================================================================
# STEP 7: Lollipop Plot — DEGs Overlapping ASD GWAS Variants
# =============================================================================
overlap_deg_info$GeneSymbol <- factor(
  overlap_deg_info$GeneSymbol,
  levels = overlap_deg_info$GeneSymbol[order(overlap_deg_info$logFC)])

png("Results/Fig9_eQTL_GWAS_overlap.png", width = 10, height = 8, units = "in", res = 300)

ggplot(overlap_deg_info,
       aes(x = logFC, y = GeneSymbol, color = category)) +
  geom_segment(aes(x = 0, xend = logFC, y = GeneSymbol, yend = GeneSymbol),
               color = "grey60", linewidth = 0.8) +
  geom_point(size = 4) +
  scale_color_manual(values = c("Upregulated"   = "red",
                                "Downregulated" = "blue"),
                     name   = "Direction") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(title    = "DEGs with eQTLs Overlapping ASD GWAS Variants (\u00b1500kb)",
       subtitle = paste0("22 genes | 35 DEG-GWAS pairs | 6.2% of DEG eGenes"),
       x        = "log2 Fold Change (MBD6het vs Control)",
       y        = "Gene Symbol") +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.y   = element_text(size = 10),
    axis.text.x   = element_text(size = 10),
    legend.title  = element_text(size = 10, face = "bold"))

dev.off()
