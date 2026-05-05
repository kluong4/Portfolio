# =============================================================================
# Script 01: Data Loading, Preprocessing, and Quality Control
# Project:   MBD6 Haploinsufficiency and ASD - Transcriptomic Analysis
# Course:    BMI 5332: Statistical Analysis of Genomic Data
#            UT Health Houston, Spring 2026
# Data:      GEO Accession: GSE314093
#            Platform: GPL13497 (Agilent-028004 SurePrint G3 Human GE 8x60K)
# =============================================================================

# ----- Load Libraries --------------------------------------------------------
library(GEOquery)
library(limma)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

# ----- Set Working Directory -------------------------------------------------
setwd("~/Desktop/BMI 5332_R/Project")

# =============================================================================
# STEP 1: Load Raw Files
# =============================================================================
# Create a data frame that stores metadata for the microarray files
targets <- data.frame(
  FileName = file.path("Data/GSE314093", c(
    "GSM9381448_CTL_1.txt.gz",
    "GSM9381450_CTL_2.txt.gz",
    "GSM9381449_MBD6hetero_1.txt.gz",
    "GSM9381451_MBD6hetero_2.txt.gz"
  )),
  Condition = c("CTL", "CTL", "MBD6hetero", "MBD6hetero"),
  stringsAsFactors = FALSE
)

# Read the raw Agilent files
raw <- read.maimages(
  targets$FileName,
  source = "agilent",
  green.only = TRUE,
  other.columns = "gIsWellAboveBG"
)

# Sanity checks
dim(raw)           # 44495 probes x 4 samples
head(raw$genes)    # Row, Col, ControlType, ProbeName, SystematicName columns
head(raw$E)        # expression matrix (raw intensities)
colnames(raw$E)    # confirm sample names


# =============================================================================
# STEP 2: Background Correction and Normalization
# =============================================================================
png(file.path("Results", "Fig1_ExpressionValues.png"), width = 1200, height = 600, res = 150)

par(mfrow = c(1, 2), mar = c(7, 5, 4, 2))

# Before normalization
boxplot(log2(raw$E),
        main  = "(A) Raw Expression Values",
        ylab  = "log2(Expression)",
        col   = c("#A465BA", "#A465BA", "#1B9E77", "#1B9E77"),
        names = c("CTL_1", "CTL_2", "MBD6het_1", "MBD6het_2"), las = 2)

# Background correction (normexp method)
bg <- backgroundCorrect(raw, method = "normexp")

# Quantile normalization across arrays
norm <- normalizeBetweenArrays(bg, method = "quantile")

# After normalization — distributions should now overlap
boxplot(log2(norm$E),
        main  = "(B) Normalized Expression Values",
        ylab  = "log2(Expression)",
        col   = c("#A465BA", "#A465BA", "#1B9E77", "#1B9E77"),
        names = c("CTL_1", "CTL_2", "MBD6het_1", "MBD6het_2"), las = 2)

dev.off()

# =============================================================================
# STEP 3: Probe Filtering
# =============================================================================
# Keep probes well above background in at least 2 samples
keep <- rowSums(norm$other$gIsWellAboveBG == 1) >= 2
norm_filtered <- norm[keep, ]

cat("Total probes before filtering: ", 44495, "\n")
cat("Probes after filtering:        ", 30691, "\n")
cat("Probes removed:                ", 44495 - 30691, "\n")

# Fetch platform annotation and add gene symbols
gpl <- getGEO("GPL13497")
annotation <- Table(gpl)

norm_filtered$genes <- merge(
  norm_filtered$genes,
  annotation[, c("ID", "GENE_SYMBOL")],
  by.x = "ProbeName",
  by.y = "ID",
  all.x = TRUE
)

# Remove probes with no gene symbol
has_symbol    <- !is.na(norm_filtered$genes$GENE_SYMBOL) &
  norm_filtered$genes$GENE_SYMBOL != ""
norm_filtered <- norm_filtered[has_symbol, ]

# Rename column for consistency
colnames(norm_filtered$genes)[
  colnames(norm_filtered$genes) == "GENE_SYMBOL"] <- "GeneSymbol"

cat("Probes before symbol filtering: ", 30691, "\n")
cat("Probes after symbol filtering:  ", 28280, "\n")
cat("Probes removed:                 ", 30691 - 28280, "\n")


# =============================================================================
# STEP 4: Quality Control Plots
# =============================================================================

# --- PCA ---------------------------------------------------------------------
png(file.path("Results", "Fig2_PCAplot.png"), width = 1000, height = 700, res = 150)

pca <- prcomp(t(log2(norm_filtered$E + 0.1)), scale. = TRUE)
pc1_var <- round(summary(pca)$importance[2, 1] * 100, 1)
pc2_var <- round(summary(pca)$importance[2, 2] * 100, 1)

plot(pca$x[, 1], pca$x[, 2],
     col  = c("#A465BA", "#A465BA", "#1B9E77", "#1B9E77"),
     pch  = 19, cex  = 2.5,
     xlab = paste0("PC1 (", pc1_var, "%)"),
     ylab = paste0("PC2 (", pc2_var, "%)"),
     main = "Principal Component Analysis of Normalized Expression Data")

text(pca$x[, 1], pca$x[, 2],
     labels = c("CTL_1", "CTL_2", "MBD6het_1", "MBD6het_2"),
     pos    = c(4, 4, 2, 2), cex = 0.9, offset = 1.0)

legend("topright", legend = c("Control", "MBD6het"),
       col = c("#A465BA", "#1B9E77"), pch = 19, cex = 0.9)

dev.off()

# --- Sample Correlation Heatmap ----------------------------------------------
png(file.path("Results", "Fig3_SampleCorrelationHeatmap.png"), width = 600, height = 500)

cor_matrix <- cor(log2(norm_filtered$E + 0.1))
rownames(cor_matrix) <- colnames(cor_matrix) <-
  c("CTL_1", "CTL_2", "MBD6het_1", "MBD6het_2")

cor_melt         <- melt(cor_matrix)
colnames(cor_melt) <- c("Sample1", "Sample2", "Correlation")

ggplot(cor_melt, aes(x = Sample1, y = Sample2, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Correlation, 3)), size = 5, color = "black") +
  scale_fill_gradient(low = "white", high = "gold", limits = c(0.97, 1.00)) +
  labs(title = "Sample Correlation Heatmap", x = "", y = "",
       fill  = "Pearson\nCorrelation") +
  theme_minimal() +
  theme(
    plot.title  = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    panel.grid  = element_blank()
  )

dev.off()

# --- MA Plots ----------------------------------------------------------------
png(file.path("Results", "Fig4_MAplots.png"), width = 600, height = 500)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (i in 1:4) {
  plotMA(norm_filtered, array = i,
         main = c("CTL_1", "CTL_2", "MBD6het_1", "MBD6het_2")[i],
         xlab = "Average log-expression (A)",
         ylab = "Log2 Fold Change (M)")
  abline(h = 0, col = "red")
}

dev.off()

# =============================================================================
# SAVE: norm_filtered object for use in downstream scripts
# =============================================================================
saveRDS(norm_filtered, file = "Data/Processed/norm_filtered.rds")
