# =============================================================================
# Script 03: Principal Component Analysis (PCA)
# Project:   Cross-ancestry T2D PRS Transferability (EUR -> CSA)
# Data:      1KGP Phase 3 QC-passed genotypes (EUR and SAS)
# =============================================================================
# OVERVIEW:
# This script performs principal component analysis (PCA) on the QC-passed
# 1KGP genotype data. PCA serves two purposes in this project:
#
#   1. ANCESTRY VERIFICATION: Confirm that EUR and SAS samples cluster
#      as expected and that no population outliers exist that could
#      confound downstream analyses. Samples deviating >6 SD from their
#      superpopulation centroid in PC space are removed.
#
#   2. COVARIATE GENERATION: The top 10 PCs are saved as covariates for
#      use in PRS association models in Script 06. Including PCs as
#      covariates controls for residual population stratification within
#      each superpopulation.
#
# OUTPUTS:
#   data/1KGP/pca_combined.rds        - PCA object (scores + loadings)
#   data/1KGP/pca_covariates_EUR.txt  - Top 10 PCs for EUR samples
#   data/1KGP/pca_covariates_SAS.txt  - Top 10 PCs for SAS samples
#   data/1KGP/samples_EUR_final.txt   - EUR samples passing outlier QC
#   data/1KGP/samples_SAS_final.txt   - SAS samples passing outlier QC
#   results/figures/pca_EUR_SAS.png   - Manuscript PCA figure
#   results/QC/pca_outliers.csv       - Outlier sample report
# =============================================================================

library(bigsnpr)
library(bigstatsr)
library(data.table)
library(ggplot2)

# -----------------------------------------------------------------------------
# Load QC-passed genotype data
# -----------------------------------------------------------------------------
# Attach the original bigSNP object
obj.bigsnp <- snp_attach("data/1KGP/1000G_phase3_common_norel.rds")
G <- obj.bigsnp$genotypes
map <- obj.bigsnp$map
fam <- obj.bigsnp$fam

# Load population panel and QC-passed sample lists
pop_panel <- fread("data/1KGP/integrated_call_samples_v3.20130502.ALL.panel")
fam_EUR <- fread("data/1KGP/1KGP_EUR_QC.fam")
fam_SAS <- fread("data/1KGP/1KGP_SAS_QC.fam")

# Get indices of QC-passed EUR and SAS samples in the full bigSNP object
idx_EUR_qc <- which(fam$sample.ID %in% fam_EUR$V2)
idx_SAS_qc <- which(fam$sample.ID %in% fam_SAS$V2)
idx_both   <- c(idx_EUR_qc, idx_SAS_qc)

cat("EUR samples loaded:", length(idx_EUR_qc), "\n")     #503
cat("SAS samples loaded:", length(idx_SAS_qc), "\n")     #484
cat("Total samples for PCA:", length(idx_both), "\n\n")  #987

# Load QC-passed SNP list (intersection of EUR and SAS passing SNPs)
snps_both <- fread("data/1KGP/snps_pass_both.txt", header = FALSE)$V1

# Get SNP indices in the bigSNP map
idx_snps <- which(map$marker.ID %in% snps_both)
cat("SNPs for PCA:", length(idx_snps), "\n\n") #1422959

# -----------------------------------------------------------------------------
# LD pruning
# -----------------------------------------------------------------------------
# snp_indLRLDR() identifies SNPs to REMOVE that are in long-range LD
# regions known to distort PCA These regions are hardcoded in bigsnpr 
# based on published exclusion lists (Price et al. 2008, Anderson et al. 2010).
# snp_clumping() then performs LD-based pruning keeping one SNP per
# LD block (r2 < 0.2 within 1000 SNP window) within EUR samples.

### Identify SNPs to exclude:
# All SNP indices NOT in QC-passed set
idx_all   <- seq_len(nrow(map))
idx_excl_nonqc <- setdiff(idx_all, idx_snps)
cat("SNPs excluded (failed QC): ", length(idx_excl_nonqc), "\n") #241893

# Long-range LD regions within QC-passed SNPs only
lrldr_local <- snp_indLRLDR(
  infos.chr = map$chromosome[idx_snps],
  infos.pos = map$physical.pos[idx_snps]
)
idx_excl_lrldr <- idx_snps[lrldr_local]
cat("SNPs excluded (long-range LD regions):", length(idx_excl_lrldr), "\n") #39578

# Combine all exclusions
idx_excl_all <- union(idx_excl_nonqc, idx_excl_lrldr)
cat("Total SNPs excluded from clumping:   ", length(idx_excl_all), "\n\n")

### LD clumping
ind_keep <- snp_clumping(
  G,
  infos.chr = map$chromosome,
  ind.row   = idx_EUR_qc,        # EUR samples for LD estimation
  exclude   = idx_excl_all,      # skip non-QC and long-range LD SNPs
  thr.r2    = 0.2,               # r2 threshold for LD pruning
  size      = 1000,              # window size in SNPs
  ncores    = nb_cores()         # use all available cores
)


# -----------------------------------------------------------------------------
# Run PCA
# -----------------------------------------------------------------------------
# big_SVD() performs a truncated singular value decomposition (SVD) on
# the genotype matrix, which is equivalent to PCA. We compute the top
# 20 PCs to allow flexible use downstream -- the top 10 are used as
# covariates but we compute 20 to check for any structure beyond PC10.

pca <- big_SVD(
  G,
  fun.scaling = snp_scaleBinom(),   # center and scale genotypes
  ind.row     = idx_both,           # EUR + SAS samples only
  ind.col     = ind_keep,           # LD-pruned QC-passed SNPs
  k           = 20,                 # compute top 20 PCs
)

# -----------------------------------------------------------------------------
# Prepare PCA results dataframe
# -----------------------------------------------------------------------------
# Attach population labels to PC scores for visualization and outlier
# detection. 
scores_df <- data.table(
  sample    = fam$sample.ID[idx_both],
  super_pop = pop_panel$super_pop[match(fam$sample.ID[idx_both],
                                        pop_panel$sample)],
  as.data.table(pca$u)
)

# Rename PC columns for clarity
pc_cols <- paste0("PC", 1:20)
setnames(scores_df, paste0("V", 1:20), pc_cols)

# PCA scores dataframe preview
print(head(scores_df[, 1:7]))

# -----------------------------------------------------------------------------
# Ancestry verification
# -----------------------------------------------------------------------------
sep_check <- scores_df[, .(
  mean_PC1 = round(mean(PC1), 4),
  mean_PC2 = round(mean(PC2), 4),
  sd_PC1   = round(sd(PC1),   4),
  n        = .N
), by = super_pop]

# Report mean PC1 and PC2 per superpopulation to confirm separation.
print(sep_check)

# -----------------------------------------------------------------------------
# Outlier detection and removal (>6 SD from superpopulation centroid)
# -----------------------------------------------------------------------------
outliers <- c()

for(pop in c("EUR", "SAS")) {
  
  pop_idx    <- which(scores_df$super_pop == pop)
  
  # Fix: use .SD with .SDcols to correctly subset data.table columns
  pop_scores <- as.matrix(scores_df[pop_idx, .SD, .SDcols = pc_cols[1:10]])
  
  # Z-scores relative to superpopulation centroid
  pop_means  <- colMeans(pop_scores)
  pop_sds    <- apply(pop_scores, 2, sd)
  z_scores   <- sweep(
    sweep(pop_scores, 2, pop_means, "-"),
    2, pop_sds, "/")
  
  # Flag samples exceeding 6 SD on any PC
  outlier_flag <- apply(abs(z_scores), 1, max) > 6
  outlier_idx  <- pop_idx[outlier_flag]
  
  cat(pop, "outliers detected:", sum(outlier_flag), "\n")
  if(sum(outlier_flag) > 0) {
    cat("  Outlier samples:", scores_df$sample[outlier_idx], "\n")
  }
  outliers <- c(outliers, outlier_idx)
}

# Remove outliers
if(length(outliers) > 0) {
  scores_df_clean <- scores_df[-outliers]
  cat("\nTotal outliers removed:", length(outliers), "\n")
} else {
  scores_df_clean <- scores_df
  cat("\nNo outliers detected -- all samples retained.\n")
}

# Save outlier report
outlier_report <- data.table(
  sample    = scores_df$sample[outliers],
  super_pop = scores_df$super_pop[outliers]
)
fwrite(outlier_report, "results/QC/pca_outliers.csv")

# -----------------------------------------------------------------------------
# Save Outputs
# -----------------------------------------------------------------------------

# final sample lists
eur_final <- scores_df_clean[super_pop == "EUR",
                             .(FID = sample, IID = sample)]
sas_final <- scores_df_clean[super_pop == "SAS",
                             .(FID = sample, IID = sample)]

fwrite(eur_final, "data/1KGP/samples_EUR_final.txt", sep = "\t")
fwrite(sas_final, "data/1KGP/samples_SAS_final.txt", sep = "\t")


# Save Top 10 PCA covariates
pcs_EUR <- scores_df_clean[super_pop == "EUR", c("sample", paste0("PC", 1:10)), with = FALSE]
pcs_SAS <- scores_df_clean[super_pop == "SAS", c("sample", paste0("PC", 1:10)), with = FALSE]

fwrite(pcs_EUR, "data/1KGP/pca_covariates_EUR.txt", sep = "\t")
fwrite(pcs_SAS, "data/1KGP/pca_covariates_SAS.txt", sep = "\t")

saveRDS(list(pca = pca, scores = scores_df_clean), "data/1KGP/pca_combined.rds")

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------
pop_colors <- c("EUR" = "#4a90d9", "SAS" = "#e07b39")

# PC1 vs PC2 
p1 <- ggplot(scores_df_clean,
             aes(x = PC1, y = PC2, color = super_pop)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = pop_colors,
                     name   = "Superpopulation",
                     labels = c("EUR" = "European (EUR)",
                                "SAS" = "South Asian (SAS)")) +
  labs(
    x = paste0("PC1 (",
               round(pca$d[1]^2 / sum(pca$d^2) * 100, 1),
               "% variance explained)"),
    y = paste0("PC2 (",
               round(pca$d[2]^2 / sum(pca$d^2) * 100, 1),
               "% variance explained)")
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(color = "gray40", size = 10),
        panel.grid.minor = element_blank())

ggsave("results/figures/fig1_pca_PC1_PC2.png", plot = p1, width = 7, height = 6, dpi = 300)

# PC3 vs PC4
p2 <- ggplot(scores_df_clean,
             aes(x = PC3, y = PC4, color = super_pop)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = pop_colors,
                     name   = "Superpopulation") +
  labs(
    title    = "PCA of 1KGP Phase 3: PC3 vs PC4",
    subtitle = "Within-population substructure",
    x = paste0("PC3 (",
               round(pca$d[3]^2 / sum(pca$d^2) * 100, 1),
               "% variance explained)"),
    y = paste0("PC4 (",
               round(pca$d[4]^2 / sum(pca$d^2) * 100, 1),
               "% variance explained)")
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(color = "gray40", size = 10),
        panel.grid.minor = element_blank())
# ggsave("results/QC/pca_PC3_PC4.png", plot = p2, width = 7, height = 6, dpi = 300)

# Scree plot 
var_explained <- data.table(
  PC     = 1:20,
  VarExp = pca$d[1:20]^2 / sum(pca$d^2) * 100
)

p3 <- ggplot(var_explained, aes(x = PC, y = VarExp)) +
  geom_col(fill = "#4a90d9", alpha = 0.8) +
  geom_line(color = "gray40", linewidth = 0.8) +
  geom_point(color = "gray40", size = 2) +
  scale_x_continuous(breaks = 1:20) +
  labs(
    title    = "PCA Scree Plot — Variance Explained per PC",
    subtitle = "Top 20 principal components",
    x        = "Principal Component",
    y        = "Variance Explained (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(color = "gray40", size = 10),
        panel.grid.minor = element_blank())
# ggsave("results/QC/pca_scree.png", plot = p3, width = 8, height = 5, dpi = 300)


# -----------------------------------------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------------------------------------
# SNPs used for PCA (LD-pruned): 225428
# Outliers removed: 6
# Final EUR samples: 503
# Final SAS samples: 478
# PC1 variance explained: 50.33%
# PC2 variance explained: 6.24%
