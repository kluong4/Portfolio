# =============================================================================
# Script 01: Quality Control of Pan-UKB Summary Statistics
# Project:   Cross-ancestry T2D PRS Transferability
# Data:      Pan-UK Biobank E11 (ICD-10 Type 2 Diabetes, both sexes)
# =============================================================================
# OVERVIEW:
# This script performs quality control on GWAS summary statistics downloaded
# from the Pan-UK Biobank. Only EUR and CSA ancestries are retained for primary
# analysis, as these are the only two ancestry groups passing Pan-UKB internal
# QC filters for the E11 phenotype (pops_pass_qc: EUR n_cases=22,634;
# CSA n_cases=1,662). AFR, EAS, and MID columns are retained in the full
# sumstats object for reference but are not used in primary analyses.
#
# OUTPUTS:
#   sumstats_E11_QC.rds                       - Full QC-passed sumstats
#   gws_EUR_loci_annotated.rds                - EUR GWS loci with CSA data
#   sumstats_EUR_formatted.rds                - EUR-only formatted for PRS
#   sumstats_CSA_formatted.rds                - CSA-only formatted for PRS
#   ldsc_EUR.sumstats.gz                      - EUR formatted for LDSC
#   ldsc_CSA.sumstats.gz                      - CSA formatted for LDSC
#   directional_consistency_by_zscore.csv.    - Key QC finding
#   directional_concordance_zscore.png        - Manuscript figure
# =============================================================================

library(data.table)
library(ggplot2)

# -----------------------------------------------------------------------------
# Find T2D Phenotype Code from phenotype_manifest file
# -----------------------------------------------------------------------------

# Download the phenotype manifest to find T2D phenotype code
system("curl -o data/sumstats/phenotype_manifest.tsv.bgz https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz")

# Find the T2D phenotype code in R
manifest <- fread("data/sumstats/phenotype_manifest.tsv.bgz")

# Search for T2D
t2d <- manifest[grepl("type 2 diabetes|T2D|E11", description, ignore.case = TRUE)]
print(t2d[, .(phenocode, description, pops, filename, aws_path)])

print(t2d[1, .(
  phenocode,
  pops_pass_qc, # EUR and CSA QC pass
  n_cases_AFR, n_controls_AFR,
  n_cases_CSA, n_controls_CSA,
  n_cases_EAS, n_controls_EAS,
  n_cases_EUR, n_controls_EUR,
  n_cases_MID, n_controls_MID
)])

# Download E11 and index file for all 5 available ancestries
system("curl -L -o data/sumstats/icd10-E11-both_sexes.tsv.bgz https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/icd10-E11-both_sexes.tsv.bgz")
system("curl -L -o data/sumstats/icd10-E11-both_sexes.tsv.bgz.tbi https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files_tabix/icd10-E11-both_sexes.tsv.bgz.tbi")


# -----------------------------------------------------------------------------
# Load summary statistics
# -----------------------------------------------------------------------------
# NOTE: p-values are stored as -log10(p) to avoid floating point underflow
# at very small p-values. Genome-wide significance (p < 5e-8) corresponds
# to -log10(p) > 7.3.

sumstats <- fread("data/sumstats/icd10-E11-both_sexes.tsv.bgz")

cat("Variants loaded:", nrow(sumstats), "\n") # 28987534
cat("Columns:", paste(colnames(sumstats), collapse = ", "), "\n\n")


# -----------------------------------------------------------------------------
# Remove low confidence variants
# -----------------------------------------------------------------------------
# Pan-UKB flags variants as low_confidence = TRUE when:
#   - Alternate allele count in cases <= 3
#   - Alternate allele count in controls <= 3
#   - Minor allele count (cases + controls combined) <= 20
# These thresholds identify variants where effect size estimates are
# unreliable due to extremely rare alleles. 
# Variants remaining: 13060391

sumstats <- sumstats[low_confidence_EUR == FALSE & low_confidence_CSA == FALSE]


# -----------------------------------------------------------------------------
# Remove strand-ambiguous SNPs
# -----------------------------------------------------------------------------
# # A/T and C/G SNPs are strand-ambiguous: it is impossible to determine
# which strand was genotyped without allele frequency information, making
# allele alignment between the summary statistics and target genotype data
# (1KGP) unreliable. Removing these prevents strand-flip errors during
# SNP harmonization in downstream PRS construction scripts.
# Variants remaining: 11239095

ambiguous <- with(sumstats,
                  (ref == "A" & alt == "T") | (ref == "T" & alt == "A") |
                    (ref == "C" & alt == "G") | (ref == "G" & alt == "C"))
sumstats <- sumstats[!ambiguous]


# -----------------------------------------------------------------------------
# Verify uniqueness with variant ID
# -----------------------------------------------------------------------------
# Create variant ID
sumstats[, varID := paste(chr, pos, ref, alt, sep = ":")]

# remove duplicates (keep highest -log10(p) EUR) - No duplicates
setorder(sumstats, -neglog10_pval_EUR)
sumstats <- sumstats[!duplicated(varID)]


### QC Summary
# SNPs after QC: 11239095
# EUR GWS loci (p < 5e-8): 8789 
# CSA GWS loci (p < 5e-8): 33


# -----------------------------------------------------------------------------
# Cross-ancestry directional concordance analysis
# -----------------------------------------------------------------------------
# At EUR genome-wide significant loci, we assess whether EUR and CSA effect
# sizes (beta coefficients) point in the same direction. For a trait with
# shared genetic architecture across ancestries (expected for T2D), we
# anticipate high directional concordance at well-powered loci. Low
# concordance at weak CSA signals is expected and reflects statistical noise
# (near-zero CSA betas have effectively random signs) rather than true
# biological heterogeneity. This analysis directly informs interpretation
# of PRS transferability findings.

# Extract EUR GWS loci that also have CSA data
gws_EUR <- sumstats[neglog10_pval_EUR > 7.3 & !is.na(beta_CSA) & !is.na(se_CSA)]

cat("EUR GWS loci with CSA data:", nrow(gws_EUR), "\n")
cat("EUR GWS loci also GWS in CSA:", 
    sum(gws_EUR$neglog10_pval_CSA > 7.3, na.rm = TRUE), "\n")
cat("Overall same direction EUR & CSA:", 
    sum(sign(gws_EUR$beta_EUR) == sign(gws_EUR$beta_CSA), na.rm = TRUE),
    "out of", sum(!is.na(gws_EUR$beta_CSA)), 
    sprintf("(%.1f%%)\n\n",
            mean(sign(gws_EUR$beta_EUR) == sign(gws_EUR$beta_CSA), 
                 na.rm = TRUE) * 100))

# Compute CSA Z-score magnitude to stratify by statistical power.
gws_EUR[, csa_absz := abs(beta_CSA / se_CSA)]

gws_EUR[, z_bin := cut(csa_absz, 
                       breaks = c(0, 1, 2, 3, 4, Inf),
                       labels = c("<1", "1-2", "2-3", "3-4", ">4"))]

# Compute directional consistency and 95% CI per Z-score bin.
consistency_by_bin <- gws_EUR[, .(
  n        = .N,
  same_dir = sum(sign(beta_EUR) == sign(beta_CSA), na.rm = TRUE)
), by = z_bin][order(z_bin)]

consistency_by_bin[, pct      := same_dir / n * 100]
consistency_by_bin[, lower_ci := (same_dir / n -
                                    1.96 * sqrt((same_dir/n * (1 - same_dir/n)) / n)) * 100]
consistency_by_bin[, upper_ci := (same_dir / n +
                                    1.96 * sqrt((same_dir/n * (1 - same_dir/n)) / n)) * 100]
consistency_by_bin[, upper_ci := pmin(upper_ci, 100)]

cat("Directional concordance by CSA Z-score bin:\n")
print(consistency_by_bin)

# Annotate gws_EUR with Z-score bins for use in downstream scripts
# (functional annotation and effect size concordance analyses)
saveRDS(gws_EUR, "data/sumstats/outputs/gws_EUR_loci_annotated.rds")


# -----------------------------------------------------------------------------
# Directional concordance figure
# -----------------------------------------------------------------------------
# support interpretation that low overall concordance is power-driven
ggplot(consistency_by_bin,
       aes(x = z_bin, y = pct, group = 1)) +
  # Reference lines
  geom_hline(yintercept = 50,  linetype = "dashed",
             color = "gray60", linewidth = 0.6) +
  geom_hline(yintercept = 100, linetype = "dotted",
             color = "gray40", linewidth = 0.6) +
  # 95% confidence interval ribbon
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              fill = "#4a90d9", alpha = 0.2) +
  geom_line(color  = "#4a90d9", linewidth = 1.2) +
  geom_point(color = "#4a90d9", size = 4) +
  # Label each point with percentage and sample size
  geom_text(aes(label = paste0(round(pct, 1), "%\n(n=", n, ")")),
            vjust = -1.2, size = 3.2, color = "gray30") +
  scale_y_continuous(limits = c(40, 108),
                     breaks = c(50, 60, 70, 80, 90, 100),
                     labels = paste0(c(50, 60, 70, 80, 90, 100), "%")) +
  labs(
    title    = "Directional Concordance of T2D Effect Sizes: EUR vs CSA",
    subtitle = "Concordance by CSA GWAS Z-score magnitude at EUR GWS loci",
    x        = "CSA |Z-score| bin",
    y        = "% Same Effect Direction (EUR vs CSA)",
    caption  = paste0("Dashed line = 50% (chance); shaded band = 95% CI\n",
                      "EUR GWS loci: -log10(p) > 7.3 (p < 5e-8); ",
                      "n = ", nrow(gws_EUR), " loci")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "gray40", size = 10),
    panel.grid.minor = element_blank()
  )

ggsave("results/figures/directional_concordance_zscore.png",
       width = 8, height = 5.5, dpi = 300)


# -----------------------------------------------------------------------------
# Save all downstream outputs
# -----------------------------------------------------------------------------

# Full QC-passed sumstats
saveRDS(sumstats, "data/sumstats/outputs/sumstats_E11_QC.rds")

# EUR-formatted sumstats for PRS construction
sumstats_EUR <- sumstats[
  !is.na(beta_EUR) & !is.na(se_EUR),
  .(chr, pos, ref, alt, varID,
    beta = beta_EUR, se = se_EUR,
    neglog10_pval = neglog10_pval_EUR,
    af_cases = af_cases_EUR, af_controls = af_controls_EUR)
]
saveRDS(sumstats_EUR, "data/sumstats/outputs/sumstats_EUR_formatted.rds")

# CSA-formatted sumstats for PRS construction
sumstats_CSA <- sumstats[
  !is.na(beta_CSA) & !is.na(se_CSA),
  .(chr, pos, ref, alt, varID,
    beta = beta_CSA, se = se_CSA,
    neglog10_pval = neglog10_pval_CSA,
    af_cases = af_cases_CSA, af_controls = af_controls_CSA)
]
saveRDS(sumstats_CSA, "data/sumstats/outputs/sumstats_CSA_formatted.rds")

# LDSC-formatted files for genetic correlation
ldsc_EUR <- sumstats[!is.na(beta_EUR) & !is.na(se_EUR),
                     .(SNP = varID, A1 = alt, A2 = ref,
                       Z = beta_EUR / se_EUR,
                       N = af_cases_EUR + af_controls_EUR)]
fwrite(ldsc_EUR, "data/sumstats/outputs/ldsc_EUR.sumstats.gz", sep = "\t")

ldsc_CSA <- sumstats[!is.na(beta_CSA) & !is.na(se_CSA),
                     .(SNP = varID, A1 = alt, A2 = ref,
                       Z = beta_CSA / se_CSA,
                       N = af_cases_CSA + af_controls_CSA)]
fwrite(ldsc_CSA, "data/sumstats/outputs/ldsc_CSA.sumstats.gz", sep = "\t")

# Directional consistency table
fwrite(consistency_by_bin, "data/sumstats/outputs/directional_consistency_by_zscore.csv")
saveRDS(consistency_by_bin, "data/sumstats/outputs/directional_consistency_by_zscore.rds")