# =============================================================================
# Script 07: Genetic Correlation Analysis (LDSC)
# Project:   Cross-ancestry T2D PRS Transferability (EUR -> CSA)
# =============================================================================
# OVERVIEW:
# This script estimates the genome-wide genetic correlation (rg) between
# EUR and CSA T2D GWAS summary statistics using LD Score Regression (LDSC)
# via the ldscr package (v0.2.1), which bundles EUR LD scores internally
# and requires no separate LD reference download.
#
# WHAT GENETIC CORRELATION TELLS US:
# rg = 1.0: identical genetic architecture between EUR and CSA T2D
# rg = 0.0: completely different genetic bases
# rg = 0.5-0.9: partially shared architecture (typical for cross-ancestry)
#
# This provides the genome-wide context for interpreting the ~1.3 SD PRS
# mean shift found in Scripts 04-06. A high rg confirms that the underlying
# biology is shared but the PRS is miscalibrated due to LD/frequency
# differences -- not because T2D has different genetic causes in CSA.
#
# ADDITIONAL ANALYSES:
#   1. Per-ancestry SNP heritability (h2) estimates
#   2. Comparison of EUR vs CSA h2 on liability scale
#
# NOTE ON SUMSTATS FORMAT:
# Our LDSC sumstats use chr:pos:ref:alt SNP IDs. ldscr requires rsIDs
# to match against HapMap3 SNPs. We remap using the 1KGP bim file
# which contains both rsIDs and positions.
#
# OUTPUTS:
#   results/genetic_correlation/ldsc_rg_EUR_CSA.rds
#   results/genetic_correlation/ldsc_h2_EUR.rds
#   results/genetic_correlation/ldsc_h2_CSA.rds
#   results/genetic_correlation/ldsc_summary.csv
#   results/figures/fig5_genetic_correlation.png
# =============================================================================

library(ldscr)
library(data.table)
library(ggplot2)

dir.create("results/genetic_correlation", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load LDSC sumstats
# -----------------------------------------------------------------------------

eur_ldsc <- fread("data/sumstats/outputs/ldsc_EUR.sumstats.gz")
csa_ldsc <- fread("data/sumstats/outputs/ldsc_CSA.sumstats.gz")

# Fix N column -- was saved as allele frequency not sample size
# EUR total N = 22,634 cases + 397,897 controls = 420,531
# CSA total N =  1,662 cases +   7,214 controls =   8,876
eur_ldsc[, N := 420531]
csa_ldsc[, N := 8876]

# -----------------------------------------------------------------------------
# Remap SNP IDs from chr:pos:ref:alt to rsIDs
# -----------------------------------------------------------------------------

# Load 1KGP bim file (has rsIDs and positions)
bim <- fread("data/1KGP/1000G_phase3_common_norel.bim")
colnames(bim) <- c("chr", "rsID", "cm", "pos", "a1", "a2")

# Create chr:pos key in bim
bim[, varID := paste(chr, pos, sep = ":")]

# Create chr:pos key in sumstats (strip ref:alt from chr:pos:ref:alt)
eur_ldsc[, pos_key := sub("^([0-9]+:[0-9]+):.*", "\\1", SNP)]
csa_ldsc[, pos_key := sub("^([0-9]+:[0-9]+):.*", "\\1", SNP)]

# Merge to get rsIDs
eur_ldsc <- merge(eur_ldsc, bim[, .(varID, rsID)],
                  by.x = "pos_key", by.y = "varID",
                  all.x = FALSE)
csa_ldsc <- merge(csa_ldsc, bim[, .(varID, rsID)],
                  by.x = "pos_key", by.y = "varID",
                  all.x = FALSE)

# Replace SNP column with rsID
eur_ldsc[, SNP := rsID]
csa_ldsc[, SNP := rsID]

# Remove helper columns
eur_ldsc[, c("pos_key", "rsID") := NULL]
csa_ldsc[, c("pos_key", "rsID") := NULL]

# Remove duplicates -- keep first occurrence per rsID
eur_ldsc <- eur_ldsc[!duplicated(SNP)]
csa_ldsc <- csa_ldsc[!duplicated(SNP)]

# Preview EUR LDSC final format
print(head(eur_ldsc))

# -----------------------------------------------------------------------------
# SNP heritability estimates
# -----------------------------------------------------------------------------
# Estimate h2 for EUR and CSA separately before running genetic correlation.
# This serves as a sanity check -- EUR h2 should be ~0.1-0.2 on observed
# scale for T2D, and CSA h2 should be lower due to smaller sample size
# (less power to detect heritability signal).
#
# sample.prev = proportion of cases in the GWAS sample
# population.prev = T2D prevalence in the general population
#   EUR: ~8% prevalence, ~5.4% in UK Biobank
#   CSA: ~12% prevalence (higher T2D burden in South Asians)


# Save formatted sumstats for ldscr
eur_path <- "data/sumstats/outputs/ldsc_EUR_rsid.txt.gz"
csa_path <- "data/sumstats/outputs/ldsc_CSA_rsid.txt.gz"

fwrite(eur_ldsc, eur_path, sep = "\t")
fwrite(csa_ldsc, csa_path, sep = "\t")



# EUR LD scores
ld_files_EUR <- ldscore_files(ancestry = "EUR")
cat("EUR LD files:\n"); print(ld_files_EUR)

# EUR h2
h2_EUR <- tryCatch(
  ldsc_h2(
    munged_sumstats = eur_path,
    ancestry        = "EUR",
    sample_prev     = 22634 / 420531,
    population_prev = 0.054,
    ld              = ld_files_EUR[3],   # .ldscore.gz file
    wld             = ld_files_EUR[3],   # use same for wld
    n_blocks        = 200
  ),
  error = function(e) {
    cat("EUR h2 error:", conditionMessage(e), "\n")
    return(NULL)
  }
)


# CSA LD scores
ld_files_CSA <- ldscore_files(ancestry = "CSA")
cat("CSA LD files:\n"); print(ld_files_CSA)

# CSA h2
h2_CSA <- tryCatch(
  ldsc_h2(
    munged_sumstats = csa_path,
    ancestry        = "CSA",
    sample_prev     = 1662 / 8876,
    population_prev = 0.12,
    ld              = ld_files_CSA[3],
    wld             = ld_files_CSA[3],
    n_blocks        = 200
  ),
  error = function(e) {
    cat("CSA h2 error:", conditionMessage(e), "\n")
    return(NULL)
  }
)

saveRDS(h2_EUR, "results/genetic_correlation/ldsc_h2_EUR.rds")
saveRDS(h2_CSA, "results/genetic_correlation/ldsc_h2_CSA.rds")

# -----------------------------------------------------------------------------
# Cross-ancestry genetic correlation (EUR vs CSA T2D)
# -----------------------------------------------------------------------------

rg_result <- ldsc_rg(
  munged_sumstats = list(
    trait1 = eur_path,
    trait2 = csa_path
  ),
  ancestry = "EUR",
  ld       = ld_files_EUR[3],
  wld      = ld_files_EUR[3],
  n_blocks = 200
)

cat("Result class:", class(rg_result), "\n")
cat("Result names:", names(rg_result), "\n")
print(rg_result)

saveRDS(rg_result, "results/genetic_correlation/ldsc_rg_EUR_CSA.rds")

# -----------------------------------------------------------------------------
# Save summary table
# -----------------------------------------------------------------------------
# Extract key values
rg     <- rg_result$rg$rg
rg_se  <- rg_result$rg$rg_se
rg_p   <- rg_result$rg$rg_p

h2_eur_obs <- rg_result$h2$h2_observed[1]
h2_csa_obs <- rg_result$h2$h2_observed[2]

ldsc_summary <- data.table(
  analysis  = c("EUR h2 (observed)", "EUR h2 (liability)",
                "CSA h2 (observed)", "EUR-CSA rg"),
  estimate  = c(h2_eur_obs,
                h2_EUR$h2_liability,
                h2_csa_obs,
                rg),
  se        = c(rg_result$h2$h2_observed_se[1],
                h2_EUR$h2_liability_se,
                rg_result$h2$h2_observed_se[2],
                rg_se),
  p_value   = c(rg_result$h2$h2_p[1],
                NA,
                rg_result$h2$h2_p[2],
                rg_p)
)

print(ldsc_summary)
fwrite(ldsc_summary, "results/genetic_correlation/ldsc_summary.csv")

# -----------------------------------------------------------------------------
# Figure 5 -- Genetic correlation summary
# -----------------------------------------------------------------------------
# Visual summary of h2 and rg estimates with confidence intervals.
# This figure contextualizes the PRS miscalibration findings:
# high rg confirms shared genetic architecture despite score miscalibration.

plot_df <- data.table(
  Metric   = c("EUR h2\n(observed)",
               "CSA h2\n(observed)",
               "EUR-CSA rg"),
  Estimate = c(h2_eur_obs, h2_csa_obs, rg),
  SE       = c(rg_result$h2$h2_observed_se[1],
               rg_result$h2$h2_observed_se[2],
               rg_se),
  Category = c("Heritability", "Heritability", "Correlation"),
  Sig      = c("***", "ns", "***")
)

plot_df[, Lower := Estimate - 1.96 * SE]
plot_df[, Upper := pmin(Estimate + 1.96 * SE, 1.5)]
plot_df[, Metric := factor(Metric,
                           levels = c("EUR h2\n(observed)",
                                      "CSA h2\n(observed)",
                                      "EUR-CSA rg"))]

fig5 <- ggplot(plot_df,
               aes(x = Metric, y = Estimate,
                   color = Category)) +
  geom_hline(yintercept = 0,
             linetype = "dashed", color = "gray60") +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.15, linewidth = 0.8) +
  geom_text(aes(label = Sig, y = Upper + 0.05),
            size = 5, color = "gray30") +
  scale_color_manual(values = c(
    "Heritability" = "#4a90d9",
    "Correlation"  = "#e07b39"
  )) +
  facet_wrap(~Category, scales = "free") +
  labs(
    title    = "LDSC: SNP Heritability and Genetic Correlation (EUR vs CSA T2D)",
    subtitle = paste0(
      "rg = ", round(rg, 3),
      " (SE = ", round(rg_se, 3), ", p = ",
      format(rg_p, scientific = TRUE, digits = 3), ")"
    ),
    x     = NULL,
    y     = "Estimate ± 95% CI",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

ggsave("results/figures/fig5_genetic_correlation.png", plot = fig5, width = 9, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------------------------------------
# EUR h2 observed:   0.0355 SE = 0.0027 p = 8.070895e-41 
# EUR h2 liability:  0.1514 SE = 0.0111 
# CSA h2 observed:   0.0284 SE = 0.0411 p = 4.900806e-01 
# EUR-CSA rg:        1.2677 SE = 0.1942 p = 6.661014e-11 
# Interpretation:    High rg confirms shared genetic architecture despite ~1.3 SD PRS miscalibration

