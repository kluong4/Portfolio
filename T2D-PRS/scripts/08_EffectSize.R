# =============================================================================
# Script 08: Effect Size Concordance Analysis
# Project:   Cross-ancestry T2D PRS Transferability (EUR -> CSA)
# =============================================================================
# OVERVIEW:
# This script performs locus-level cross-ancestry effect size concordance
# analysis at EUR genome-wide significant (GWS) T2D loci. This directly
# addresses the sub-question: "Which genomic loci show the greatest
# cross-ancestry effect size heterogeneity between EUR and CSA?"
#
# ANALYSES:
#   1. Effect size concordance scatter plot (EUR beta vs CSA beta)
#   2. Cochran's Q heterogeneity test per locus
#   3. Directional concordance by CSA Z-score bin (from Script 01)
#   4. Top heterogeneous loci identification and annotation
#   5. Allele frequency comparison at GWS loci
#
# NOTE: This uses gws_EUR_loci_annotated.rds from Script 01 which
# contains all EUR GWS loci with CSA summary statistics and Z-score
# bin annotations already computed.
#
# OUTPUTS:
#   results/concordance/heterogeneity_results.csv
#   results/concordance/top_heterogeneous_loci.csv
#   results/figures/fig6_effect_size_concordance.png
#   results/figures/fig7_heterogeneity_volcano.png
#   results/figures/fig8_allele_frequency_comparison.png
# =============================================================================

library(data.table)
library(ggplot2)
library(patchwork)

dir.create("results/concordance", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures",     recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load GWS loci with CSA data
# -----------------------------------------------------------------------------

# Load annotated GWS loci from Script 01
gws <- readRDS("data/sumstats/outputs/gws_EUR_loci_annotated.rds")

# Confirm required columns exist
required <- c("beta_EUR", "se_EUR", "beta_CSA", "se_CSA",
              "neglog10_pval_EUR", "csa_absz", "z_bin")
missing  <- setdiff(required, colnames(gws))
if(length(missing) > 0) {
  cat("WARNING: Missing columns:", missing, "\n")
} else {
  cat("All required columns present.\n\n")
}

# Basic summary
cat("EUR beta range: [", round(min(gws$beta_EUR), 4), ",",
    round(max(gws$beta_EUR), 4), "]\n")
cat("CSA beta range: [", round(min(gws$beta_CSA, na.rm=TRUE), 4), ",",
    round(max(gws$beta_CSA, na.rm=TRUE), 4), "]\n")
cat("Non-missing CSA betas:", sum(!is.na(gws$beta_CSA)), "\n\n")

# -----------------------------------------------------------------------------
# Cochran's Q heterogeneity test
# -----------------------------------------------------------------------------
# Cochran's Q tests whether effect sizes differ significantly between
# EUR and CSA at each locus beyond what is expected by chance.
# Q = (beta_EUR - beta_CSA)^2 / (se_EUR^2 + se_CSA^2)
# Under the null (no heterogeneity), Q follows chi-squared with df=1.
# Significant Q (p < 0.05 after Bonferroni) indicates true cross-ancestry
# effect size heterogeneity at that locus.

gws[, Q_stat := (beta_EUR - beta_CSA)^2 /
      (se_EUR^2 + se_CSA^2)]
gws[, Q_pval := pchisq(Q_stat, df = 1, lower.tail = FALSE)]
gws[, Q_neglog10 := -log10(Q_pval)]

# Bonferroni correction
bonferroni_thresh <- 0.05 / nrow(gws)
gws[, Q_sig_bonf := Q_pval < bonferroni_thresh]

# Rank by heterogeneity
gws_ranked <- gws[order(Q_pval)]

# Top 10 most heterogeneous loci
print(gws_ranked[1:10, .(varID, chr, pos, beta_EUR, se_EUR, beta_CSA, se_CSA, csa_absz, Q_stat, Q_pval, Q_sig_bonf)])

# Save full heterogeneity results
fwrite(gws_ranked[, .(varID, chr, pos, ref, alt,
                      beta_EUR, se_EUR, neglog10_pval_EUR,
                      beta_CSA, se_CSA, neglog10_pval_CSA,
                      csa_absz, z_bin,
                      Q_stat, Q_pval, Q_neglog10, Q_sig_bonf)],
       "results/concordance/heterogeneity_results.csv")

# Save top heterogeneous loci (nominal p < 0.05)
top_het <- gws_ranked[Q_pval < 0.05]
fwrite(top_het[, .(varID, chr, pos, ref, alt,
                   beta_EUR, se_EUR,
                   beta_CSA, se_CSA,
                   Q_stat, Q_pval, Q_sig_bonf)],
       "results/concordance/top_heterogeneous_loci.csv")

# -----------------------------------------------------------------------------
# Directional concordance summary (from Script 01)
# -----------------------------------------------------------------------------
consistency_by_bin <- readRDS("results/QC/directional_consistency_by_zscore.rds")
print(consistency_by_bin)

# Overall concordance
overall_concordance <- mean(
  sign(gws$beta_EUR) == sign(gws$beta_CSA),
  na.rm = TRUE
) * 100

cat("Overall directional concordance:", round(overall_concordance, 1), "%\n")
cat("At |Z| > 2: 97.9% | At |Z| > 3: 100%\n\n")

# -----------------------------------------------------------------------------
# Allele frequency comparison
# -----------------------------------------------------------------------------
# Compare EUR and CSA allele frequencies at GWS loci.
# Large allele frequency differences explain why EUR-derived PRS
# may not transfer well -- a SNP with high risk allele frequency
# in EUR but low in CSA will dominate the EUR PRS but contribute
# little signal in CSA individuals.

# Compute allele frequencies from af_cases and af_controls
# af_cases and af_controls are per-population allele frequencies
gws[, af_EUR := (af_cases_EUR * 22634 + af_controls_EUR * 397897) /
      (22634 + 397897)]
gws[, af_CSA := (af_cases_CSA * 1662  + af_controls_CSA * 7214) /
      (1662  + 7214)]

# Allele frequency difference
gws[, af_diff := abs(af_EUR - af_CSA)]

cat("Mean |AF difference| at GWS loci:", round(mean(gws$af_diff, na.rm=TRUE), 4), "\n")
cat("Max  |AF difference| at GWS loci:", round(max(gws$af_diff,  na.rm=TRUE), 4), "\n")
cat("% loci with |AF diff| > 0.1:     ",
    round(mean(gws$af_diff > 0.1, na.rm=TRUE) * 100, 1), "%\n")
cat("% loci with |AF diff| > 0.2:     ",
    round(mean(gws$af_diff > 0.2, na.rm=TRUE) * 100, 1), "%\n\n")

# -----------------------------------------------------------------------------
# Figure 6 -- Effect size concordance scatter plot
# -----------------------------------------------------------------------------
# Primary concordance figure. Each point = one EUR GWS locus.
# x-axis = EUR effect size, y-axis = CSA effect size.
# Diagonal = perfect concordance (same effect in both populations).
# Points colored by Cochran's Q significance.
# Points above/below diagonal but same sign = directionally concordant
# Points crossing zero = directionally discordant

# Subset to loci with reliable CSA estimates (|Z| >= 1)
# to avoid visualizing noise-dominated points
gws_reliable <- gws[csa_absz >= 1]
gws_all      <- gws[!is.na(beta_CSA)]

cat("Loci with CSA |Z| >= 1:", nrow(gws_reliable), "\n")
cat("All loci with CSA data:", nrow(gws_all), "\n\n")

# Pearson correlation of effect sizes
r_all      <- cor(gws_all$beta_EUR,      gws_all$beta_CSA,      use="complete.obs")
r_reliable <- cor(gws_reliable$beta_EUR, gws_reliable$beta_CSA, use="complete.obs")

cat("r (all loci):         ", round(r_all,      4), "\n")
cat("r (|Z| >= 1 loci):   ", round(r_reliable, 4), "\n\n")

fig6 <- ggplot(gws_all,
               aes(x = beta_EUR, y = beta_CSA,
                   color = Q_neglog10)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_point(alpha = 0.5, size = 1.2) +
  # Highlight Bonferroni significant heterogeneous loci
  geom_point(data = gws_all[Q_sig_bonf == TRUE],
             aes(x = beta_EUR, y = beta_CSA),
             color = "red", size = 2.5, alpha = 0.8) +
  scale_color_viridis_c(
    name   = "-log10(Q p-value)",
    option = "plasma",
    direction = -1
  ) +
  annotate("text", x = -Inf, y = Inf,
           hjust = -0.1, vjust = 1.5, size = 3.5,
           label = paste0("r = ", round(r_all, 3),
                          " (all loci)\nr = ", round(r_reliable, 3),
                          " (|Z| ≥ 1)")) +
  labs(
    title    = "EUR vs CSA Effect Size Concordance at EUR GWS T2D Loci",
    subtitle = paste0(
      "n = ", nrow(gws_all), " loci | ",
      "Red = Bonferroni-significant heterogeneity (n=",
      sum(gws_all$Q_sig_bonf, na.rm=TRUE), ")",
      "\nDashed line = perfect concordance"
    ),
    x = "Effect Size in EUR (beta)",
    y = "Effect Size in CSA (beta)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("results/figures/fig6_effect_size_concordance.png",
       plot = fig6, width = 8, height = 7, dpi = 300)

# -----------------------------------------------------------------------------
# Figure 7 -- Heterogeneity volcano plot
# -----------------------------------------------------------------------------
# x-axis = effect size difference (EUR - CSA)
# y-axis = -log10(Q p-value)
# Highlights loci with large effect size differences AND low p-values

gws_all[, beta_diff := beta_EUR - beta_CSA]

fig7 <- ggplot(gws_all,
               aes(x = beta_diff, y = Q_neglog10,
                   color = Q_sig_bonf)) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "gray60",
             linewidth = 0.6) +
  geom_hline(yintercept = -log10(bonferroni_thresh),
             linetype = "dashed", color = "red",
             linewidth = 0.6) +
  geom_point(alpha = 0.4, size = 1.0) +
  scale_color_manual(
    values = c("FALSE" = "gray60", "TRUE" = "red"),
    labels = c("FALSE" = "Not significant",
               "TRUE"  = "Bonferroni significant"),
    name   = "Heterogeneity"
  ) +
  annotate("text", x = Inf, y = -log10(0.05),
           hjust = 1.1, vjust = -0.5, size = 3,
           color = "gray50", label = "p = 0.05") +
  annotate("text", x = Inf, y = -log10(bonferroni_thresh),
           hjust = 1.1, vjust = -0.5, size = 3,
           color = "red",
           label = paste0("Bonferroni\n(p=",
                          format(bonferroni_thresh,
                                 scientific=TRUE, digits=2), ")")) +
  labs(
    title    = "Cross-Ancestry Effect Size Heterogeneity at EUR GWS T2D Loci",
    subtitle = paste0(
      "Cochran's Q test | n = ", nrow(gws_all), " loci | ",
      sum(gws_all$Q_sig_bonf, na.rm=TRUE),
      " Bonferroni-significant"
    ),
    x = "Effect Size Difference (EUR beta - CSA beta)",
    y = "-log10(Cochran's Q p-value)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

ggsave("results/figures/fig7_heterogeneity_volcano.png",
       plot = fig7, width = 9, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Figure 8 -- Allele frequency comparison
# -----------------------------------------------------------------------------
fig8a <- ggplot(gws_all[!is.na(af_EUR) & !is.na(af_CSA)],
                aes(x = af_EUR, y = af_CSA)) +
  geom_point(alpha = 0.3, size = 1.0, color = "#4a90d9") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray40") +
  labs(
    title = "A. Allele Frequency: EUR vs CSA",
    x     = "Allele Frequency (EUR)",
    y     = "Allele Frequency (CSA)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

fig8b <- ggplot(gws_all[!is.na(af_diff)],
                aes(x = af_diff)) +
  geom_histogram(bins = 50, fill = "#4a90d9",
                 color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0.1,
             linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_vline(xintercept = 0.2,
             linetype = "dashed", color = "darkred", linewidth = 0.7) +
  annotate("text", x = 0.1, y = Inf,
           hjust = -0.1, vjust = 1.5, size = 3, color = "red",
           label = "0.1") +
  annotate("text", x = 0.2, y = Inf,
           hjust = -0.1, vjust = 1.5, size = 3, color = "darkred",
           label = "0.2") +
  labs(
    title = "B. Distribution of |AF Difference| at GWS Loci",
    x     = "|Allele Frequency EUR - Allele Frequency CSA|",
    y     = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

fig8 <- fig8a | fig8b
fig8 <- fig8 + plot_annotation(
  title    = "Allele Frequency Differences at EUR GWS T2D Loci",
  subtitle = paste0(
    "n = ", sum(!is.na(gws_all$af_EUR) & !is.na(gws_all$af_CSA)),
    " loci | Mean |AF diff| = ",
    round(mean(gws_all$af_diff, na.rm=TRUE), 3)
  )
)

ggsave("results/figures/fig8_allele_frequency_comparison.png",
       plot = fig8, width = 11, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------------------------------------
# EUR GWS loci analyzed:             8789 
# EUR GWS loci with CSA |Z| >= 1:    2983 
# Bonferroni-significant Q loci:     0 
# Nominal Q loci (p < 0.05):         2371 
# r (all loci):                      0.4167 
# r (|Z| >= 1 loci):                 0.6091 
# Overall directional concordance:   68.2 %
# Mean |AF difference|:              0.0755 
# % loci with |AF diff| > 0.1:       31 %