# =============================================================================
# Script 06: PRS Evaluation
# Project:   Cross-ancestry T2D PRS Transferability (EUR -> CSA)
# =============================================================================
# OVERVIEW:
# This script performs the primary evaluation of PRS transferability across
# all scenarios from Scripts 04 (C+T) and 05 (LDpred2-inf). Since 1KGP
# lacks T2D phenotype data, evaluation uses:
#
#   1. SCORE DISTRIBUTION ANALYSIS
#      Density plots and summary statistics of PRS distributions in EUR vs
#      SAS under EUR-parameter standardization. Mean shift quantifies
#      miscalibration that would affect clinical deployment.
#
#   2. DECILE ANALYSIS
#      Proportion of SAS individuals in top/bottom PRS deciles defined by
#      EUR distribution. Clinically: what fraction of SAS individuals would
#      be flagged as high-risk under a EUR-derived score?
#
#   3. CROSS-ANCESTRY SCORE CORRELATION
#      Pearson correlation between EUR-derived and CSA-derived PRS within
#      SAS individuals. Low correlation means the two scores rank the same
#      people differently -- a clinically important finding.
#
#   4. METHOD COMPARISON (C+T vs LDpred2-inf)
#      Direct comparison of mean shifts across both methods to confirm
#      miscalibration is method-independent.
#
#   5. VARIANCE EXPLAINED BY PRS
#      Within each population, R² of PRS regressed on PC1 to assess whether
#      score variation reflects ancestry (bad) or genetic variation (good).
#
# OUTPUTS:
#   results/PRS/evaluation/prs_summary_stats.csv
#   results/PRS/evaluation/decile_analysis.csv
#   results/PRS/evaluation/score_correlation.csv
#   results/PRS/evaluation/method_comparison.csv
#   results/figures/fig1_prs_distributions.png
#   results/figures/fig2_decile_analysis.png
#   results/figures/fig3_score_correlation.png
#   results/figures/fig4_method_comparison.png
# =============================================================================

library(data.table)
library(ggplot2)
library(patchwork)

dir.create("results/PRS/evaluation", recursive = TRUE, showWarnings = FALSE)

# Consistent color scheme across all figures
pop_colors <- c("EUR" = "#4a90d9", "SAS" = "#e07b39")
legend_labels <- c("EUR" = "European (EUR)", "SAS" = "South Asian (SAS)")

# -----------------------------------------------------------------------------
# STEP 1: Load all PRS scores
# -----------------------------------------------------------------------------

cat("Loading PRS scores...\n\n")

# C+T scores
ct_EUR    <- fread("results/PRS/CT/prs_CT_EUR_in_EUR.txt")
ct_SAS    <- fread("results/PRS/CT/prs_CT_EUR_in_SAS.txt")
ct_SASld  <- fread("results/PRS/CT/prs_CT_EURLDsas_in_SAS.txt")
ct_CSA    <- fread("results/PRS/CT/prs_CT_CSA_in_SAS.txt")

# LDpred2-inf scores
ld_EUR    <- fread("results/PRS/LDpred2/prs_LDpred2_EUR_in_EUR.txt")
ld_SAS    <- fread("results/PRS/LDpred2/prs_LDpred2_EUR_in_SAS.txt")
ld_CSA    <- fread("results/PRS/LDpred2/prs_LDpred2_CSA_in_SAS.txt")

# PCA covariates (for ancestry regression analysis)
pca_EUR   <- fread("data/1KGP/pca_covariates_EUR.txt")
pca_SAS   <- fread("data/1KGP/pca_covariates_SAS.txt")

# Verify column names
cat("C+T columns:    ", colnames(ct_EUR),  "\n")
cat("LDpred2 columns:", colnames(ld_EUR),  "\n\n")

# -----------------------------------------------------------------------------
# STEP 2: Summary statistics
# -----------------------------------------------------------------------------
# Compute mean, SD, and mean shift for all scenarios and methods.
# These numbers go directly into manuscript Table 1.

cat("Computing summary statistics...\n\n")

compute_summary <- function(prs_eur, prs_sas, method, scenario) {
  data.table(
    method        = method,
    scenario      = scenario,
    EUR_mean      = mean(prs_eur, na.rm = TRUE),
    EUR_sd        = sd(prs_eur,   na.rm = TRUE),
    SAS_mean      = mean(prs_sas, na.rm = TRUE),
    SAS_sd        = sd(prs_sas,   na.rm = TRUE),
    mean_shift    = mean(prs_sas, na.rm = TRUE) -
      mean(prs_eur, na.rm = TRUE),
    n_EUR         = sum(!is.na(prs_eur)),
    n_SAS         = sum(!is.na(prs_sas))
  )
}

summary_stats <- rbind(
  compute_summary(ct_EUR$PRS_std_EUR,  ct_SAS$PRS_std_EUR,
                  "C+T",       "EUR_effects_EUR_LD"),
  compute_summary(ct_EUR$PRS_std_EUR,  ct_SASld$PRS_std_EUR,
                  "C+T",       "EUR_effects_SAS_LD"),
  compute_summary(ct_EUR$PRS_std_EUR,  ct_CSA$PRS_std_EUR,
                  "C+T",       "CSA_effects_SAS_LD"),
  compute_summary(ld_EUR$PRS_std_EUR,  ld_SAS$PRS_std_EUR,
                  "LDpred2-inf", "EUR_effects_EUR_LD"),
  compute_summary(ld_EUR$PRS_std_EUR,  ld_CSA$PRS_std_EUR,
                  "LDpred2-inf", "CSA_effects_SAS_LD")
)

cat("PRS Summary Statistics:\n")
print(summary_stats)
cat("\n")

fwrite(summary_stats, "results/PRS/evaluation/prs_summary_stats.csv")
cat("Summary stats saved.\n\n")

# -----------------------------------------------------------------------------
# STEP 3: Decile analysis
# -----------------------------------------------------------------------------
# Define risk deciles from EUR PRS distribution.
# Compute proportion of SAS individuals in each decile.
# Clinical relevance: a well-calibrated PRS should place ~10% of any
# population in each decile. Enrichment in top decile = over-classification
# of SAS as high risk; enrichment in bottom decile = under-classification.

cat("Running decile analysis...\n\n")

decile_analysis <- function(prs_eur, prs_sas, method, scenario) {
  
  # Define decile cutoffs from EUR distribution
  decile_cuts <- quantile(prs_eur, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
  
  # Assign EUR and SAS individuals to deciles
  eur_decile <- cut(prs_eur, breaks = decile_cuts, include.lowest = TRUE, labels = 1:10)
  sas_decile <- cut(prs_sas, breaks = decile_cuts, include.lowest = TRUE, labels = 1:10)
  
  # Proportion in each decile
  eur_prop <- as.numeric(table(eur_decile)) / length(prs_eur[!is.na(prs_eur)])
  sas_prop <- as.numeric(table(sas_decile)) / length(prs_sas[!is.na(prs_sas)])
  
  data.table(
    method   = method,
    scenario = scenario,
    decile   = 1:10,
    EUR_prop = eur_prop,
    SAS_prop = sas_prop,
    # Enrichment ratio: SAS/EUR proportion per decile
    # >1 = SAS over-represented, <1 = under-represented
    enrichment = sas_prop / eur_prop
  )
}

decile_results <- rbind(
  decile_analysis(ct_EUR$PRS_std_EUR, ct_SAS$PRS_std_EUR,
                  "C+T", "EUR_effects_EUR_LD"),
  decile_analysis(ld_EUR$PRS_std_EUR, ld_SAS$PRS_std_EUR,
                  "LDpred2-inf", "EUR_effects_EUR_LD")
)

cat("Decile analysis results:\n")
print(decile_results)
cat("\n")

# Key clinical metrics
for(method in c("C+T", "LDpred2-inf")) {
  d <- decile_results[decile_results$method == method & 
                        decile_results$scenario == "EUR_effects_EUR_LD", ]
  cat(method, ":\n")
  cat("  % SAS in top decile (decile 10):",
      round(d$SAS_prop[10] * 100, 1), "%",
      "(EUR:", round(d$EUR_prop[10] * 100, 1), "%)\n")
  cat("  % SAS in bottom decile (decile 1):",
      round(d$SAS_prop[1] * 100, 1), "%",
      "(EUR:", round(d$EUR_prop[1] * 100, 1), "%)\n")
  cat("  Top decile enrichment ratio:",
      round(d$enrichment[10], 2), "\n\n")
}

fwrite(decile_results, "results/PRS/evaluation/decile_analysis.csv")

# -----------------------------------------------------------------------------
# STEP 4: Cross-ancestry score correlation
# -----------------------------------------------------------------------------
# Pearson correlation between EUR-derived and CSA-derived PRS within SAS.
# High correlation: both scores rank SAS individuals similarly (consistent)
# Low correlation: scores identify different people as high risk (divergent)
# This is the key locus-level transferability metric available without
# phenotype data.

cat("Computing cross-ancestry score correlations...\n\n")

# Merge C+T EUR-derived and CSA-derived scores for SAS individuals
ct_sas_merged <- merge(
  ct_SAS[,  .(IID, PRS_ct_EUR  = PRS_std_EUR)],
  ct_CSA[,  .(IID, PRS_ct_CSA  = PRS_std_EUR)],
  by = "IID"
)

# Merge LDpred2 EUR-derived and CSA-derived scores
ld_sas_merged <- merge(
  ld_SAS[, .(IID, PRS_ld_EUR = PRS_std_EUR.V1)],
  ld_CSA[, .(IID, PRS_ld_CSA = PRS_std_EUR.V1)],
  by = "IID"
)

# Pearson correlations
ct_cor  <- cor(ct_sas_merged$PRS_ct_EUR, ct_sas_merged$PRS_ct_CSA,
               use = "complete.obs")
ld_cor  <- cor(ld_sas_merged$PRS_ld_EUR, ld_sas_merged$PRS_ld_CSA,
               use = "complete.obs")

cat("C+T EUR vs CSA PRS correlation in SAS:     ", round(ct_cor, 4), "\n")
cat("LDpred2 EUR vs CSA PRS correlation in SAS: ", round(ld_cor, 4), "\n\n")

cor_results <- data.table(
  method       = c("C+T", "LDpred2-inf"),
  cor_EUR_CSA  = c(ct_cor, ld_cor),
  n_SAS        = c(nrow(ct_sas_merged), nrow(ld_sas_merged))
)

fwrite(cor_results, "results/PRS/evaluation/score_correlation.csv")

# -----------------------------------------------------------------------------
# STEP 5: Method comparison table
# -----------------------------------------------------------------------------

method_comparison <- summary_stats[, .(
  method, scenario, mean_shift, EUR_sd, SAS_sd
)]
fwrite(method_comparison, "results/PRS/evaluation/method_comparison.csv")

# -----------------------------------------------------------------------------
# STEP 6: Figure 1 -- PRS distributions (main manuscript figure)
# -----------------------------------------------------------------------------
# Four-panel figure showing PRS distributions across all scenarios.
# Panel A: C+T EUR effects EUR LD
# Panel B: C+T EUR effects SAS LD
# Panel C: C+T CSA effects SAS LD
# Panel D: LDpred2-inf EUR effects EUR LD

cat("Generating Figure 1: PRS distributions...\n")

make_density_plot <- function(prs_eur, prs_sas, title_str) {
  
  sas_mean <- round(mean(prs_sas, na.rm = TRUE), 3)
  
  df <- data.table(
    PRS = c(prs_eur, prs_sas),
    Pop = c(rep("EUR", length(prs_eur)),
            rep("SAS", length(prs_sas)))
  )
  
  ggplot(df, aes(x = PRS, fill = Pop, color = Pop)) +
    geom_density(alpha = 0.4, linewidth = 0.7) +
    scale_fill_manual(values = pop_colors, labels = legend_labels) +
    scale_color_manual(values = pop_colors, labels = legend_labels) +
    geom_vline(xintercept = 0,
               linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_vline(xintercept = mean(prs_sas, na.rm = TRUE),
               linetype = "dashed", color = "#e07b39", linewidth = 0.7) +
    annotate("text",
             x = mean(prs_sas, na.rm = TRUE) + 0.05,
             y = Inf, vjust = 1.5, hjust = 0, size = 2.8,
             color = "#e07b39",
             label = paste0("SAS mean = ", sas_mean)) +
    labs(title = title_str,
         x = "Standardized PRS",
         y = "Density",
         fill  = NULL,
         color = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title       = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank(),
          legend.position  = "bottom")
}

p1a <- make_density_plot(
  ct_EUR$PRS_std_EUR, ct_SAS$PRS_std_EUR,
  "A. C+T: EUR effects, EUR LD"
)

p1b <- make_density_plot(
  ct_EUR$PRS_std_EUR, ct_SASld$PRS_std_EUR,
  "B. C+T: EUR effects, SAS LD"
)

p1c <- make_density_plot(
  ct_EUR$PRS_std_EUR, ct_CSA$PRS_std_EUR,
  "C. C+T: CSA effects, SAS LD"
)

p1d <- make_density_plot(
  ld_EUR$PRS_std_EUR, ld_SAS$PRS_std_EUR,
  "D. LDpred2-inf: EUR effects, EUR LD"
)

# Add shared legend
legend_plot <- ggplot(
  data.table(PRS = c(0,0), Pop = c("EUR","SAS")),
  aes(x = PRS, fill = Pop, color = Pop)
) +
  geom_density() +
  scale_fill_manual(values  = pop_colors,
                    labels  = c("EUR" = "European (EUR)",
                                "SAS" = "South Asian (SAS)")) +
  scale_color_manual(values = pop_colors) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title    = element_blank())

fig2 <- (p1a | p1b | p1c | p1d) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave("results/figures/fig2_prs_distributions2.png",
       plot = fig2, width = 18, height = 4.5, dpi = 300)
cat("Saved: results/figures/fig1_prs_distributions.png\n")

# -----------------------------------------------------------------------------
# STEP 7: Figure 4 -- Decile analysis
# -----------------------------------------------------------------------------
# Bar chart showing proportion of SAS individuals in each EUR-defined
# decile. Red dashed line at 10% = expected under perfect calibration.

cat("Generating Figure 4: Decile analysis...\n")

decile_plot_df <- melt(
  decile_results[scenario == "EUR_effects_EUR_LD"],
  id.vars       = c("method", "scenario", "decile"),
  measure.vars  = c("EUR_prop", "SAS_prop"),
  variable.name = "Population",
  value.name    = "Proportion"
)
decile_plot_df[, Population := ifelse(Population == "EUR_prop", "EUR", "SAS")]
decile_plot_df[, Proportion_pct := Proportion * 100]

fig4 <- ggplot(decile_plot_df,
               aes(x = factor(decile), y = Proportion_pct,
                   fill = Population)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_hline(yintercept = 10, linetype = "dashed",
             color = "gray40", linewidth = 0.7) +
  scale_fill_manual(values  = pop_colors,
                    labels  = c("EUR" = "European (EUR)",
                                "SAS" = "South Asian (SAS)")) +
  facet_wrap(~method, ncol = 2) +
  labs(
    x        = "PRS Decile",
    y        = "% of Individuals",
    fill     = "Population"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

ggsave("results/figures/fig4_decile_analysis.png",
       plot = fig4, width = 14, height = 4, dpi = 300)
cat("Saved: results/figures/fig2_decile_analysis.png\n")

# -----------------------------------------------------------------------------
# STEP 8: Figure 5 -- Score correlation
# -----------------------------------------------------------------------------
# Scatter plot of EUR-derived vs CSA-derived PRS within SAS individuals.
# Each point is one SAS individual. Diagonal = perfect agreement.

cat("Generating Figure 3: Score correlation...\n")

p5a <- ggplot(ct_sas_merged,
              aes(x = PRS_ct_EUR, y = PRS_ct_CSA)) +
  geom_point(alpha = 0.4, size = 1.5, color = "#e07b39") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray40") +
  geom_smooth(method = "lm", se = TRUE,
              color = "#4a90d9", linewidth = 0.8) +
  annotate("text", x = -Inf, y = Inf,
           hjust = -0.1, vjust = 1.5, size = 3.5,
           label = paste0("r = ", round(ct_cor, 3))) +
  labs(title    = "A. C+T",
       x        = "EUR-derived PRS (EUR effects)",
       y        = "CSA-derived PRS (CSA effects)") +
  theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

p5b <- ggplot(ld_sas_merged,
              aes(x = PRS_ld_EUR, y = PRS_ld_CSA)) +
  geom_point(alpha = 0.4, size = 1.5, color = "#e07b39") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray40") +
  geom_smooth(method = "lm", se = TRUE,
              color = "#4a90d9", linewidth = 0.8) +
  annotate("text", x = -Inf, y = Inf,
           hjust = -0.1, vjust = 1.5, size = 3.5,
           label = paste0("r = ", round(ld_cor, 3))) +
  labs(title    = "B. LDpred2-inf",
       x        = "EUR-derived PRS (EUR effects)",
       y        = "CSA-derived PRS (CSA effects)") +
  theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

fig5 <- p5a | p5b
fig5 <- fig5 + plot_annotation(
  title    = "EUR vs CSA PRS Correlation in SAS Individuals",
  subtitle = paste0(
    "Each point = one SAS individual (n=", nrow(ct_sas_merged), ")",
    "\nDiagonal = perfect agreement | Blue line = linear fit"
  )
)

ggsave("results/figures/fig3_score_correlation.png",
       plot = fig3, width = 10, height = 5, dpi = 300)
cat("Saved: results/figures/fig3_score_correlation.png\n")

# -----------------------------------------------------------------------------
# STEP 9: Method comparison
# -----------------------------------------------------------------------------
# Summary bar chart comparing mean shifts across all scenarios and methods.
# This is the key synthesis figure showing that miscalibration is
# consistent regardless of PRS method or effect size ancestry.

cat("Generating Figure 4: Method comparison...\n")

# Add readable labels
method_comparison[, scenario_label := fcase(
  scenario == "EUR_effects_EUR_LD", "EUR effects\nEUR LD",
  scenario == "EUR_effects_SAS_LD", "EUR effects\nSAS LD",
  scenario == "CSA_effects_SAS_LD", "CSA effects\nSAS LD"
)]

fig4 <- ggplot(method_comparison,
               aes(x = scenario_label, y = mean_shift,
                   fill = method)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.6) +
  geom_hline(yintercept = 0,
             linetype = "solid", color = "gray40", linewidth = 0.5) +
  geom_hline(yintercept = 1.3165,
             linetype = "dashed", color = "gray60", linewidth = 0.5) +
  scale_fill_manual(values = c("C+T" = "#4a90d9",
                               "LDpred2-inf" = "#2c5f8a")) +
  annotate("text", x = Inf, y = 1.3165,
           hjust = 1.1, vjust = -0.5, size = 3,
           color = "gray50",
           label = "C+T EUR/EUR benchmark") +
  labs(
    title    = "PRS Mean Shift: EUR vs SAS Across Methods and Scenarios",
    subtitle = "Mean shift = SAS mean PRS - EUR mean PRS (EUR-standardized units)",
    x        = "PRS Scenario",
    y        = "Mean Shift (SDs)",
    fill     = "Method"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

ggsave("results/figures/fig4_method_comparison.png",
       plot = fig4, width = 9, height = 5.5, dpi = 300)

# -----------------------------------------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------------------------------------
# PRS Summary Statistics:
print(summary_stats[, .(method, scenario, mean_shift, n_EUR, n_SAS)])

#Score Correlations (EUR vs CSA in SAS):
print(cor_results)

# Decile Analysis (top decile % in SAS):
for(method in c("C+T", "LDpred2-inf")) {
  d <- decile_results[decile_results$method == method &
                        decile_results$scenario == "EUR_effects_EUR_LD", ]
  cat(method, "top decile SAS %:",
      round(d$SAS_prop[10] * 100, 1), "%",
      "vs EUR 10%\n")
}

