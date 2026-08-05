# =============================================================================
# Script 10: Final Publication-Ready Figures
# Project:   Cross-ancestry T2D PRS Transferability (EUR -> CSA)
# =============================================================================
# OVERVIEW:
# This script assembles all results into final publication-ready figures
# for three deliverables:
#   1. MANUSCRIPT -- high resolution, detailed, with statistics
#   2. CONFERENCE TALK -- simplified, large fonts, key message visible
#   3. POSTER -- compact multi-panel, self-contained
#
# FIGURE INVENTORY:
#   Fig 1: PRS distribution comparison (C+T EUR effects) -- MAIN RESULT
#   Fig 2: Decile analysis -- CLINICAL IMPACT
#   Fig 3: Method comparison bar chart -- ROBUSTNESS
#   Fig 4: Directional concordance by Z-score bin -- QC FINDING
#   Fig 5: Effect size concordance scatter -- LOCUS LEVEL
#   Fig 6: Genetic correlation (LDSC) -- ARCHITECTURE
#   Fig 7: GO pathway enrichment -- BIOLOGY
#   Fig 8: Poster composite -- ALL KEY RESULTS
# =============================================================================

library(data.table)
library(ggplot2)
library(patchwork)

# Consistent color scheme throughout
pop_colors <- c("EUR" = "#4a90d9", "SAS" = "#e07b39")
method_colors <- c("C+T" = "#4a90d9", "LDpred2-inf" = "#2c5f8a")


# -----------------------------------------------------------------------------
# Figure 1
# -----------------------------------------------------------------------------

pca_data <- readRDS("data/1KGP/pca_combined.rds")
pca <- pca_data$pca
scores_df_clean <- pca_data$scores

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

# -----------------------------------------------------------------------------
# STEP 1: Load all results
# -----------------------------------------------------------------------------

cat("Loading all results...\n\n")

# PRS scores
ct_EUR   <- fread("results/PRS/CT/prs_CT_EUR_in_EUR.txt")
ct_SAS   <- fread("results/PRS/CT/prs_CT_EUR_in_SAS.txt")
ct_SASld <- fread("results/PRS/CT/prs_CT_EURLDsas_in_SAS.txt")
ct_CSA   <- fread("results/PRS/CT/prs_CT_CSA_in_SAS.txt")
ld_EUR   <- fread("results/PRS/LDpred2/prs_LDpred2_EUR_in_EUR.txt")
ld_SAS   <- fread("results/PRS/LDpred2/prs_LDpred2_EUR_in_SAS.txt")
ld_CSA   <- fread("results/PRS/LDpred2/prs_LDpred2_CSA_in_SAS.txt")

# Evaluation results
summary_stats  <- fread("results/PRS/evaluation/prs_summary_stats.csv")
decile_results <- fread("results/PRS/evaluation/decile_analysis.csv")

# Directional concordance from Script 01
consistency_by_bin <- readRDS("results/QC/directional_consistency_by_zscore.rds")

# Effect size concordance from Script 08
het <- fread("results/concordance/heterogeneity_results.csv")

# Genetic correlation from Script 07
rg_result <- readRDS("results/genetic_correlation/ldsc_rg_EUR_CSA.rds")
h2_EUR    <- readRDS("results/genetic_correlation/ldsc_h2_EUR.rds")

# GO enrichment from Script 09
go_all <- fread("results/annotation/go_enrichment_all.csv")

cat("All results loaded.\n\n")

# LDpred2 standardization params
ldpred2_params <- readRDS(
  "results/PRS/LDpred2/ldpred2_eur_standardization_params.rds"
)
ldpred2_eur_mean <- ldpred2_params$eur_mean
ldpred2_eur_sd   <- ldpred2_params$eur_sd

# Restandardize LDpred2 scores
prs_EUR_ld <- (ld_EUR$PRS_raw - ldpred2_eur_mean) / ldpred2_eur_sd
prs_SAS_ld <- (ld_SAS$PRS_raw - ldpred2_eur_mean) / ldpred2_eur_sd

# -----------------------------------------------------------------------------
# STEP 2: MANUSCRIPT Figure 1 -- Main PRS Distribution Result
# -----------------------------------------------------------------------------
# Two-panel figure: C+T and LDpred2-inf side by side
# This is the primary result figure for the manuscript

cat("Generating Manuscript Figure 1...\n")

make_dist_panel <- function(prs_eur, prs_sas, title_str,
                            method_str, n_eur, n_sas) {
  
  sas_mean <- round(mean(prs_sas, na.rm=TRUE), 2)
  shift    <- round(mean(prs_sas, na.rm=TRUE) -
                      mean(prs_eur, na.rm=TRUE), 2)
  
  df <- data.table(
    PRS = c(prs_eur, prs_sas),
    Pop = c(rep("EUR", length(prs_eur)),
            rep("SAS", length(prs_sas)))
  )
  
  ggplot(df, aes(x = PRS, fill = Pop, color = Pop)) +
    geom_density(alpha = 0.4, linewidth = 0.8) +
    scale_fill_manual(values  = pop_colors,
                      labels  = c("EUR" = "European (EUR; n=503)",
                                  "SAS" = "South Asian (SAS; n=478)")) +
    scale_color_manual(values = pop_colors) +
    geom_vline(xintercept = 0,
               linetype = "dashed", color = "gray50",
               linewidth = 0.6) +
    geom_vline(xintercept = mean(prs_sas, na.rm=TRUE),
               linetype = "dashed", color = "#e07b39",
               linewidth = 0.8) +
    annotate("text",
             x = mean(prs_sas, na.rm=TRUE) + 0.08,
             y = Inf, vjust = 1.8, hjust = 0, size = 3.5,
             color = "#e07b39",
             label = paste0("SAS mean = +", sas_mean, " SDs")) +
    labs(title    = title_str,
         subtitle = paste0("Mean shift = +", shift, " SDs"),
         x        = "Standardized PRS (EUR mean=0, SD=1)",
         y        = "Density",
         fill     = "Population",
         color    = "Population") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(color = "#e07b39", size = 10),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_blank()
    )
}

p_ct <- make_dist_panel(
  ct_EUR$PRS_std_EUR, ct_SAS$PRS_std_EUR,
  "A. Clumping & Thresholding (C+T)",
  "C+T", nrow(ct_EUR), nrow(ct_SAS)
)

p_ld <- make_dist_panel(
  prs_EUR_ld, prs_SAS_ld,
  "B. LDpred2-inf",
  "LDpred2-inf", nrow(ld_EUR), nrow(ld_SAS)
)

ms_fig1 <- p_ct | p_ld
ms_fig1 <- ms_fig1 + patchwork::plot_annotation(
  title    = "Type 2 Diabetes PRS Miscalibration in South Asian Individuals",
  subtitle = paste0(
    "EUR-derived PRS standardized using EUR parameters (mean=0, SD=1) ",
    "and applied to SAS individuals\n",
    "Both methods independently confirm ~1.3 SD systematic over-scoring ",
    "of SAS individuals"
  ),
  caption  = paste0(
    "Pan-UK Biobank E11 summary statistics (EUR n_cases=22,634; ",
    "CSA n_cases=1,662) | 1000 Genomes Phase 3"
  )
)

ggsave("results/figures/manuscript/ms_fig1_prs_distribution.png",
       plot = ms_fig1, width = 12, height = 6, dpi = 300)
cat("Saved: manuscript/ms_fig1_prs_distribution.png\n")

# -----------------------------------------------------------------------------
# STEP 3: MANUSCRIPT Figure 2 -- Decile Analysis
# -----------------------------------------------------------------------------

cat("Generating Manuscript Figure 2...\n")

decile_dt <- as.data.table(decile_results)

decile_plot_df <- melt(
  decile_dt[scenario == "EUR_effects_EUR_LD"],
  id.vars       = c("method", "decile"),
  measure.vars  = c("EUR_prop", "SAS_prop"),
  variable.name = "Population",
  value.name    = "Proportion"
)
decile_plot_df[, Population    := ifelse(Population == "EUR_prop",
                                         "EUR", "SAS")]
decile_plot_df[, Proportion_pct := Proportion * 100]
decile_plot_df[, method := factor(method,
                                  levels = c("C+T", "LDpred2-inf"))]

ms_fig2 <- ggplot(decile_plot_df,
                  aes(x = factor(decile), y = Proportion_pct,
                      fill = Population)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
  geom_hline(yintercept = 10,
             linetype = "dashed", color = "gray40",
             linewidth = 0.7) +
  scale_fill_manual(values  = pop_colors,
                    labels  = c("EUR" = "European (EUR)",
                                "SAS" = "South Asian (SAS)")) +
  facet_wrap(~method, ncol = 2) +
  annotate("text", x = 9.5, y = 11.5, hjust = 1,
           size = 3, color = "gray40",
           label = "Expected (10%)") +
  labs(
    title    = "PRS Risk Decile Distribution: EUR vs SAS",
    subtitle = paste0(
      "Deciles defined by EUR PRS distribution | ",
      "~40% of SAS individuals fall in top decile vs 10% expected"
    ),
    x        = "PRS Decile (1 = lowest risk, 10 = highest risk)",
    y        = "% of Individuals",
    fill     = "Population",
    caption  = "EUR-derived PRS applied universally — clinical deployment simulation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray40", size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

ggsave("results/figures/manuscript/ms_fig2_decile_analysis.png",
       plot = ms_fig2, width = 11, height = 5.5, dpi = 300)
cat("Saved: manuscript/ms_fig2_decile_analysis.png\n")

# -----------------------------------------------------------------------------
# STEP 4: MANUSCRIPT Figure 3 -- Directional Concordance
# -----------------------------------------------------------------------------

cat("Generating Manuscript Figure 3...\n")

consist_dt <- as.data.table(consistency_by_bin)
consist_dt[, lower_ci := (same_dir/n -
                            1.96*sqrt((same_dir/n*(1-same_dir/n))/n))*100]
consist_dt[, upper_ci := pmin((same_dir/n +
                                 1.96*sqrt((same_dir/n*(1-same_dir/n))/n))*100, 100)]
consist_dt[, pct := same_dir/n*100]

ms_fig3 <- ggplot(consist_dt,
                  aes(x = z_bin, y = pct, group = 1)) +
  geom_hline(yintercept = 50,
             linetype = "dashed", color = "gray60",
             linewidth = 0.6) +
  geom_hline(yintercept = 100,
             linetype = "dotted", color = "gray40",
             linewidth = 0.6) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              fill = "#4a90d9", alpha = 0.2) +
  geom_line(color = "#4a90d9", linewidth = 1.2) +
  geom_point(color = "#4a90d9", size = 4) +
  geom_text(aes(label = paste0(round(pct, 1), "%\n(n=", n, ")")),
            vjust = -1.0, size = 3.2, color = "gray30") +
  scale_y_continuous(limits = c(40, 110),
                     breaks = c(50, 60, 70, 80, 90, 100),
                     labels = paste0(c(50,60,70,80,90,100), "%")) +
  labs(
    title    = "Effect Size Directional Concordance: EUR vs CSA by CSA Statistical Power",
    subtitle = paste0(
      "At EUR GWS T2D loci (n=8,789) | ",
      "Concordance increases from ~64% (noise) to 100% (strong signal)\n",
      "Low overall concordance (68%) explained entirely by statistical noise, ",
      "not true heterogeneity"
    ),
    x        = "CSA |Z-score| Bin (proxy for statistical power)",
    y        = "% Same Effect Direction (EUR vs CSA)",
    caption  = "Shaded band = 95% CI | Dashed = 50% (chance level)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray40", size = 9),
    panel.grid.minor = element_blank()
  )

ggsave("results/figures/manuscript/ms_fig3_directional_concordance.png",
       plot = ms_fig3, width = 9, height = 6, dpi = 300)
cat("Saved: manuscript/ms_fig3_directional_concordance.png\n")

# -----------------------------------------------------------------------------
# STEP 5: MANUSCRIPT Figure 4 -- Effect Size Concordance
# -----------------------------------------------------------------------------

cat("Generating Manuscript Figure 4...\n")

het_dt      <- as.data.table(het)
het_all     <- het_dt[!is.na(beta_CSA)]
het_reliable <- het_all[csa_absz >= 1]

r_all      <- cor(het_all$beta_EUR,      het_all$beta_CSA,      use="complete.obs")
r_reliable <- cor(het_reliable$beta_EUR, het_reliable$beta_CSA, use="complete.obs")

ms_fig4 <- ggplot(het_all,
                  aes(x = beta_EUR, y = beta_CSA,
                      color = -log10(Q_pval))) +
  geom_hline(yintercept = 0, linetype="dotted", color="gray70") +
  geom_vline(xintercept = 0, linetype="dotted", color="gray70") +
  geom_abline(slope=1, intercept=0,
              linetype="dashed", color="gray40", linewidth=0.7) +
  geom_point(alpha = 0.4, size = 1.0) +
  scale_color_viridis_c(
    name      = "-log10(Q p-value)\nheterogeneity",
    option    = "plasma",
    direction = -1
  ) +
  annotate("text", x = -Inf, y = Inf,
           hjust = -0.1, vjust = 1.8, size = 3.8,
           label = paste0(
             "r = ", round(r_all, 3), " (all loci, n=",
             nrow(het_all), ")\n",
             "r = ", round(r_reliable, 3),
             " (CSA |Z|≥1, n=", nrow(het_reliable), ")"
           )) +
  labs(
    title    = "EUR vs CSA Effect Size Concordance at EUR GWS T2D Loci",
    subtitle = paste0(
      "No loci reach Bonferroni-significant heterogeneity | ",
      "Moderate correlation consistent with shared architecture (rg≈1.27)\n",
      "Systematic attenuation of CSA effects reflects winner's curse ",
      "in underpowered CSA GWAS"
    ),
    x        = "Effect Size in EUR (β)",
    y        = "Effect Size in CSA (β)",
    caption  = "Dashed line = perfect concordance (slope=1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray40", size = 9),
    panel.grid.minor = element_blank()
  )

ggsave("results/figures/manuscript/ms_fig4_effect_size_concordance.png",
       plot = ms_fig4, width = 8, height = 7, dpi = 300)
cat("Saved: manuscript/ms_fig4_effect_size_concordance.png\n")

# -----------------------------------------------------------------------------
# STEP 6: MANUSCRIPT Figure 5 -- Genetic Correlation Summary
# -----------------------------------------------------------------------------

cat("Generating Manuscript Figure 5...\n")

rg     <- rg_result$rg$rg
rg_se  <- rg_result$rg$rg_se
rg_p   <- rg_result$rg$rg_p
h2_eur_obs <- rg_result$h2$h2_observed[1]
h2_csa_obs <- rg_result$h2$h2_observed[2]

plot_df <- data.table(
  Metric   = factor(
    c("EUR h2\n(observed)", "CSA h2\n(observed)", "EUR-CSA\nrg"),
    levels = c("EUR h2\n(observed)", "CSA h2\n(observed)", "EUR-CSA\nrg")
  ),
  Estimate = c(h2_eur_obs, h2_csa_obs, rg),
  SE       = c(rg_result$h2$h2_observed_se[1],
               rg_result$h2$h2_observed_se[2],
               rg_se),
  Category = c("Heritability", "Heritability", "Correlation"),
  Label    = c(
    paste0("h²=", round(h2_eur_obs, 3), "\np=8.1e-41"),
    paste0("h²=", round(h2_csa_obs, 3), "\np=0.49 (ns)"),
    paste0("rg=", round(rg, 3), "\np=6.7e-11")
  )
)
plot_df[, Lower := Estimate - 1.96 * SE]
plot_df[, Upper := pmin(Estimate + 1.96 * SE, 1.5)]

ms_fig5 <- ggplot(plot_df,
                  aes(x = Metric, y = Estimate, color = Category)) +
  geom_hline(yintercept = 0,
             linetype = "dashed", color = "gray60") +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.15, linewidth = 1.0) +
  geom_text(aes(label = Label, y = Upper + 0.05),
            size = 3.2, color = "gray30", vjust = 0) +
  scale_color_manual(values = c(
    "Heritability" = "#4a90d9",
    "Correlation"  = "#e07b39"
  )) +
  facet_wrap(~Category, scales = "free") +
  labs(
    x        = NULL,
    y        = "Estimate ± 95% CI"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray40", size = 9),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  ) +
  coord_cartesian(clip = "off") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

ggsave("results/figures/genetic_correlation.png",
       plot = ms_fig5, width = 9, height = 5.5, dpi = 300)
cat("Saved: manuscript/ms_fig5_genetic_correlation.png\n")

# -----------------------------------------------------------------------------
# STEP 7: TALK Figures -- Simplified, large fonts
# -----------------------------------------------------------------------------

cat("\nGenerating Conference Talk figures...\n")

talk_theme <- theme_minimal(base_size = 18) +
  theme(
    plot.title       = element_text(face = "bold", size = 20),
    plot.subtitle    = element_text(color = "gray40", size = 14),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 14)
  )

# Talk Figure 1 -- The key clinical finding (decile analysis, C+T only)
talk_decile <- decile_plot_df[method == "C+T"]

talk_fig1 <- ggplot(talk_decile,
                    aes(x = factor(decile), y = Proportion_pct,
                        fill = Population)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
  geom_hline(yintercept = 10,
             linetype = "dashed", color = "gray40",
             linewidth = 1.0) +
  scale_fill_manual(values  = pop_colors,
                    labels  = c("EUR" = "European",
                                "SAS" = "South Asian")) +
  annotate("text", x = 10, y = 42, hjust = 1, size = 5,
           color = "#e07b39", fontface = "bold",
           label = "40% of SAS in\ntop risk decile!") +
  annotate("text", x = 5.5, y = 11.5, size = 4,
           color = "gray40", label = "Expected: 10%") +
  labs(
    title    = "EUR-Derived T2D PRS Systematically Misclassifies South Asians",
    subtitle = "40% of SAS individuals fall in top risk decile vs 10% expected",
    x        = "PRS Risk Decile (1=lowest, 10=highest)",
    y        = "% of Individuals",
    fill     = NULL
  ) +
  talk_theme

ggsave("results/figures/talk/talk_fig1_decile.png",
       plot = talk_fig1, width = 12, height = 7, dpi = 200)
cat("Saved: talk/talk_fig1_decile.png\n")

# Talk Figure 2 -- PRS distribution (clean single panel)
talk_df <- data.table(
  PRS = c(ct_EUR$PRS_std_EUR, ct_SAS$PRS_std_EUR),
  Pop = c(rep("EUR", nrow(ct_EUR)), rep("SAS", nrow(ct_SAS)))
)

talk_fig2 <- ggplot(talk_df, aes(x = PRS, fill = Pop, color = Pop)) +
  geom_density(alpha = 0.4, linewidth = 1.0) +
  scale_fill_manual(values  = pop_colors,
                    labels  = c("EUR" = "European", "SAS" = "South Asian")) +
  scale_color_manual(values = pop_colors) +
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_vline(xintercept = mean(ct_SAS$PRS_std_EUR),
             linetype = "dashed", color = "#e07b39", linewidth = 1.0) +
  annotate("text",
           x = mean(ct_SAS$PRS_std_EUR) + 0.1,
           y = 0.38, hjust = 0, size = 5.5,
           color = "#e07b39", fontface = "bold",
           label = paste0("+", round(mean(ct_SAS$PRS_std_EUR), 2),
                          " SDs shift")) +
  labs(
    title    = "EUR-Derived T2D PRS: South Asians Systematically Over-Scored",
    subtitle = "EUR PRS standardized to mean=0, SD=1 and applied to South Asian individuals",
    x        = "Standardized PRS",
    y        = "Density",
    fill     = NULL,
    color    = NULL
  ) +
  talk_theme

ggsave("results/figures/talk/talk_fig2_distribution.png",
       plot = talk_fig2, width = 12, height = 7, dpi = 200)
cat("Saved: talk/talk_fig2_distribution.png\n")

# Talk Figure 3 -- Method comparison (both methods confirm finding)
method_df <- summary_stats[
  scenario == "EUR_effects_EUR_LD",
  .(method, mean_shift)
]

talk_fig3 <- ggplot(method_df,
                    aes(x = method, y = mean_shift, fill = method)) +
  geom_col(alpha = 0.85, width = 0.5) +
  geom_text(aes(label = paste0("+", round(mean_shift, 2), " SDs")),
            vjust = -0.5, size = 6, fontface = "bold") +
  scale_fill_manual(values = method_colors) +
  ylim(0, 1.7) +
  labs(
    title    = "Miscalibration Is Method-Independent",
    subtitle = "Both C+T and LDpred2-inf confirm ~1.3 SD systematic over-scoring",
    x        = "PRS Method",
    y        = "Mean Shift in SAS (SDs)",
    fill     = NULL
  ) +
  talk_theme +
  theme(legend.position = "none")

ggsave("results/figures/talk/talk_fig3_methods.png",
       plot = talk_fig3, width = 8, height = 7, dpi = 200)
cat("Saved: talk/talk_fig3_methods.png\n")

# -----------------------------------------------------------------------------
# STEP 8: POSTER -- Composite figure
# -----------------------------------------------------------------------------

cat("\nGenerating Poster composite figure...\n")

poster_theme <- theme_minimal(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 10),
    plot.subtitle    = element_text(color = "gray40", size = 8),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.4, "cm")
  )

# Panel A: PRS distribution (C+T)
pa <- ggplot(talk_df, aes(x = PRS, fill = Pop, color = Pop)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  scale_fill_manual(values  = pop_colors,
                    labels  = c("EUR"="EUR (n=503)",
                                "SAS"="SAS (n=478)")) +
  scale_color_manual(values = pop_colors) +
  geom_vline(xintercept = 0, linetype="dashed",
             color="gray50", linewidth=0.5) +
  geom_vline(xintercept = mean(ct_SAS$PRS_std_EUR),
             linetype="dashed", color="#e07b39", linewidth=0.7) +
  annotate("text",
           x = mean(ct_SAS$PRS_std_EUR) + 0.1,
           y = Inf, vjust=1.5, hjust=0, size=2.8,
           color="#e07b39",
           label=paste0("+", round(mean(ct_SAS$PRS_std_EUR),2), " SDs")) +
  labs(title="A. PRS Distribution (C+T)",
       x="Standardized PRS", y="Density",
       fill=NULL, color=NULL) +
  poster_theme

# Panel B: Decile analysis (C+T only)
pb <- ggplot(talk_decile,
             aes(x=factor(decile), y=Proportion_pct, fill=Population)) +
  geom_col(position="dodge", alpha=0.85, width=0.7) +
  geom_hline(yintercept=10, linetype="dashed",
             color="gray40", linewidth=0.5) +
  scale_fill_manual(values = pop_colors,
                    labels = c("EUR"="EUR","SAS"="SAS")) +
  labs(title="B. Risk Decile Distribution",
       x="PRS Decile", y="% Individuals",
       fill=NULL) +
  poster_theme

# Panel C: Directional concordance
pc_data <- consist_dt
pc <- ggplot(pc_data, aes(x=z_bin, y=pct, group=1)) +
  geom_hline(yintercept=50, linetype="dashed",
             color="gray60", linewidth=0.5) +
  geom_ribbon(aes(ymin=lower_ci, ymax=upper_ci),
              fill="#4a90d9", alpha=0.2) +
  geom_line(color="#4a90d9", linewidth=0.9) +
  geom_point(color="#4a90d9", size=2.5) +
  geom_text(aes(label=paste0(round(pct,0),"%")),
            vjust=-1.0, size=2.5, color="gray30") +
  scale_y_continuous(limits=c(40,110),
                     breaks=c(50,70,90)) +
  labs(title="C. Directional Concordance by CSA Power",
       x="CSA |Z-score| Bin",
       y="% Same Direction") +
  poster_theme

# Panel D: Effect size concordance (subsample for clarity)
set.seed(42)
het_sub <- het_all[sample(.N, min(2000, .N))]
r_sub   <- cor(het_sub$beta_EUR, het_sub$beta_CSA, use="complete.obs")

pd <- ggplot(het_sub,
             aes(x=beta_EUR, y=beta_CSA,
                 color=-log10(Q_pval))) +
  geom_abline(slope=1, intercept=0,
              linetype="dashed", color="gray40", linewidth=0.5) +
  geom_point(alpha=0.3, size=0.8) +
  scale_color_viridis_c(option="plasma", direction=-1,
                        name="-log10(Q)") +
  annotate("text", x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
           size=2.8, label=paste0("r=", round(r_sub,3))) +
  labs(title="D. Effect Size Concordance",
       x="EUR beta", y="CSA beta") +
  poster_theme +
  theme(legend.position="none")

# Panel E: Genetic correlation
pe_df <- data.table(
  Metric   = factor(c("EUR h²", "CSA h²", "rg"),
                    levels=c("EUR h²","CSA h²","rg")),
  Estimate = c(h2_eur_obs, h2_csa_obs, rg),
  SE       = c(rg_result$h2$h2_observed_se[1],
               rg_result$h2$h2_observed_se[2], rg_se),
  Cat      = c("h²","h²","rg")
)
pe_df[, Lower := Estimate - 1.96*SE]
pe_df[, Upper := pmin(Estimate + 1.96*SE, 1.5)]

pe <- ggplot(pe_df, aes(x=Metric, y=Estimate, color=Cat)) +
  geom_hline(yintercept=0, linetype="dashed", color="gray60") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper),
                width=0.15, linewidth=0.7) +
  scale_color_manual(values=c("h²"="#4a90d9","rg"="#e07b39")) +
  labs(title="E. LDSC Heritability & rg",
       x=NULL, y="Estimate ± 95% CI") +
  poster_theme +
  theme(legend.position="none")

# Panel F: GO enrichment (top 10)
go_dt  <- as.data.table(go_all)
top_go <- go_dt[order(p.adjust)][1:10]
top_go[, GR_num := sapply(GeneRatio, function(x) {
  parts <- as.numeric(strsplit(x,"/")[[1]])
  parts[1]/parts[2]
})]
top_go[, Desc_short := substr(Description, 1, 35)]

pf <- ggplot(top_go,
             aes(x=GR_num,
                 y=reorder(Desc_short, -p.adjust),
                 color=p.adjust, size=Count)) +
  geom_point(alpha=0.8) +
  scale_color_viridis_c(option="plasma", direction=-1,
                        name="Adj.p") +
  scale_size_continuous(range=c(2,6), name="Genes") +
  labs(title="F. GO Pathway Enrichment",
       x="Gene Ratio", y=NULL) +
  poster_theme +
  theme(axis.text.y=element_text(size=7))

# Assemble poster
poster <- (pa | pb | pc) / (pd | pe | pf) +
  patchwork::plot_annotation(
    title    = "EUR-Derived T2D Polygenic Risk Score Miscalibration in South Asian Populations",
    subtitle = paste0(
      "Pan-UK Biobank E11 GWAS | 1000 Genomes Phase 3 | ",
      "EUR n=503, SAS n=478 | Both C+T and LDpred2-inf methods"
    ),
    theme = theme(
      plot.title    = element_text(face="bold", size=14, hjust=0.5),
      plot.subtitle = element_text(size=9, hjust=0.5, color="gray40")
    )
  )

ggsave("results/figures/poster/poster_composite.png",
       plot = poster, width = 18, height = 12, dpi = 300)
cat("Saved: poster/poster_composite.png\n\n")

# -----------------------------------------------------------------------------
# STEP 9: Figure inventory summary
# -----------------------------------------------------------------------------

cat("Generating figure inventory...\n\n")

inventory <- data.table(
  Figure   = c("ms_fig1", "ms_fig2", "ms_fig3",
               "ms_fig4", "ms_fig5",
               "talk_fig1", "talk_fig2", "talk_fig3",
               "poster_composite"),
  Type     = c(rep("Manuscript", 5),
               rep("Talk", 3),
               "Poster"),
  Content  = c("PRS distributions (C+T + LDpred2-inf)",
               "Decile analysis",
               "Directional concordance by Z-score bin",
               "Effect size concordance scatter",
               "Genetic correlation (LDSC)",
               "Decile analysis -- clinical impact",
               "PRS distribution -- key result",
               "Method comparison -- robustness",
               "All key results composite"),
  Key_stat = c("+1.32 SDs mean shift",
               "40% SAS in top decile",
               "100% concordance at |Z|>3",
               "r=0.42 (all), r=0.61 (powered)",
               "rg=1.27, p=6.7e-11",
               "40% vs 10% expected",
               "+1.32 SDs mean shift",
               "C+T +1.32, LDpred2 +1.27",
               "Complete summary")
)

fwrite(inventory, "results/figures/figure_inventory.csv")
cat("Figure inventory saved.\n\n")
print(inventory)

# -----------------------------------------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------------------------------------
cat("============================================================\n")
cat("SCRIPT 10 COMPLETE -- ALL ANALYSES DONE\n")
cat("============================================================\n")
cat("\nMANUSCRIPT FIGURES (results/figures/manuscript/):\n")
cat("  ms_fig1: PRS distributions (main result)\n")
cat("  ms_fig2: Decile analysis (clinical impact)\n")
cat("  ms_fig3: Directional concordance (QC finding)\n")
cat("  ms_fig4: Effect size concordance (locus level)\n")
cat("  ms_fig5: Genetic correlation (architecture)\n")
cat("\nCONFERENCE TALK FIGURES (results/figures/talk/):\n")
cat("  talk_fig1: Decile analysis (lead with clinical impact)\n")
cat("  talk_fig2: PRS distribution (core finding)\n")
cat("  talk_fig3: Method comparison (robustness)\n")
cat("\nPOSTER (results/figures/poster/):\n")
cat("  poster_composite: All 6 key panels in one figure\n")
cat("\nKEY FINDINGS SUMMARY:\n")
cat("  PRS mean shift (C+T):       +1.32 SDs\n")
cat("  PRS mean shift (LDpred2):   +1.27 SDs\n")
cat("  SAS in top risk decile:      39.5% vs 10% expected\n")
cat("  EUR-CSA genetic correlation: rg = 1.27 (p = 6.7e-11)\n")
cat("  Directional concordance:     68% overall, 100% at |Z|>3\n")
cat("  Effect size correlation:     r = 0.42 (all), 0.61 (powered)\n")
cat("  Bonferroni heterogeneous:    0 loci\n")
cat("  Mean |AF difference|:        0.076\n")
cat("============================================================\n")