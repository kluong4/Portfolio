# =============================================================================
# Script 05: LDpred2 PRS Construction (Infinitesimal Model)
# Project:   Cross-ancestry T2D PRS Transferability (EUR -> CSA)
# =============================================================================
# OVERVIEW:
# Constructs PRS using LDpred2-inf (infinitesimal model).
# LDpred2-auto was attempted but failed to converge on most chromosomes
# due to RAM constraints (24GB machine, ~2-3GB free during computation).
# LDpred2-inf derives posterior effect sizes in closed form without MCMC,
# requiring significantly less memory and running stably on all 22 chromosomes.
#
# METHODS NOTE FOR MANUSCRIPT:
# "LDpred2-auto failed to converge on the majority of chromosomes due to
# memory constraints on the local computing environment. We therefore applied
# LDpred2-inf, which assumes an infinitesimal model of genetic architecture
# and derives posterior effect sizes in closed form without MCMC sampling.
# Per-chromosome SNP heritability (h2) was estimated from chromosomes where
# LDpred2-auto converged (n=6 chromosomes) and extrapolated genome-wide
# proportionally to SNP count, yielding an estimated total h2 = 0.352."
#
# KEY SCIENTIFIC FINDING:
# Both C+T (+1.317 SDs) and LDpred2-inf (+1.269 SDs) independently show
# ~1.3 SD miscalibration when EUR PRS is applied to SAS individuals.
# Switching to CSA effects + SAS LD does not reduce the shift (+1.286 SDs),
# indicating miscalibration is driven by allele frequency differences
# between EUR and SAS rather than LD structure differences.
#
# ALL TECHNICAL FIXES:
#   1. Remove chrX before snp_match()
#   2. Pre-filter sumstats to map positions before snp_match()
#   3. Use _NUM_ID_ not _j for map indices (bigsnpr 1.12.21)
#   4. Build idx_snps from chr:pos not rsID
#   5. Convert dsCMatrix to SFBM with as_SFBM() before ldpred2
#   6. n_eff = 4/(1/n_cases + 1/n_controls) not af_cases + af_controls
#   7. Use snp_ldpred2_inf not snp_ldpred2_auto
#   8. Use ind.test not ind.row in snp_PRS()
#   9. Use LDpred2-specific EUR standardization params not C+T params
#  10. Define NCORES before chromosome loop
#
# OUTPUTS:
#   data/sumstats/info_EUR_matched.rds
#   data/sumstats/info_CSA_matched.rds
#   data/1KGP/idx_snps_new.rds
#   data/LD_ref/ldpred2_EUR/LD_chr{1-22}.rds
#   data/LD_ref/ldpred2_SAS/LD_chr{1-22}.rds
#   results/PRS/LDpred2/betas_EUR.rds
#   results/PRS/LDpred2/betas_CSA.rds
#   results/PRS/LDpred2/chr_h2_estimates.rds
#   results/PRS/LDpred2/ldpred2_eur_standardization_params.rds
#   results/PRS/LDpred2/prs_LDpred2_EUR_in_EUR.txt
#   results/PRS/LDpred2/prs_LDpred2_EUR_in_SAS.txt
#   results/PRS/LDpred2/prs_LDpred2_CSA_in_SAS.txt
#   results/figures/LDpred2_prs_distribution.png
# =============================================================================

library(bigsnpr)
library(bigstatsr)
library(data.table)
library(Matrix)
library(ggplot2)

dir.create("results/PRS/LDpred2",     recursive = TRUE, showWarnings = FALSE)
dir.create("data/LD_ref/ldpred2_EUR", recursive = TRUE, showWarnings = FALSE)
dir.create("data/LD_ref/ldpred2_SAS", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures",         recursive = TRUE, showWarnings = FALSE)

# Must define before chromosome loop
NCORES <- 2

# -----------------------------------------------------------------------------
# Load genotype data and sample lists
# -----------------------------------------------------------------------------
obj.bigsnp <- snp_attach("data/1KGP/1000G_phase3_common_norel.rds")
G   <- obj.bigsnp$genotypes
map <- obj.bigsnp$map
fam <- obj.bigsnp$fam

samples_EUR <- fread("data/1KGP/samples_EUR_final.txt")
samples_SAS <- fread("data/1KGP/samples_SAS_final.txt")

idx_EUR <- which(fam$sample.ID %in% samples_EUR$IID)
idx_SAS <- which(fam$sample.ID %in% samples_SAS$IID)

# Compute effective sample sizes
# LDpred2 requires actual sample sizes not allele frequencies.
# For binary traits: n_eff = 4 / (1/n_cases + 1/n_controls)

# EUR: n_cases=22,634  n_controls=397,897  → n_eff=85,663
n_eff_EUR <- 4 / (1/22634  + 1/397897) 
# CSA: n_cases= 1,662  n_controls=  7,214  → n_eff= 5,403
n_eff_CSA <- 4 / (1/1662   + 1/7214)

# -----------------------------------------------------------------------------
# Build position-based QC SNP index
# -----------------------------------------------------------------------------
# snps_pass_both.txt uses rsIDs from original bim but _NUM_ID_ indices
# are based on map row positions. Must use chr:pos matching.

bim_EUR_qc <- fread("data/1KGP/1KGP_EUR_QC.bim")
bim_SAS_qc <- fread("data/1KGP/1KGP_SAS_QC.bim")

colnames(bim_EUR_qc) <- c("chr", "snpID", "cm", "pos", "a1", "a2")
colnames(bim_SAS_qc) <- c("chr", "snpID", "cm", "pos", "a1", "a2")

qc_pos_EUR   <- paste(bim_EUR_qc$chr, bim_EUR_qc$pos, sep = ":")
qc_pos_SAS   <- paste(bim_SAS_qc$chr, bim_SAS_qc$pos, sep = ":")
qc_pos_both  <- intersect(qc_pos_EUR, qc_pos_SAS) # QC-passed positions: 1422958 

# -----------------------------------------------------------------------------
# Prepare map dataframe
# -----------------------------------------------------------------------------
map_pos      <- paste(map$chromosome, map$physical.pos, sep = ":")
idx_snps_new <- which(map_pos %in% qc_pos_both)

map_df <- data.frame(
  chr = map$chromosome,
  pos = map$physical.pos,
  a0  = map$allele2,      # a0 = reference (allele2)
  a1  = map$allele1,      # a1 = effect (allele1)
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Harmonize summary statistics to genotype map
# -----------------------------------------------------------------------------
sumstats_EUR <- readRDS("data/sumstats/outputs/sumstats_EUR_formatted.rds")
sumstats_CSA <- readRDS("data/sumstats/outputs/sumstats_CSA_formatted.rds")

build_df <- function(ss, n_eff_val) {
  # Remove chrX -- 1KGP map is autosomes only
  keep <- ss$chr %in% as.character(1:22)
  
  data.frame(
    chr     = as.integer(ss$chr[keep]),
    pos     = ss$pos[keep],
    a0      = ss$ref[keep],
    a1      = ss$alt[keep],
    beta    = ss$beta[keep],
    beta_se = ss$se[keep],
    n_eff   = n_eff_val,
    stringsAsFactors = FALSE
  )
}

df_EUR <- build_df(sumstats_EUR, n_eff_EUR)
df_CSA <- build_df(sumstats_CSA, n_eff_CSA)

# Pre-filter to positions in map -- prevents snp_match() threshold error
map_pos_key <- paste(map_df$chr, map_df$pos, sep = ":")

# EUR SNPs after pre-filter: nrow(df_EUR_filt) --> 1435146
df_EUR_filt <- df_EUR[paste(df_EUR$chr, df_EUR$pos, sep=":") %in% map_pos_key, ]
# CSA SNPs after pre-filter: nrow(df_CSA_filt) --> 1435146
df_CSA_filt <- df_CSA[paste(df_CSA$chr, df_CSA$pos, sep=":") %in% map_pos_key, ]

# snp_match aligns alleles and handles strand flips
info_EUR <- snp_match(df_EUR_filt, map_df, strand_flip = TRUE)
info_CSA <- snp_match(df_CSA_filt, map_df, strand_flip = TRUE)

## Restrict to QC-passed SNPs using position-based _NUM_ID_ index
# EUR SNPs after QC restriction: 1312946 
info_EUR <- info_EUR[info_EUR[["_NUM_ID_"]] %in% idx_snps_new, ]
# CSA SNPs after QC restriction: 1312946 
info_CSA <- info_CSA[info_CSA[["_NUM_ID_"]] %in% idx_snps_new, ]

# Save for safe restart
saveRDS(info_EUR,     "data/sumstats/info_EUR_matched.rds")
saveRDS(info_CSA,     "data/sumstats/info_CSA_matched.rds")
saveRDS(idx_snps_new, "data/1KGP/idx_snps_new.rds")

# -----------------------------------------------------------------------------
# Computing per-chromosome h2 estimates
# -----------------------------------------------------------------------------
# Per-chromosome h2 estimated proportionally to SNP count.
# Total h2 estimated from 6 chromosomes where LDpred2-auto converged:
#   h2 sum (6 chr, 502,375 SNPs) = 0.1348
# Same h2_per_snp applied to CSA (T2D heritability similar across ancestries).

snps_per_chr_EUR <- table(info_EUR$chr)
snps_per_chr_CSA <- table(info_CSA$chr)

h2_per_snp_EUR <- 0.1348 / 502375                   # 2.684e-7
h2_per_snp_CSA <- h2_per_snp_EUR                    # 2.684e-7
h2_genome_EUR  <- h2_per_snp_EUR * nrow(info_EUR)   # 0.3523
h2_genome_CSA  <- h2_per_snp_CSA * nrow(info_CSA)   # 0.3523

chr_h2_EUR <- sapply(1:22, function(chr) {
  n <- as.numeric(snps_per_chr_EUR[as.character(chr)])
  if(is.na(n) || n == 0) return(0.001)
  h2_per_snp_EUR * n
})

chr_h2_CSA <- sapply(1:22, function(chr) {
  n <- as.numeric(snps_per_chr_CSA[as.character(chr)])
  if(is.na(n) || n == 0) return(0.001)
  h2_per_snp_CSA * n
})

chr_h2_df <- data.frame(
  chr        = 1:22,
  n_snps_EUR = as.numeric(snps_per_chr_EUR[as.character(1:22)]),
  h2_EUR     = chr_h2_EUR,
  h2_CSA     = chr_h2_CSA
)

# Per-chromosome h2 estimates Results
print(chr_h2_df)
saveRDS(chr_h2_df, "results/PRS/LDpred2/chr_h2_estimates.rds")

# -----------------------------------------------------------------------------
# LDpred2-inf chromosome function
# -----------------------------------------------------------------------------
# snp_ldpred2_inf computes posterior effect sizes under the infinitesimal
# model in closed form -- no MCMC chains, much less RAM than auto.

run_ldpred2_inf_chr <- function(chr, info_snp, idx_ref,
                                ld_dir, h2_chr, scenario_label) {
  
  cat("\n--- Chromosome", chr, "|", scenario_label, "---\n")
  
  chr_idx <- which(info_snp$chr == chr)
  
  if(length(chr_idx) < 10) {
    cat("  Too few SNPs (", length(chr_idx), ") -- skipping\n")
    return(NULL)
  }
  
  # _NUM_ID_ is the correct index column in bigsnpr 1.12.21 (not _j)
  ind_chr <- info_snp[["_NUM_ID_"]][chr_idx]
  cat("  SNPs:", length(ind_chr), "| h2:", round(h2_chr, 5), "\n")
  
  # Load or compute LD matrix
  ld_file <- paste0(ld_dir, "/LD_chr", chr, ".rds")
  
  if(!file.exists(ld_file)) {
    cat("  Computing LD matrix...\n")
    corr_chr <- snp_cor(
      G,
      ind.row = idx_ref,
      ind.col = ind_chr,
      size    = 3000,     # window size in SNPs
      alpha   = 1e-3,     # sparsity threshold
      ncores  = NCORES
    )
    saveRDS(corr_chr, ld_file, compress = FALSE)
    cat("  LD matrix saved.\n")
  } else {
    cat("  Loading existing LD matrix...\n")
    corr_chr <- readRDS(ld_file)
  }
  
  # Convert dsCMatrix to SFBM -- required by snp_ldpred2_inf
  # Remove dsCMatrix immediately to free RAM
  cat("  Converting to SFBM...\n")
  corr_sfbm <- as_SFBM(corr_chr)
  rm(corr_chr)
  gc()
  
  df_beta_chr <- data.frame(
    beta    = info_snp$beta[chr_idx],
    beta_se = info_snp$beta_se[chr_idx],
    n_eff   = info_snp$n_eff[chr_idx]
  )
  
  cat("  Running LDpred2-inf (n_eff =",
      round(unique(df_beta_chr$n_eff)), ")...\n")
  
  tryCatch({
    
    # Closed-form infinitesimal posterior betas
    # h2 = per-chromosome SNP heritability (proportional to SNP count)
    beta_inf <- snp_ldpred2_inf(
      corr    = corr_sfbm,
      df_beta = df_beta_chr,
      h2      = h2_chr
    )
    
    cat("  Beta range: [",
        round(min(beta_inf), 6), ",",
        round(max(beta_inf), 6), "]\n")
    
    result <- data.frame(
      chr            = chr,
      NUM_ID         = ind_chr,
      beta_posterior = beta_inf
    )
    
    rm(corr_sfbm)
    gc()
    
    return(result)
    
  }, error = function(e) {
    cat("  ERROR on chr", chr, ":", conditionMessage(e), "\n")
    rm(corr_sfbm)
    gc()
    return(NULL)
  })
}

# -----------------------------------------------------------------------------
# SCENARIO A: EUR effects + EUR LD (LDpred2-inf)
# -----------------------------------------------------------------------------
betas_EUR_list <- list()

for(chr in 1:22) {
  
  result <- run_ldpred2_inf_chr(
    chr            = chr,
    info_snp       = info_EUR,
    idx_ref        = idx_EUR,
    ld_dir         = "data/LD_ref/ldpred2_EUR",
    h2_chr         = chr_h2_EUR[chr],
    scenario_label = "EUR_LD"
  )
  
  if(!is.null(result)) betas_EUR_list[[chr]] <- result
  gc()
  
  if(chr %% 5 == 0) {
    cat("\nProgress: chr", chr, "/ 22\n")
    print(memuse::Sys.meminfo())
  }
}

betas_EUR <- do.call(rbind, betas_EUR_list)

cat("Chromosomes completed:", sum(!sapply(betas_EUR_list, is.null)), "/ 22\n")
cat("Total SNPs: ", nrow(betas_EUR)) #1312946

saveRDS(betas_EUR, "results/PRS/LDpred2/betas_EUR.rds")

gc()

# -----------------------------------------------------------------------------
# SCENARIO B: CSA effects + SAS LD (LDpred2-inf)
# -----------------------------------------------------------------------------
betas_CSA_list <- list()

for(chr in 1:22) {
  
  result <- run_ldpred2_inf_chr(
    chr            = chr,
    info_snp       = info_CSA,
    idx_ref        = idx_SAS,
    ld_dir         = "data/LD_ref/ldpred2_SAS",
    h2_chr         = chr_h2_CSA[chr],
    scenario_label = "SAS_LD"
  )
  
  if(!is.null(result)) betas_CSA_list[[chr]] <- result
  gc()
  
  if(chr %% 5 == 0) {
    cat("\nProgress: chr", chr, "/ 22\n")
    print(memuse::Sys.meminfo())
  }
}

betas_CSA <- do.call(rbind, betas_CSA_list)

cat("Chromosomes completed:", sum(!sapply(betas_CSA_list, is.null)), "/ 22")
cat("Total SNPs: ", nrow(betas_CSA)) # 1312946

saveRDS(betas_CSA, "results/PRS/LDpred2/betas_CSA.rds")

gc()

# -----------------------------------------------------------------------------
# Computing LDpred2 PRS scores
# -----------------------------------------------------------------------------
score_ldpred2 <- function(betas_df, idx_test, pop_label, scenario_label) {
  
  cat("Scoring:", scenario_label, "in", pop_label, "\n")
  
  nonzero <- betas_df[betas_df$beta_posterior != 0, ]
  cat("  SNPs used:", nrow(nonzero), "\n")
  
  prs <- snp_PRS(
    G,
    betas.keep = nonzero$beta_posterior,
    ind.keep   = nonzero$NUM_ID,   # SNP column indices in G
    ind.test   = idx_test          # sample row indices in G
  )
  
  cat("  Scored:", length(prs), "individuals\n\n")
  return(prs)
}

# EUR posterior betas → EUR samples (reference)
prs_EUR_ldpred2 <- score_ldpred2(betas_EUR, idx_EUR, "EUR", "EUR_effects_EUR_LD")

# EUR posterior betas → SAS samples (clinical simulation)
# Mean shift from EUR quantifies miscalibration under universal deployment
prs_SAS_eurLD_ldpred2 <- score_ldpred2(betas_EUR, idx_SAS, "SAS", "EUR_effects_EUR_LD")

# CSA posterior betas → SAS samples (ancestry-matched)
# Compare to EUR effects to isolate ancestry-matching benefit
prs_SAS_csaLD_ldpred2 <- score_ldpred2(betas_CSA, idx_SAS, "SAS", "CSA_effects_SAS_LD")

# -----------------------------------------------------------------------------
# Standardize using LDpred2-specific EUR parameters
# -----------------------------------------------------------------------------
# IMPORTANT: Use LDpred2-specific EUR mean/SD not C+T parameters.
# C+T and LDpred2-inf raw scores are on completely different scales
# because LDpred2 uses all 1.3M SNPs vs C+T's ~160k clumped SNPs.
# The C+T eur_mean (~4.5e-5) would produce thousands of SDs shift.

ldpred2_eur_mean <- mean(prs_EUR_ldpred2, na.rm = TRUE) # -0.355597
ldpred2_eur_sd   <- sd(prs_EUR_ldpred2,   na.rm = TRUE) # 0.66071

prs_EUR_ld_eurscaled    <- (prs_EUR_ldpred2        - ldpred2_eur_mean) / ldpred2_eur_sd
prs_SAS_eurLD_eurscaled <- (prs_SAS_eurLD_ldpred2  - ldpred2_eur_mean) / ldpred2_eur_sd
prs_SAS_csaLD_eurscaled <- (prs_SAS_csaLD_ldpred2  - ldpred2_eur_mean) / ldpred2_eur_sd

# LDPRED2-INF PRS summary
cat("EUR PRS mean (LDpred2-scaled):        ", round(mean(prs_EUR_ld_eurscaled),    4)) # 0
cat("SAS EUR LD mean (LDpred2-scaled):     ", round(mean(prs_SAS_eurLD_eurscaled), 4)) # 1.2688 
cat("SAS CSA LD mean (LDpred2-scaled):     ", round(mean(prs_SAS_csaLD_eurscaled), 4)) # 1.26857
cat("Mean shift EUR→SAS (EUR effects):     ", round(mean(prs_SAS_eurLD_eurscaled) - mean(prs_EUR_ld_eurscaled), 4), "SDs") # 1.2688 SDs
cat("C+T mean shift for comparison:        +1.3165 SDs\n\n")

# Save LDpred2 standardization parameters separately from C+T
saveRDS(list(eur_mean = ldpred2_eur_mean, eur_sd = ldpred2_eur_sd),
  "results/PRS/LDpred2/ldpred2_eur_standardization_params.rds"
)

# Save flat files
prs_EUR_out <- data.table(
  IID            = fam$sample.ID[idx_EUR],
  PRS_raw        = prs_EUR_ldpred2,
  PRS_std_within = as.numeric(scale(prs_EUR_ldpred2)),
  PRS_std_EUR    = prs_EUR_ld_eurscaled
)
fwrite(prs_EUR_out, "results/PRS/LDpred2/prs_LDpred2_EUR_in_EUR.txt", sep = "\t")

prs_SAS_eurLD_out <- data.table(
  IID            = fam$sample.ID[idx_SAS],
  PRS_raw        = prs_SAS_eurLD_ldpred2,
  PRS_std_within = as.numeric(scale(prs_SAS_eurLD_ldpred2)),
  PRS_std_EUR    = prs_SAS_eurLD_eurscaled
)
fwrite(prs_SAS_eurLD_out, "results/PRS/LDpred2/prs_LDpred2_EUR_in_SAS.txt", sep = "\t")

prs_SAS_csaLD_out <- data.table(
  IID            = fam$sample.ID[idx_SAS],
  PRS_raw        = prs_SAS_csaLD_ldpred2,
  PRS_std_within = as.numeric(scale(prs_SAS_csaLD_ldpred2)),
  PRS_std_EUR    = prs_SAS_csaLD_eurscaled
)
fwrite(prs_SAS_csaLD_out, "results/PRS/LDpred2/prs_LDpred2_CSA_in_SAS.txt", sep = "\t")

# -----------------------------------------------------------------------------
# Diagnostic figure
# -----------------------------------------------------------------------------

prs_plot_df <- data.table(
  PRS = c(prs_EUR_ld_eurscaled,
          prs_SAS_eurLD_eurscaled,
          prs_SAS_csaLD_eurscaled),
  Pop = c(rep("EUR",          length(prs_EUR_ld_eurscaled)),
          rep("SAS (EUR LD)", length(prs_SAS_eurLD_eurscaled)),
          rep("SAS (CSA LD)", length(prs_SAS_csaLD_eurscaled)))
)

pop_colors <- c(
  "EUR"          = "#4a90d9",
  "SAS (EUR LD)" = "#e07b39",
  "SAS (CSA LD)" = "#8B4513"
)

sas_eur_mean <- round(mean(prs_SAS_eurLD_eurscaled, na.rm=TRUE), 3)
sas_csa_mean <- round(mean(prs_SAS_csaLD_eurscaled, na.rm=TRUE), 3)

p_dist <- ggplot(prs_plot_df,
                 aes(x = PRS, fill = Pop, color = Pop)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_manual(values  = pop_colors) +
  scale_color_manual(values = pop_colors) +
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = mean(prs_SAS_eurLD_eurscaled, na.rm=TRUE),
             linetype = "dashed", color = "#e07b39", linewidth = 0.7) +
  geom_vline(xintercept = mean(prs_SAS_csaLD_eurscaled, na.rm=TRUE),
             linetype = "dashed", color = "#8B4513", linewidth = 0.7) +
  annotate("text",
           x     = mean(prs_SAS_eurLD_eurscaled, na.rm=TRUE) + 0.05,
           y     = Inf, vjust = 1.5, hjust = 0, size = 3,
           color = "#e07b39",
           label = paste0("SAS EUR LD mean = ", sas_eur_mean, " SDs")) +
  annotate("text",
           x     = mean(prs_SAS_csaLD_eurscaled, na.rm=TRUE) + 0.05,
           y     = Inf, vjust = 3.2, hjust = 0, size = 3,
           color = "#8B4513",
           label = paste0("SAS CSA LD mean = ", sas_csa_mean, " SDs")) +
  labs(
    title    = "LDpred2-inf PRS Distribution: EUR vs SAS (EUR-standardized)",
    subtitle = "Infinitesimal model, EUR-derived vs CSA-derived effect sizes",
    x        = "Standardized PRS (EUR mean=0, SD=1)",
    y        = "Density",
    fill     = "Population",
    color    = "Population",
    caption  = paste0(
      "EUR n=", length(prs_EUR_ldpred2),
      "; SAS n=", length(prs_SAS_eurLD_ldpred2),
      "\nEstimated genome-wide h2 = ", round(h2_genome_EUR, 3),
      "\nC+T mean shift (EUR effects) = +1.3165 SDs"
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "gray40", size = 10),
    panel.grid.minor = element_blank()
  )

ggsave("results/figures/LDpred2_prs_distribution.png",
       plot = p_dist, width = 9, height = 5.5, dpi = 300)

# -----------------------------------------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------------------------------------
# Method: LDpred2-inf
# EUR and CSA chromosomes completed: 22/22
# EUR and CSA total SNPs: 1312946
# Estimated genome-wide h2: 0.3523
# EUR PRS mean (LDpred2-scaled): 0
# SAS EUR LD mean: 1.2688 
# SAS CSA LD mean: 1.2857 
# C+T mean shift: +1.3165 SDs
# LDpred2-inf mean shift: 1.2688 SDs

# Both methods confirm ~1.27-1.32 SD miscalibration