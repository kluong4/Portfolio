# =============================================================================
# Script 04: Clumping and Thresholding (C+T) PRS Construction
# Project:   Cross-ancestry T2D PRS Transferability (EUR -> CSA)
# Data:      Pan-UKB E11 summary statistics + 1KGP EUR/SAS genotypes
# =============================================================================
# OVERVIEW:
# This script constructs polygenic risk scores using the Clumping and
# Thresholding (C+T) method. PLINK2 is used for clumping (memory efficient,
# reads from disk) and scoring. bigsnpr::snp_clumping() was not used because
# it requires loading full LD matrices into RAM (~15-20GB), exceeding local
# machine capacity.
#
# KEY FINDING: EUR-derived PRS applied to SAS individuals with EUR
# standardization shows a mean shift of +1.32 SDs, indicating substantial
# score miscalibration that would affect clinical deployment.
#
# FOUR PRS SCORES CONSTRUCTED:
#   1. EUR effects + EUR LD → EUR samples (reference, mean=0 by definition)
#   2. EUR effects + EUR LD → SAS samples (clinical simulation)
#   3. EUR effects + SAS LD → SAS samples (isolates LD reference effect)
#   4. CSA effects + SAS LD → SAS samples (fully ancestry-matched)
#
# OUTPUTS:
#   data/1KGP/1KGP_{EUR,SAS}_QC_varid.{bed,bim,fam} - varid-updated genotypes
#   data/1KGP/clump_{scenario}.clumps               - clumped SNP lists
#   data/sumstats/scores/score_{scenario}_thr*.txt  - PLINK2 score files
#   results/PRS/CT/prs_CT_{scenario}.txt            - primary PRS output
#   results/PRS/CT/all_prs_CT.rds                   - all thresholds, all scenarios
#   results/PRS/CT/ct_clump_summary.csv             - SNP counts per threshold
#   results/figures/CT_prs_distribution_preview.png - diagnostic figure
# =============================================================================

library(bigsnpr)
library(bigstatsr)
library(data.table)
library(ggplot2)

plink2 <- "plink2"

# -----------------------------------------------------------------------------
# STEP 1: Load genotype data and summary statistics
# -----------------------------------------------------------------------------

# Loading genotype data
obj.bigsnp <- snp_attach("data/1KGP/1000G_phase3_common_norel.rds")
G   <- obj.bigsnp$genotypes
map <- obj.bigsnp$map
fam <- obj.bigsnp$fam

# Load final QC-passed sample lists (post-PCA outlier removal)
samples_EUR <- fread("data/1KGP/samples_EUR_final.txt")
samples_SAS <- fread("data/1KGP/samples_SAS_final.txt")

idx_EUR <- which(fam$sample.ID %in% samples_EUR$IID) #EUR samples: 503 
idx_SAS <- which(fam$sample.ID %in% samples_SAS$IID) #SAS samples: 478

# Load formatted summary statistics (11239095 SNPs each)
# chr, pos, ref, alt, varID, beta, se, neglog10_pval, af_cases, af_controls
sumstats_EUR <- readRDS("data/sumstats/outputs/sumstats_EUR_formatted.rds")
sumstats_CSA <- readRDS("data/sumstats/outputs/sumstats_CSA_formatted.rds")

# -----------------------------------------------------------------------------
# STEP 2: Build sumstats files for PLINK2
# -----------------------------------------------------------------------------
# PLINK2 --clump requires SNP and P columns.
# PLINK2 --score requires SNP, effect allele, and beta columns.
# SNP IDs must match bim file variant IDs exactly (chr:pos:ref:alt format).
# P values must be raw p-values (not -log10 transformed).
# We build a single file per ancestry with all needed columns.

sumstats_EUR_plink <- data.table(
  SNP  = paste(sumstats_EUR$chr, sumstats_EUR$pos,
               sumstats_EUR$ref, sumstats_EUR$alt, sep = ":"),
  P    = 10^(-sumstats_EUR$neglog10_pval),  # convert -log10(p) to raw p
  BETA = sumstats_EUR$beta,
  SE   = sumstats_EUR$se,
  A1   = sumstats_EUR$alt                   # effect allele = alternate allele
)

sumstats_CSA_plink <- data.table(
  SNP  = paste(sumstats_CSA$chr, sumstats_CSA$pos,
               sumstats_CSA$ref, sumstats_CSA$alt, sep = ":"),
  P    = 10^(-sumstats_CSA$neglog10_pval),
  BETA = sumstats_CSA$beta,
  SE   = sumstats_CSA$se,
  A1   = sumstats_CSA$alt
)

fwrite(sumstats_EUR_plink, "data/sumstats/outputs/sumstats_EUR_forclump.txt", sep = "\t")
fwrite(sumstats_CSA_plink, "data/sumstats/outputs/sumstats_CSA_forclump.txt", sep = "\t")

# -----------------------------------------------------------------------------
# STEP 3: Update bim variant IDs to chr:pos:ref:alt format
# -----------------------------------------------------------------------------
# 1KGP bim files use rsIDs which cannot be matched to our chr:pos:ref:alt
# sumstats format. PLINK2 --set-all-var-ids updates IDs in place.
# Format string @:#:\$r:\$a = chr:pos:ref:alt
# --new-id-max-allele-len 1000 handles multi-character indel alleles.
#
# Verified: 1,312,855 / 1,422,959 variants match (92%) after ID update.

for(pop in c("EUR", "SAS")) {
  
  out_file <- paste0("data/1KGP/1KGP_", pop, "_QC_varid")
  
  # Skip if already done (avoid re-running unnecessarily)
  if(file.exists(paste0(out_file, ".bed"))) {
    cat(pop, "varid bim already exists -- skipping.\n")
    next
  }
  
  system(paste(
    plink2,
    "--bfile",               paste0("data/1KGP/1KGP_", pop, "_QC"),
    "--set-all-var-ids",     "@:#:\\$r:\\$a",
    "--new-id-max-allele-len", "1000",
    "--make-bed",
    "--no-psam-pheno",
    "--out",                 out_file
  ))
  
  bim_check <- fread(paste0(out_file, ".bim"))
  colnames(bim_check) <- c("chr","snpID","cm","pos","a1","a2")
  cat(pop, "bim ID example:", bim_check$snpID[1], "\n")
}

# -----------------------------------------------------------------------------
# STEP 4: Create --keep files compatible with PLINK2 fam format
# -----------------------------------------------------------------------------
# The varid-updated fam files have FID=0 (integer) but our sample lists
# have FID=sample_ID (character). PLINK2 --keep with IID-only format
# (single column, no header) avoids the FID mismatch.
# Verified: 478/484 SAS and 502/503 EUR IIDs match the fam file.

fwrite(data.table(IID = samples_EUR$IID), "data/1KGP/keep_EUR_final.txt", col.names = FALSE)
fwrite(data.table(IID = samples_SAS$IID), "data/1KGP/keep_SAS_final.txt", col.names = FALSE)

# -----------------------------------------------------------------------------
# STEP 5: PLINK2 clumping
# -----------------------------------------------------------------------------
# PLINK2 --clump parameters:
#   --clump-p1 1  : all SNPs eligible as index SNPs
#   --clump-p2 1  : all SNPs eligible as clump members
#   --clump-r2 0.1: remove SNPs with r2 > 0.1 with index SNP
#   --clump-kb 250: 250kb window around each index SNP
#
# Running with p-threshold = 1 retains all SNPs so we can apply
# different thresholds in post-processing without re-clumping.
# PLINK2 uses ~2-4GB RAM vs bigsnpr's 15-20GB, making it feasible
# on a 24GB local machine.

run_plink_clumping <- function(pop_bfile, sumstats_file, label) {
  
  cat("Clumping:", label, "\n")
  out_prefix <- paste0("data/1KGP/clump_", label)
  
  # Skip if already done
  if(file.exists(paste0(out_prefix, ".clumps"))) {
    cat("  Clump file already exists -- loading from disk.\n")
    clumped <- fread(paste0(out_prefix, ".clumps"))
    cat("  SNPs loaded:", nrow(clumped), "\n\n")
    return(clumped)
  }
  
  system(paste(
    plink2,
    "--bfile",      pop_bfile,
    "--clump",      sumstats_file,
    "--clump-p1",   "1",
    "--clump-p2",   "1",
    "--clump-r2",   "0.1",
    "--clump-kb",   "250",
    "--no-psam-pheno",
    "--out",        out_prefix
  ))
  
  clump_file <- paste0(out_prefix, ".clumps")
  if(!file.exists(clump_file)) {
    cat("WARNING: Clump file not created. Check",
        out_prefix, ".log\n")
    return(NULL)
  }
  
  clumped <- fread(clump_file)
  cat("  SNPs after clumping:", nrow(clumped), "\n\n")
  return(clumped)
}

# SCENARIO A: EUR effects, EUR LD
clump_A <- run_plink_clumping(
  "data/1KGP/1KGP_EUR_QC_varid",
  "data/sumstats/outputs/sumstats_EUR_forclump.txt",
  "EUR_effects_EUR_LD"
)

# SCENARIO B: EUR effects, SAS LD
clump_B <- run_plink_clumping(
  "data/1KGP/1KGP_SAS_QC_varid",
  "data/sumstats/outputs/sumstats_EUR_forclump.txt",
  "EUR_effects_SAS_LD"
)

# SCENARIO C: CSA effects, SAS LD
clump_C <- run_plink_clumping(
  "data/1KGP/1KGP_SAS_QC_varid",
  "data/sumstats/outputs/sumstats_CSA_forclump.txt",
  "CSA_effects_SAS_LD"
)

# -----------------------------------------------------------------------------
# STEP 6: Match clumped SNPs back to sumstats to retrieve effect sizes
# -----------------------------------------------------------------------------
# The clump output contains index SNP IDs only. We match these back to
# the full sumstats to retrieve BETA and A1 columns needed for scoring.
# We also count SNPs surviving each p-value threshold after clumping.

thresholds_p <- c(5e-8, 1e-6, 1e-4, 0.001, 0.01, 0.05, 0.1, 0.5, 1.0)

prepare_score_file <- function(clumped, sumstats_plink, label) {
  
  if(is.null(clumped)) return(NULL)
  
  # Match clumped SNP IDs to sumstats
  # clumped$ID contains the index SNP IDs from PLINK2 --clump output
  matched <- sumstats_plink[SNP %in% clumped$ID]
  cat(label, "-- clumped SNPs matched to sumstats:", nrow(matched), "\n")
  cat("  Columns:", paste(colnames(matched), collapse = ", "), "\n")
  
  # Count SNPs per p-value threshold
  # This tells us how many SNPs contribute to PRS at each threshold
  snp_counts <- sapply(thresholds_p, function(thr) {
    sum(matched$P <= thr, na.rm = TRUE)
  })
  
  count_df <- data.table(
    label       = label,
    threshold_p = thresholds_p,
    n_snps      = snp_counts
  )
  
  cat("  SNPs per threshold:\n")
  print(count_df[, .(threshold_p, n_snps)])
  cat("\n")
  
  return(list(matched = matched, counts = count_df))
}

score_A <- prepare_score_file(clump_A, sumstats_EUR_plink, "EUR_effects_EUR_LD")
score_B <- prepare_score_file(clump_B, sumstats_EUR_plink, "EUR_effects_SAS_LD")
score_C <- prepare_score_file(clump_C, sumstats_CSA_plink, "CSA_effects_SAS_LD")

# Save clumping summary table for manuscript Table 1
clump_summary <- rbind(score_A$counts, score_B$counts, score_C$counts)
fwrite(clump_summary, "results/PRS/CT/ct_clump_summary.csv")

# -----------------------------------------------------------------------------
# STEP 7: Compute PRS with PLINK2 --score
# -----------------------------------------------------------------------------
# PLINK2 --score computes: PRS_i = sum(beta_j * dosage_ij) / n_snps
# Score file format (no header, tab-separated):
#   col 1: SNP ID, col 2: effect allele, col 3: beta
# no-mean-imputation: missing genotypes contribute 0 rather than
# the mean dosage, which is more conservative for PRS computation.
#
# A column verification step confirms score files have exactly 3 columns
# before calling PLINK2 to catch write errors early.

compute_prs_plink2 <- function(score_data, bfile, keep_file,
                               threshold_p, scenario_label, pop_label) {
  
  # Filter to SNPs passing p-value threshold
  snps_thr <- score_data$matched[P <= threshold_p]
  
  if(nrow(snps_thr) == 0) {
    cat("No SNPs pass threshold", threshold_p,
        "for", scenario_label, "\n")
    return(NULL)
  }
  
  # Build consistent file name strings
  # formatC handles scientific notation; gsub cleans for filenames
  thr_str    <- gsub("\\.", "_", gsub("e-0", "e-",
                                      formatC(threshold_p, format = "g")))
  score_file <- paste0("data/sumstats/scores/score_",
                       scenario_label, "_thr", thr_str, ".txt")
  out_prefix <- paste0("results/PRS/CT/prs_",
                       scenario_label, "_", pop_label,
                       "_thr", thr_str)
  
  # Write score file: 3 columns, tab-separated, no header, no quotes
  score_out <- data.table(
    SNP  = snps_thr$SNP,
    A1   = snps_thr$A1,
    BETA = snps_thr$BETA
  )
  fwrite(score_out, score_file,
         sep = "\t", col.names = FALSE, quote = FALSE)
  
  # Verify 3 columns before calling PLINK2
  check <- fread(score_file, header = FALSE, nrows = 3)
  if(ncol(check) != 3) {
    cat("ERROR: Score file has", ncol(check),
        "columns (expected 3) for", scenario_label, "\n")
    return(NULL)
  }
  
  # Run PLINK2 --score
  cmd <- paste(
    plink2,
    "--bfile",       bfile,
    "--keep",        keep_file,
    "--score",       score_file, "1 2 3 no-mean-imputation",
    "--no-psam-pheno",
    "--out",         out_prefix
  )
  system(cmd)
  
  # Read .sscore output file
  sscore_file <- paste0(out_prefix, ".sscore")
  if(!file.exists(sscore_file)) {
    log_file  <- paste0(out_prefix, ".log")
    if(file.exists(log_file)) {
      log_lines <- readLines(log_file)
      err_lines <- grep("Error|error", log_lines, value = TRUE)
      cat("PLINK2 error for", scenario_label, "thr", threshold_p, ":\n")
      cat(err_lines, sep = "\n")
    }
    return(NULL)
  }
  
  prs_out <- fread(sscore_file)
  cat("Scored:", nrow(prs_out), "samples |",
      scenario_label, "| thr:", threshold_p,
      "| SNPs:", nrow(snps_thr), "\n")
  return(prs_out)
}

# Compute PRS at all thresholds for all 4 scenarios
all_prs <- list()

for(thr in thresholds_p) {
  
  thr_str <- gsub("\\.", "_", gsub("e-0", "e-",
                                   formatC(thr, format = "g")))
  cat("--- Threshold p <", thr, "---\n")
  
  # Score 1: EUR effects, EUR LD → EUR samples
  # This is the reference -- by definition mean=0 after EUR standardization
  prs1 <- compute_prs_plink2(
    score_A,
    "data/1KGP/1KGP_EUR_QC_varid",
    "data/1KGP/keep_EUR_final.txt",
    thr, "EUR_effects_EUR_LD", "EUR"
  )
  
  # Score 2: EUR effects, EUR LD → SAS samples
  # Clinical simulation: EUR PRS applied universally to SAS individuals
  # Mean shift from 0 quantifies miscalibration
  prs2 <- compute_prs_plink2(
    score_A,
    "data/1KGP/1KGP_SAS_QC_varid",
    "data/1KGP/keep_SAS_final.txt",
    thr, "EUR_effects_EUR_LD", "SAS"
  )
  
  # Score 3: EUR effects, SAS LD → SAS samples
  # Isolates contribution of LD reference ancestry to miscalibration
  # Compare to Score 2: difference = LD reference effect
  prs3 <- compute_prs_plink2(
    score_B,
    "data/1KGP/1KGP_SAS_QC_varid",
    "data/1KGP/keep_SAS_final.txt",
    thr, "EUR_effects_SAS_LD", "SAS"
  )
  
  # Score 4: CSA effects, SAS LD → SAS samples
  # Fully ancestry-matched PRS
  # Compare to Score 2: total improvement from ancestry matching
  prs4 <- compute_prs_plink2(
    score_C,
    "data/1KGP/1KGP_SAS_QC_varid",
    "data/1KGP/keep_SAS_final.txt",
    thr, "CSA_effects_SAS_LD", "SAS"
  )
  
  all_prs[[thr_str]] <- list(
    EUR_in_EUR      = prs1,
    EUR_in_SAS      = prs2,
    EURLDsas_in_SAS = prs3,
    CSA_in_SAS      = prs4
  )
}

# Save all PRS results (all thresholds, all scenarios)
saveRDS(all_prs, "results/PRS/CT/all_prs_CT.rds")

# -----------------------------------------------------------------------------
# STEP 8: Standardize and save primary PRS scores
# -----------------------------------------------------------------------------
# Primary analysis uses p < 0.05 threshold (balance of signal and noise).
# Two standardization approaches:
#   within-population: mean=0, sd=1 within each group separately
#   EUR-parameter:     subtract EUR mean, divide by EUR sd
#                      simulates clinical deployment using EUR reference

primary_thr <- gsub("\\.", "_", gsub("e-0", "e-", formatC(0.05, format = "g")))

prs_EUR_raw       <- all_prs[[primary_thr]]$EUR_in_EUR$SCORE1_AVG
prs_SAS_raw       <- all_prs[[primary_thr]]$EUR_in_SAS$SCORE1_AVG
prs_SAS_sasld_raw <- all_prs[[primary_thr]]$EURLDsas_in_SAS$SCORE1_AVG
prs_SAS_csa_raw   <- all_prs[[primary_thr]]$CSA_in_SAS$SCORE1_AVG

# EUR standardization parameters (used for clinical simulation)
eur_mean <- mean(prs_EUR_raw, na.rm = TRUE)
eur_sd   <- sd(prs_EUR_raw,   na.rm = TRUE)

prs_EUR_eurscaled    <- (prs_EUR_raw       - eur_mean) / eur_sd
prs_SAS_eurscaled    <- (prs_SAS_raw       - eur_mean) / eur_sd
prs_SAS_sasld_scaled <- (prs_SAS_sasld_raw - eur_mean) / eur_sd
prs_SAS_csa_scaled   <- (prs_SAS_csa_raw   - eur_mean) / eur_sd

# PRIMARY PRS SUMMARY (p < 0.05)
cat("EUR PRS mean (EUR-scaled):        ", round(mean(prs_EUR_eurscaled,    na.rm = TRUE), 4), "\n") # 0
cat("SAS PRS mean EUR effects EUR LD:  ", round(mean(prs_SAS_eurscaled,    na.rm = TRUE), 4), "\n") # 1.3165
cat("SAS PRS mean EUR effects SAS LD:  ", round(mean(prs_SAS_sasld_scaled, na.rm = TRUE), 4), "\n") # 1.1866
cat("SAS PRS mean CSA effects SAS LD:  ", round(mean(prs_SAS_csa_scaled,   na.rm = TRUE), 4), "\n") # 5.5477
cat("Mean shift EUR→SAS (EUR effects): ",
    round(mean(prs_SAS_eurscaled, na.rm = TRUE) -
            mean(prs_EUR_eurscaled, na.rm = TRUE), 4), "\n\n") # 1.3165

# Save flat files for each scenario (used by Scripts 06 and 10)
scenarios <- list(
  list(name = "EUR_in_EUR", prs  = all_prs[[primary_thr]]$EUR_in_EUR),
  list(name = "EUR_in_SAS", prs  = all_prs[[primary_thr]]$EUR_in_SAS),
  list(name = "EURLDsas_in_SAS", prs  = all_prs[[primary_thr]]$EURLDsas_in_SAS),
  list(name = "CSA_in_SAS", prs  = all_prs[[primary_thr]]$CSA_in_SAS)
)

for(sc in scenarios) {
  if(is.null(sc$prs)) next
  out <- data.table(
    FID            = sc$prs$`#FID`,
    IID            = sc$prs$IID,
    PRS_raw        = sc$prs$SCORE1_AVG,
    PRS_std_within = as.numeric(scale(sc$prs$SCORE1_AVG)),
    PRS_std_EUR    = (sc$prs$SCORE1_AVG - eur_mean) / eur_sd
  )
  outfile <- paste0("results/PRS/CT/prs_CT_", sc$name, ".txt")
  fwrite(out, outfile, sep = "\t")
  cat("Saved:", outfile, "\n")
}

# Save standardization parameters for use in Script 06
saveRDS(list(eur_mean = eur_mean, eur_sd = eur_sd), "results/PRS/CT/eur_standardization_params.rds")

# -----------------------------------------------------------------------------
# STEP 9: Diagnostic figure
# -----------------------------------------------------------------------------
# Preview of the main PRS distribution figure (Script 10 produces the
# final version with all 4 scenarios). This figure shows the key finding:
# SAS mean shift of +1.32 SDs under EUR standardization indicating
# substantial PRS miscalibration when EUR reference is applied universally.

prs_plot_df <- data.table(
  PRS  = c(prs_EUR_eurscaled, prs_SAS_eurscaled),
  Pop  = c(rep("EUR", length(prs_EUR_eurscaled)),
           rep("SAS", length(prs_SAS_eurscaled)))
)

sas_mean_label <- round(mean(prs_SAS_eurscaled, na.rm = TRUE), 3)

p_dist <- ggplot(prs_plot_df, aes(x = PRS, fill = Pop, color = Pop)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  scale_fill_manual(values  = c("EUR" = "#4a90d9", "SAS" = "#e07b39")) +
  scale_color_manual(values = c("EUR" = "#4a90d9", "SAS" = "#e07b39")) +
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = mean(prs_SAS_eurscaled, na.rm = TRUE),
             linetype = "dashed", color = "#e07b39", linewidth = 0.8) +
  annotate("text",
           x     = mean(prs_SAS_eurscaled, na.rm = TRUE) + 0.05,
           y     = Inf, vjust = 1.5, hjust = 0, size = 3.2,
           color = "#e07b39",
           label = paste0("SAS mean = ", sas_mean_label, " SDs")) +
  labs(
    title    = "C+T PRS Distribution: EUR vs SAS (EUR-standardized)",
    subtitle = "EUR effect sizes, EUR LD clumping, p < 0.05 threshold",
    x        = "Standardized PRS (EUR mean=0, SD=1)",
    y        = "Density",
    caption  = paste0("EUR n=", sum(!is.na(prs_EUR_eurscaled)),
                      "; SAS n=", sum(!is.na(prs_SAS_eurscaled)),
                      "\nMean shift = +", sas_mean_label, " SDs")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "gray40", size = 10),
    panel.grid.minor = element_blank()
  )

ggsave("results/PRS/CT/CT_prs_distribution_preview.png", plot = p_dist, width = 8, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------------------------------------
cat("Clumped SNPs (EUR/EUR LD):        ", nrow(clump_A), "\n") #158487
cat("Clumped SNPs (EUR/SAS LD):        ", nrow(clump_B), "\n") #174047
cat("Clumped SNPs (CSA/SAS LD):        ", nrow(clump_C), "\n") #174401
cat("Primary threshold:                 p < 0.05\n")           #p < 0.05
cat("EUR PRS mean (EUR-scaled):        ",
    round(mean(prs_EUR_eurscaled,    na.rm = TRUE), 4), "\n")  #0
cat("SAS PRS mean EUR effects EUR LD:  ",
    round(mean(prs_SAS_eurscaled,    na.rm = TRUE), 4), "\n")  #1.3165
cat("SAS PRS mean EUR effects SAS LD:  ",
    round(mean(prs_SAS_sasld_scaled, na.rm = TRUE), 4), "\n")  #1.1866
cat("SAS PRS mean CSA effects SAS LD:  ",
    round(mean(prs_SAS_csa_scaled,   na.rm = TRUE), 4), "\n")  #5.5477
cat("KEY FINDING -- Mean shift (SAS-EUR):",
    round(mean(prs_SAS_eurscaled, na.rm = TRUE) -
            mean(prs_EUR_eurscaled, na.rm = TRUE), 4), "SDs\n") #1.3165 SDs