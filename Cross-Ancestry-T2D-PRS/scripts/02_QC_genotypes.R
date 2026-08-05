# =============================================================================
# Script 02: Genotype QC of 1000 Genomes Project Phase 3 Data
# Project:   Cross-ancestry T2D PRS Transferability (EUR -> CSA)
# Data:      1000 Genomes Project Phase 3 (bigsnpr format)
# =============================================================================
# OVERVIEW:
# This script loads the 1KGP Phase 3 genotype data downloaded in bigSNP
# format, subsets to EUR and SAS superpopulations, 
# performs sample- and SNP-level QC within each superpopulation
# separately, and saves QC-passed genotype data for use in PCA (Script 03)
# and PRS construction (Scripts 04-05).
#
# WHY SEPARATE QC PER SUPERPOPULATION:
# Running QC jointly across ancestries would confound MAF and HWE filters
# with population structure. For example, a SNP monomorphic in EUR but
# polymorphic in SAS would be incorrectly removed if MAF is computed
# jointly. Performing QC within each superpopulation ensures filters
# reflect within-population allele frequencies.
#
# OUTPUTS:
#   samples_EUR.txt           - EUR sample IDs passing QC
#   samples_SAS.txt           - SAS sample IDs passing QC
#   snps_pass_EUR.txt         - SNPs passing QC in EUR
#   snps_pass_SAS.txt         - SNPs passing QC in SAS
#   snps_pass_both.txt        - SNPs passing QC in BOTH (used downstream)
#   1KGP_EUR_QC.bed/bim/fam   - EUR PLINK files after QC
#   1KGP_SAS_QC.bed/bim/fam   - SAS PLINK files after QC
#   qc_summary.rds            - QC metrics for reporting
#   genotype_qc_summary.csv   - QC summary for manuscript
# =============================================================================

library(bigsnpr)
library(bigstatsr)
library(data.table)
library(ggplot2)

# Set PLINK2 command -- verified accessible from PATH
plink2 <- "plink2"

# -----------------------------------------------------------------------------
# STEP 1: Load 1KGP genotype data and population panel
# -----------------------------------------------------------------------------
# bigsnpr stores genotype data in two files:
#   .rds  -- R object containing map (variant info) and fam (sample info)
#   .bk   -- memory-mapped binary matrix of genotypes (too large for RAM)
# The attach_bigSNP() function links to the .bk file without loading it
# fully into RAM -- essential for a 4 GB file on a 24 GB machine.

# Loading 1KGP genotype data
obj.bigsnp <- snp_attach("data/1KGP/1000G_phase3_common_norel.rds")

# Extract components for easy access
G   <- obj.bigsnp$genotypes  # genotype matrix (memory-mapped)
map <- obj.bigsnp$map         # variant info (chr, pos, ref, alt)
fam <- obj.bigsnp$fam         # sample info (family ID, sample ID)

cat("Variants in dataset:", nrow(map), "\n")      # 1664852
cat("Samples in dataset: ", nrow(fam), "\n\n")    # 2490 

# Load population panel
pop_panel <- fread("data/1KGP/integrated_call_samples_v3.20130502.ALL.panel")

# Verify sample IDs match between fam file and population panel
cat("Samples in fam file:   ", nrow(fam), "\n")
cat("Samples in pop panel:  ", nrow(pop_panel), "\n")
cat("Matching sample IDs:   ", 
    sum(fam$sample.ID %in% pop_panel$sample), "\n\n")

# -----------------------------------------------------------------------------
# Identify EUR and SAS samples
# -----------------------------------------------------------------------------
# We subset to EUR (n=503) and SAS (n=489) superpopulations only.

idx_EUR <- which(fam$sample.ID %in% 
                   pop_panel[super_pop == "EUR", sample])
idx_SAS <- which(fam$sample.ID %in% 
                   pop_panel[super_pop == "SAS", sample])

cat("EUR samples identified:", length(idx_EUR), "\n")
cat("SAS samples identified:", length(idx_SAS), "\n\n")

# Save sample lists for PLINK2 subsetting
eur_samples <- data.table(
  FID = fam$family.ID[idx_EUR],
  IID = fam$sample.ID[idx_EUR]
)
sas_samples <- data.table(
  FID = fam$family.ID[idx_SAS],
  IID = fam$sample.ID[idx_SAS]
)

fwrite(eur_samples, "data/1KGP/outputs/samples_EUR.txt", sep = "\t", col.names = TRUE)
fwrite(sas_samples, "data/1KGP/outputs/samples_SAS.txt", sep = "\t", col.names = TRUE)

# -----------------------------------------------------------------------------
# SNP-level QC per superpopulation using PLINK2
# -----------------------------------------------------------------------------
# Filters applied:
#   --maf 0.01        : Remove SNPs with MAF < 1% within superpopulation
#                       Very rare variants have unreliable effect estimates
#                       and contribute noise to PRS
#   --geno 0.05       : Remove SNPs with >5% missing genotypes
#                       High missingness suggests genotyping failure
#   --hwe 1e-6        : Remove SNPs failing Hardy-Weinberg Equilibrium test
#                       (p < 1e-6). Extreme HWE deviation suggests
#                       genotyping error. Threshold is lenient to avoid
#                       removing truly selected loci.
#
# QC is run separately per superpopulation so that MAF and HWE are
# computed within the relevant ancestry group.
for(pop in c("EUR", "SAS")) {
  
  cat("Processing", pop, "...\n")
  
  system(paste(
    plink2,
    "--bfile",  "data/1KGP/1000G_phase3_common_norel",
    "--keep",   paste0("data/1KGP/samples_", pop, ".txt"),
    "--maf",    "0.01",
    "--geno",   "0.05",
    "--hwe",    "1e-6",
    "--write-snplist allow-dups",    # added allow-dups to handle duplicate IDs
    "--no-psam-pheno",
    "--out",    paste0("data/1KGP/snps_pass_", pop)
  ))
  
  # Count SNPs passing QC
  snps_pass <- fread(
    paste0("data/1KGP/snps_pass_", pop, ".snplist"),
    header = FALSE
  )
  cat(pop, "SNPs passing QC:", nrow(snps_pass), "\n\n")
}

# -----------------------------------------------------------------------------
# Find SNPs passing QC in BOTH superpopulations
# -----------------------------------------------------------------------------
# Downstream PRS construction requires the same SNP set across both
# populations so scores are directly comparable. We take the intersection
# of SNPs passing QC in both EUR and SAS. This is the SNP set used in
# all subsequent scripts.

snps_EUR <- fread("data/1KGP/snps_pass_EUR.snplist", header = FALSE)$V1
snps_SAS <- fread("data/1KGP/snps_pass_SAS.snplist", header = FALSE)$V1

snps_both <- intersect(snps_EUR, snps_SAS)
cat("SNPs passing QC in EUR:        ", length(snps_EUR), "\n")
cat("SNPs passing QC in SAS:        ", length(snps_SAS), "\n")
cat("SNPs passing QC in BOTH:       ", length(snps_both), "\n\n")

# Save intersection SNP list
fwrite(data.table(SNP = snps_both), 
       "data/1KGP/snps_pass_both.txt",
       col.names = FALSE)

# -----------------------------------------------------------------------------
# STEP 5: Sample-level QC per superpopulation using PLINK2
# -----------------------------------------------------------------------------
# Filters applied:
#   --mind 0.02       : Remove samples with >2% missing genotypes
#                       High per-sample missingness suggests poor DNA
#                       quality or failed genotyping
#
# Note: Relatedness QC is not performed here because 1KGP already provides
# a "common_norel" dataset with related individuals removed (indicated by
# the filename). Heterozygosity outlier removal is performed in Script 03
# after PCA, where population structure is already accounted for.
for(pop in c("EUR", "SAS")) {
  
  cat("Processing", pop, "...\n")
  
  system(paste(
    plink2,
    "--bfile",      "data/1KGP/1000G_phase3_common_norel",
    "--keep",       paste0("data/1KGP/samples_", pop, ".txt"),
    "--extract",    "data/1KGP/snps_pass_both.txt",
    "--mind",       "0.02",
    "--make-bed",
    "--no-psam-pheno",
    "--out",        paste0("data/1KGP/1KGP_", pop, "_QC")
  ))
  
  # Report samples retained
  fam_qc <- fread(paste0("data/1KGP/1KGP_", pop, "_QC.fam"))
  cat(pop, "samples after QC:", nrow(fam_qc), "\n\n")
}

# -----------------------------------------------------------------------------
# QC Summary
# -----------------------------------------------------------------------------
# Collect all QC metrics into a summary table for manuscript reporting.
# These numbers go into the Methods/Results section under
# "Genotype Quality Control."

fam_EUR_qc <- fread("data/1KGP/1KGP_EUR_QC.fam")
fam_SAS_qc <- fread("data/1KGP/1KGP_SAS_QC.fam")
bim_qc     <- fread("data/1KGP/1KGP_EUR_QC.bim")  # same SNPs in both

cat("EUR samples passing QC:", nrow(fam_EUR_qc), "\n")
cat("SAS samples passing QC:", nrow(fam_SAS_qc), "\n")
cat("SNPs passing QC (both):", length(snps_both), "\n")

qc_summary <- data.table(
  Population         = c("EUR", "SAS"),
  Samples_before_QC  = c(length(idx_EUR), length(idx_SAS)),
  Samples_after_QC   = c(nrow(fam_EUR_qc), nrow(fam_SAS_qc)),
  SNPs_before_QC     = nrow(map),
  SNPs_EUR_pass      = length(snps_EUR),
  SNPs_SAS_pass      = length(snps_SAS),
  SNPs_both_pass     = length(snps_both)
)

# Save QC summary
print(qc_summary)
saveRDS(qc_summary, "data/1KGP/qc_summary.rds")
fwrite(qc_summary,  "results/QC/genotype_qc_summary.csv")