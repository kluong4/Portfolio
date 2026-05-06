# =============================================================================
# Script 02: Population-Level Variation Analysis of SLC30A8
# Project:   Population and Evolutionary Analysis of SLC30A8 in Type 2 Diabetes
# Course:    BMI 5330: Introduction to Bioinformatics
#            UT Health Houston, Fall 2025
#
# Input:     Data/SLC30A8_Data_SLC30A8_variants.bed    (from 01_extract_SLC30A8_variants.sh)
# Output:    Data/ProcessedData_RDS/variants.rds
#            Results/SLC30A8_variants.tsv
#            Results/Fig1_allele_frequency_boxplot.png
# =============================================================================

# ----- Load Libraries --------------------------------------------------------
library(stringr)
library(ggplot2)
library(tidyr)
library(dplyr)

# ----- Set Working Directory -------------------------------------------------
setwd("~/Desktop/Prev Courses/25Fall/BMI 5330/Project/SLC30A8-T2D-analysis")


# =============================================================================
# STEP 1: Load Variants and Extract Allele Frequencies
# =============================================================================
variants <- read.table("Data/SLC30A8_Data/SLC30A8_variants.bed", header = FALSE)
colnames(variants) <- c("chr", "start", "end", "id", "info", "dot", "strand")

# Helper function: extract a named field from the INFO column
extract_AF_field <- function(info, field) {
  pattern <- paste0(field, "=([0-9.]+)")
  as.numeric(str_match(info, pattern)[, 2])
}

# Overall allele frequency + per-population frequencies
variants$AF <- sapply(variants$info, extract_AF_field, field = "AF")
variants$ASN_AF <- sapply(variants$info, extract_AF_field, field = "ASN_AF")
variants$AFR_AF <- sapply(variants$info, extract_AF_field, field = "AFR_AF")
variants$AMR_AF <- sapply(variants$info, extract_AF_field, field = "AMR_AF")
variants$EUR_AF <- sapply(variants$info, extract_AF_field, field = "EUR_AF")

# =============================================================================
# STEP 2: Summary Statistics Per Population
# =============================================================================
variants_population <- pivot_longer(
  variants,
  cols = c(ASN_AF, AFR_AF, AMR_AF, EUR_AF),
  names_to = "Population",
  values_to = "AF_pop"
)
variants_population <- variants_population[!is.na(variants_population$AF_pop), ]
variants_population$Population <- gsub("_AF$", "", variants_population$Population)

# Average allele frequency per population
pop_summary <- variants_population %>%
  group_by(Population) %>%
  summarize(Avg_AF = mean(AF_pop, na.rm = TRUE), .groups = "drop")
print(pop_summary)

# =============================================================================
# STEP 3: Boxplot — Allele Frequency Distribution Across Populations
# =============================================================================
ggplot(variants_population, aes(x = Population, y = AF_pop, fill = Population)) +
  geom_boxplot(alpha = 0.6) +
  scale_y_log10() +
  labs(
    title = "Boxplot of Allele Frequency Distributions Across Populations",
    x = "Population",
    y = "Allele Frequency (log10 scale)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave("Results/Fig1_AF_boxplot.png", width = 8, height = 6, dpi = 300)


# =============================================================================
# STEP 4: Extract Additional INFO Fields (VT, AA) and Clean Up
# =============================================================================
variants <- variants %>%
  mutate(
    VT = str_extract(info, "(?<=VT=)[^;]+"),
    AA = str_extract(info, "(?<=AA=)[^;]+")
  ) %>%
  select(-info) %>%
  select(id, chr, start, end, dot, strand,
         VT, AA,
         AF, AFR_AF, AMR_AF, ASN_AF, EUR_AF)

# Variant type breakdown
print(table(variants$VT))


# =============================================================================
# STEP 5: Save Outputs
# =============================================================================
# Save cleaned variants table
write.table(variants, "Results/SLC30A8_variants.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Save RDS for use in Script 03
saveRDS(variants, "Data/ProcessedData_RDS/variants.rds")
