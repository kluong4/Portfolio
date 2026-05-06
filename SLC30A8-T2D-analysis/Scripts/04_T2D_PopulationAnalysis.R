# =============================================================================
# Script 04: Population-Level Expression Variation of T2D-Associated Genes
# Project:   Population and Evolutionary Analysis of SLC30A8 in Type 2 Diabetes
# Course:    BMI 5330: Introduction to Bioinformatics
#            UT Health Houston, Fall 2025
#
# Input:     Data/T2D_Data/20130606_sample_info.xlsx              (1000 Genomes sample metadata)
#            Data/T2D_Data/header.txt                             (GEUVADIS expression matrix header)
#            Data/T2D_Data/E-GEUV-1-query-results.fpkmss.tsv      (GEUVADIS RNA-seq expression data)
#            Data/T2D_Data/MONDO_0005148_associations_export.tsv  (GWAS Catalog T2D associations)
# Output:    Results/Top_T2D_Genes_Population_Expression_Variation.tsv
#
# Data sources:
#   GEUVADIS RNA-seq dataset: https://www.ebi.ac.uk/gxa/experiments/E-GEUV-1/Downloads
#   1000 Genomes sample info: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/
#   GWAS Catalog T2D associations (MONDO:0005148): https://www.ebi.ac.uk/gwas/efotraits/MONDO_0005148
# =============================================================================

# ----- Load Libraries --------------------------------------------------------
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)

# ----- Set Working Directory -------------------------------------------------
setwd("~/Desktop/Prev Courses/25Fall/BMI 5330/Project/SLC30A8-T2D-analysis")

# ----- Configuration ---------------------------------------------------------
POPULATIONS <- c("CEU", "FIN", "GBR", "TSI", "YRI")
GENE_OF_INTEREST <- "SLC30A8"   # change to plot a different gene


# =============================================================================
# STEP 1: Load Sample Metadata and Expression Matrix
# =============================================================================
sample_info <- read_excel("Data/T2D_Data/20130606_sample_info.xlsx", sheet = "Sample Info")

# Parse sample IDs from the header file
header   <- readLines("Data/T2D_Data/header.txt")
header   <- unlist(strsplit(header, "\t"))
header   <- header[-c(1, 2)] # remove first 2 non-sample columns
sampleID <- gsub(".*,\\s*", "", header) # extract clean sample IDs

# Match samples to population labels
matrix_header <- sample_info[sample_info$Sample %in% sampleID, ]

cat("Samples in expression matrix: ", length(sampleID), "\n")
cat("Samples matched to metadata:  ", nrow(matrix_header), "\n")

# Population breakdown:
print(table(matrix_header$Population))

# Load expression matrix and clean column names
expr_matrix <- read.delim(
  "Data/T2D_Data/E-GEUV-1-query-results.fpkmss.tsv",
  comment.char = "#",
  check.names  = FALSE
)
colnames(expr_matrix) <- gsub(".*,\\s*", "", colnames(expr_matrix))

cat("\nExpression matrix dimensions: ",
    nrow(expr_matrix), "genes x", ncol(expr_matrix), "columns\n")


# =============================================================================
# STEP 2: Compute Per-Population Mean Expression and Variance
# =============================================================================
pop_vector <- matrix_header$Population[match(sampleID, matrix_header$Sample)]
pop_groups <- split(sampleID, pop_vector)

# Mean expression per population per gene
pop_means <- sapply(pop_groups, function(samples) {
  rowMeans(expr_matrix[, samples], na.rm = TRUE)
})
pop_means <- cbind(expr_matrix[, c("Gene ID", "Gene Name")], pop_means)

# Keep only rows where all population columns are non-NaN
pop_means <- pop_means[
  apply(pop_means[, POPULATIONS], 1, function(x) all(!is.nan(x))), ]

# log2(FPKM + 1) transform
pop_means[, POPULATIONS] <- log2(pop_means[, POPULATIONS] + 1)

# Compute variance across populations
pop_means$Variance <- apply(pop_means[, POPULATIONS], 1, var, na.rm = TRUE)

cat("Genes with complete population data: ", nrow(pop_means), "\n")


# =============================================================================
# STEP 3: Filter to T2D-Associated Genes (GWAS Catalog)
# =============================================================================
gwas_df <- read.delim("Data/T2D_Data/MONDO_0005148_associations_export.tsv",
                      stringsAsFactors = FALSE)

# Extract unique mapped gene names from GWAS results
t2d_gene_list <- gwas_df %>%
  separate_rows(mappedGenes, sep = ",\\s*") %>%
  select(mappedGenes, pValue)
t2d_gene_list <- t2d_gene_list[!duplicated(t2d_gene_list$mappedGenes), ]

cat("Unique T2D GWAS genes: ", nrow(t2d_gene_list), "\n")

# Filter expression data to T2D genes
t2d_genes <- pop_means[pop_means$`Gene Name` %in% t2d_gene_list$mappedGenes, ]
t2d_genes$pValue <- t2d_gene_list$pValue[
  match(t2d_genes$`Gene Name`, t2d_gene_list$mappedGenes)]
t2d_genes <- t2d_genes[!duplicated(t2d_genes$`Gene Name`), ]

cat("T2D genes found in expression matrix: ", nrow(t2d_genes), "\n")


# =============================================================================
# STEP 4: Rank Genes by Population-Level Expression Variance
# =============================================================================
top_var_genes_pop <- t2d_genes %>%
  arrange(desc(Variance)) %>%
  select(`Gene ID`, `Gene Name`, pValue,
         all_of(POPULATIONS), Variance)

# Top 3 T2D genes by expression variance
head(top_var_genes_pop, 3)

# SLC30A8 Rank
print(top_var_genes_pop[top_var_genes_pop$`Gene Name` == GENE_OF_INTEREST, ])


# =============================================================================
# STEP 5: Save Ranked Gene Table
# =============================================================================
write.table(top_var_genes_pop,
            "Results/Top_T2D_Genes_Population_Expression_Variation.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# Boxplot — Gene of Interest Expression Across Populations
# =============================================================================
gene_expr <- expr_matrix %>%
  filter(`Gene Name` == GENE_OF_INTEREST) %>%
  select(-`Gene ID`, -`Gene Name`) %>%
  t() %>%
  as.data.frame()

colnames(gene_expr) <- "Expression"
gene_expr$Sample <- rownames(gene_expr)

# Add population labels
gene_expr$Population <- matrix_header$Population[
  match(gene_expr$Sample, matrix_header$Sample)]

gene_expr <- gene_expr %>% filter(!is.na(Population))

ggplot(gene_expr, aes(x = Population, y = Expression, fill = Population)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = paste0(GENE_OF_INTEREST, " Expression Across Populations"),
    x     = "Population",
    y     = "log2(FPKM + 1)"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# ggsave(paste0("Results/Fig3_", GENE_OF_INTEREST, "_expression_boxplot.png"), width = 8, height = 6, dpi = 300)