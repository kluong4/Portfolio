# =============================================================================
# Script 01: Cleaning Data
# Project:   TCGA-LIHC DEG Analysis
# Course:    INB 321G: Computational Biology
#            UT Austin, Spring 2025
# Input:     Data/gdc_sample_sheet_LIHC_StageII.csv
#            Data/gdc_sample_sheet_LIHC_StageII.csv
# Output:    Data/CleanGroup1.csv
#            Data/CleanGroup2.csv
# =============================================================================

# ----- Load Libraries --------------------------------------------------------
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(biomaRt)
library(ggplot2)  
library(ggrepel)  
library(scatterD3)

# ----- Set Working Directory -------------------------------------------------
setwd("/Volumes/DONGHOA/Kim's UT 25 BIO Undergrad/4 Spring 2025/INB 321G/project/TCGA-LIHC-DEGanalysis")

# ----- Load Unclean Data -----------------------------------------------------
# Cases pre-selected in GDC web portal
group1 <- read.csv("Data/gdc_sample_sheet_LIHC_StageII.csv")
group2 <- read.csv("Data/gdc_sample_sheet_LIHC_StageIII.csv")

# ----- Initial GDC Query -----------------------------------------------------
query_LIHC <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access = "open"
)

# Extract metadata table from query results
samDf <- query_LIHC$results[[1]]

# ----- Match Samples to Selected Groups --------------------------------------
# Identify which queried samples belong to Stage II or Stage III input lists
samDf$inGrp1 <- samDf$sample.submitter_id%in%group1$Sample.ID
samDf$inGrp2 <- samDf$sample.submitter_id%in%group2$Sample.ID

# Logical vector: TRUE if sample is in either group
samDf$sampleLogic <- samDf$inGrp1 | samDf$inGrp2

# Extract matching TCGA case barcodes
barcodeSelected <- samDf$cases[samDf$sampleLogic]

# ----- Refined GDC Query -----------------------------------------------------
query_LIHC <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access = "open",
  barcode = barcodeSelected
)

# Download and prepare data into SummarizedExperiment object
dds1 <- GDCprepare(
  query = query_LIHC,
  directory = data_dir
)

# ----- Clean Tumor Stage Labels ----------------------------------------------
# Inspect original stage labels
table(dds1$ajcc_pathologic_stage)

# Simplify stage labels (e.g., "Stage IIIA", "Stage IIIB" → "III")
dds1$ajcc_pathologic_stage <- gsub('.*III.*', 'III', dds1$ajcc_pathologic_stage)
dds1$ajcc_pathologic_stage <- gsub('.* II.*', 'II', dds1$ajcc_pathologic_stage)

# Verify cleaned labels
table(dds1$ajcc_pathologic_stage)

# ----- Filter Samples by Stage ------------------------------------------------
# Replace NA stage values with "NA" string for filtering
dds1$ajcc_pathologic_stage[is.na(dds1$ajcc_pathologic_stage)] <- "NA"

# Keep only Stage II and Stage III samples
keepLogic <- dds1$ajcc_pathologic_stage == "III" | dds1$ajcc_pathologic_stage == "II"
dds1 <- dds1[, keepLogic]

# Confirm counts
table(dds1$ajcc_pathologic_stage)

# ----- Extract Sample IDs for Each Group -------------------------------------
group1_keep <- dds1$sample[(dds1$ajcc_pathologic_stage == "II")]
group2_keep <- dds1$sample[(dds1$ajcc_pathologic_stage == "III")]

# Subset original group files to retain only valid samples
newgroup1 <- group1[group1$Sample.ID %in% group1_keep, ]
newgroup2 <- group2[group2$Sample.ID %in% group2_keep, ]

# ----- Manual Adjustment -----------------------------------------------------
# Remove one sample to balance group sizes (e.g., ensure n = 50)
newgroup1 <- newgroup1[!newgroup1$Sample.ID %in% "TCGA-ZS-A9CD-01A",]

# ----- Save Cleaned Groups ---------------------------------------------------
write.csv(newgroup1, "Data/CleanGroup1.csv")
write.csv(newgroup2, "Data/CleanGroup2.csv")
