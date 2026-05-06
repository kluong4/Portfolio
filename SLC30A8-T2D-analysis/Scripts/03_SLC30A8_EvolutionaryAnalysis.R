# =============================================================================
# Script 03: Evolutionary Conservation Analysis of SLC30A8
# Project:   Population and Evolutionary Analysis of SLC30A8 in Type 2 Diabetes
# Course:    BMI 5330: Introduction to Bioinformatics
#            UT Health Houston, Fall 2025
#
# Input:     Data/ProcessedData_RDS/variants.rds        (from 02_population_variation.R)
#            Data/SLC30A8_Data/100vertsl.wig                 (PhyloP 100-way vertebrate scores)
#            Data/SLC30A8_Data/SLC30A8_exons.bed             (from 01_extract_SLC30A8_variants.sh)
# Output:    Results/Fig2_conservation_plot.png
#            Results/SLC30A8_variants.tsv
#
# PhyloP data source:
#   UCSC Genome Browser — PhyloP 100-way vertebrate conservation (GRCh38/hg38)
#   https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/
# =============================================================================

# ----- Load Libraries --------------------------------------------------------
library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)
library(readr)

# ----- Set Working Directory -------------------------------------------------
setwd("~/Desktop/Prev Courses/25Fall/BMI 5330/Project/SLC30A8-T2D-analysis")


# ----- Configuration ---------------------------------------------------------
PHYLOP_THRESHOLD <- 2.27    # PhyloP score cutoff for high conservation

# ----- Load Data from Script 02 ----------------------------------------------
variants <- readRDS("Data/ProcessedData_RDS/variants.rds")


# =============================================================================
# STEP 1: Load PhyloP Conservation Scores
# =============================================================================
phyloP <- fread(
  "Data/SLC30A8_Data/100vertsl.wig",
  skip = "variableStep",
  col.names = c("position", "phyloP")
)

cat("PhyloP positions loaded: ", nrow(phyloP), "\n")
cat("Position range: ", min(phyloP$position), "–", max(phyloP$position), "\n")


# =============================================================================
# STEP 2: Join PhyloP Scores to Variants
# =============================================================================
variants_phyloP <- variants %>%
  left_join(phyloP %>% select(position, phyloP),
            by = c("start" = "position"))
variants$phyloP <- variants_phyloP$phyloP


# =============================================================================
# STEP 3: Flag Highly Conserved Variants
# =============================================================================
variants <- variants %>% mutate(highly_conserved = phyloP > PHYLOP_THRESHOLD)

# Conservation summary (threshold =  2.27 )
table(variants$highly_conserved, useNA = "ifany")

# Subset highly conserved variants
high_conserved_variants <- variants %>%
  filter(phyloP > PHYLOP_THRESHOLD) %>%
  arrange(start)

cat("Highly conserved variants: ", nrow(high_conserved_variants), "\n")


# =============================================================================
# STEP 4: Load Exon Coordinates for Gene Structure Annotation
# =============================================================================
exon_coord <- read_delim("Data/SLC30A8_Data/SLC30A8_exons.bed",
                         delim = "\t", escape_double = FALSE,
                         col_names = FALSE, trim_ws = TRUE)
colnames(exon_coord) <- c("chr","start","end","info","dot","strand")


# =============================================================================
# STEP 5: Conservation Plot — PhyloP Scores + Variants + Exon Structure
# =============================================================================
phyloP_max <- max(phyloP$phyloP, na.rm = TRUE)

ggplot() +
  
  # Conservation score line
  geom_line(data = phyloP,
            aes(x = position, y = phyloP),
            color = "lightblue") +
  
  # All variants — shape and color by conservation status
  geom_point(data = variants %>% filter(!is.na(phyloP)),
             aes(x     = start,
                 y     = phyloP,
                 shape = highly_conserved,
                 color = highly_conserved),
             size = ifelse(variants$highly_conserved[!is.na(variants$phyloP)], 5, 3)) +
  
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 8),
    labels = c("FALSE" = "Not Highly Conserved", "TRUE" = "Highly Conserved"),
    name   = "Conservation Status") +
  
  scale_color_manual(
    values = c("FALSE" = "grey40", "TRUE" = "red"),
    labels = c("FALSE" = "Not Highly Conserved", "TRUE" = "Highly Conserved"),
    name   = "Conservation Status") +
  
  # Exon structure track (ENST00000456015.7)
  geom_rect(data = exon_coord,
            aes(xmin = start, xmax = end,
                ymin = phyloP_max + 1, ymax = phyloP_max + 2),
            fill = "black") +
  
  geom_hline(yintercept = phyloP_max + 1.5,
             color = "black", linetype = "solid") +
  
  geom_text(aes(x     = max(phyloP$position) - 3000,
                y     = phyloP_max + 0.5,
                label = "ENST00000456015.7"),
            size = 5) +
  
  # Zero reference line
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  
  labs(
    title = "SLC30A8 Conservation and Variants",
    x     = "Genomic Position (chr8)",
    y     = "PhyloP Conservation Score"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Results/Fig2_ConservationPlot.png", width = 15, height = 5, dpi = 300)


# =============================================================================
# STEP 6: Save Annotated Variant Tables
# =============================================================================
# Full annotated variants table (with phyloP, VT, AA — no highly_conserved flag)
variants_final <- variants %>% select(-highly_conserved)

write.table(variants_final, "Results/SLC30A8_variants.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Highly conserved variants only
high_conserved_final <- high_conserved_variants %>% select(-highly_conserved)
#write.table(high_conserved_final, "Results/SLC30A8_high_conserved_variants.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

rsID_HighlyConserved <- paste(high_conserved_final$id, collapse = ", ")
#write(rs_string, file = "variants_rsIDs.txt")
