# =============================================================================
# Script 09: Functional Annotation of Heterogeneous Loci
# Project:   Cross-ancestry T2D PRS Transferability (EUR -> CSA)
# =============================================================================
# OVERVIEW:
# This script performs functional annotation of EUR GWS T2D loci,
# with emphasis on loci showing cross-ancestry effect size heterogeneity.
# Annotation provides biological context for WHY certain loci show
# differential effects between EUR and CSA populations.
#
# ANALYSES:
#   1. Genic feature annotation (promoter, exon, intron, intergenic)
#      using annotatr
#   2. CpG island annotation
#   3. Gene ontology (GO) and KEGG pathway enrichment of nearest genes
#      using clusterProfiler
#   4. Comparison of functional categories between heterogeneous vs
#      non-heterogeneous loci
#   5. Summary visualization
#
# PACKAGES:
#   annotatr   -- regulatory/genic feature annotation
#   clusterProfiler -- pathway enrichment
#   TxDb.Hsapiens.UCSC.hg19.knownGene -- gene models (GRCh37/hg19)
#   org.Hs.eg.db -- human gene annotation database
#
# NOTE: Pan-UKB uses GRCh37 coordinates so we use hg19 annotation.
#
# OUTPUTS:
#   results/annotation/gws_annotated.csv
#   results/annotation/go_enrichment_all.csv
#   results/annotation/go_enrichment_heterogeneous.csv
#   results/annotation/kegg_enrichment.csv
#   results/figures/fig9_functional_annotation.png
#   results/figures/fig10_pathway_enrichment.png
# =============================================================================

library(data.table)
library(ggplot2)
library(patchwork)

dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures",    recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# STEP 1: Install and load annotation packages
# -----------------------------------------------------------------------------

cat("Loading annotation packages...\n\n")

# Install if needed
pkg_needed <- c("annotatr", "clusterProfiler", "org.Hs.eg.db",
                "TxDb.Hsapiens.UCSC.hg19.knownGene",
                "GenomicRanges", "IRanges")

for(pkg in pkg_needed) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    if(pkg %in% c("annotatr", "clusterProfiler", "org.Hs.eg.db",
                  "TxDb.Hsapiens.UCSC.hg19.knownGene",
                  "GenomicRanges", "IRanges")) {
      BiocManager::install(pkg, ask = FALSE)
    }
  }
}

# Install BiocManager if needed
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

library(annotatr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)

cat("All packages loaded.\n\n")

# -----------------------------------------------------------------------------
# STEP 2: Load heterogeneity results
# -----------------------------------------------------------------------------

cat("Loading heterogeneity results...\n\n")

het <- fread("results/concordance/heterogeneity_results.csv")

cat("Total GWS loci:", nrow(het), "\n")
cat("Columns:", colnames(het), "\n\n")

# Define heterogeneous loci -- nominal Q p < 0.05
# (Bonferroni threshold too conservative given low CSA power)
het[, is_heterogeneous := Q_pval < 0.05]
cat("Heterogeneous loci (Q p < 0.05):", sum(het$is_heterogeneous), "\n")
cat("Non-heterogeneous loci:         ", sum(!het$is_heterogeneous), "\n\n")

# -----------------------------------------------------------------------------
# STEP 3: Create GRanges object for annotation
# -----------------------------------------------------------------------------
# annotatr requires a GRanges object with chromosome, start, end.
# Pan-UKB uses GRCh37 (hg19) coordinates -- chromosome names need
# "chr" prefix for UCSC-style annotation.

cat("Creating GRanges object...\n")

# Filter to autosomes with valid positions
het_clean <- het[!is.na(chr) & !is.na(pos) & chr %in% 1:22]
cat("Loci on autosomes:", nrow(het_clean), "\n\n")

# Create GRanges -- SNPs are single-base positions (start = end = pos)
gws_gr <- GRanges(
  seqnames = paste0("chr", het_clean$chr),
  ranges   = IRanges(start = het_clean$pos,
                     end   = het_clean$pos),
  varID    = het_clean$varID,
  beta_EUR = het_clean$beta_EUR,
  beta_CSA = het_clean$beta_CSA,
  Q_pval   = het_clean$Q_pval,
  is_het   = het_clean$is_heterogeneous
)

cat("GRanges created:", length(gws_gr), "ranges\n\n")

# -----------------------------------------------------------------------------
# STEP 4: Annotate with annotatr
# -----------------------------------------------------------------------------
# We annotate against hg19 genic features:
#   hg19_genes_promoters    -- 1-2kb upstream of TSS
#   hg19_genes_5UTRs        -- 5' untranslated regions
#   hg19_genes_exons        -- coding exons
#   hg19_genes_introns      -- intronic regions
#   hg19_genes_3UTRs        -- 3' untranslated regions
#   hg19_genes_intergenic   -- between genes
#   hg19_cpg_islands        -- CpG islands
#   hg19_cpg_shores         -- CpG shores (2kb flanking islands)

cat("Building annotation database...\n")

annot_types <- c(
  "hg19_genes_promoters",
  "hg19_genes_5UTRs",
  "hg19_genes_exons",
  "hg19_genes_introns",
  "hg19_genes_3UTRs",
  "hg19_genes_intergenic",
  "hg19_cpg_islands",
  "hg19_cpg_shores"
)

annotations <- build_annotations(genome = "hg19", annotations = annot_types)
cat("Annotation database built.\n\n")

# Annotate GWS loci
cat("Annotating GWS loci...\n")
gws_annotated <- annotate_regions(
  regions     = gws_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet         = FALSE
)

# Convert to data frame
gws_annot_df <- as.data.frame(gws_annotated)
cat("Annotation complete.\n")
cat("Annotated rows:", nrow(gws_annot_df), "\n\n")

# Simplify annotation type labels
gws_annot_df$annot_simple <- gsub("hg19_genes_", "", gws_annot_df$annot.type)
gws_annot_df$annot_simple <- gsub("hg19_cpg_",   "CpG_", gws_annot_df$annot_simple)

# Count annotation categories
annot_counts <- as.data.frame(table(gws_annot_df$annot_simple))
colnames(annot_counts) <- c("Feature", "Count")
annot_counts <- annot_counts[order(-annot_counts$Count), ]
cat("Annotation category counts:\n")
print(annot_counts)
cat("\n")

# -----------------------------------------------------------------------------
# STEP 5: Compare annotation between heterogeneous and non-heterogeneous loci
# -----------------------------------------------------------------------------

cat("Comparing annotation between heterogeneous vs non-heterogeneous loci...\n\n")

# Proportion of each feature type in heterogeneous vs non-heterogeneous
het_annot    <- gws_annot_df[gws_annot_df$is_het == TRUE,  ]
nonhet_annot <- gws_annot_df[gws_annot_df$is_het == FALSE, ]

het_prop    <- prop.table(table(het_annot$annot_simple))    * 100
nonhet_prop <- prop.table(table(nonhet_annot$annot_simple)) * 100

# Combine into comparison table
feature_comparison <- data.table(
  Feature     = names(het_prop),
  Het_pct     = as.numeric(het_prop),
  NonHet_pct  = as.numeric(nonhet_prop[names(het_prop)])
)
feature_comparison[, Enrichment := Het_pct / NonHet_pct]
feature_comparison <- feature_comparison[order(-Enrichment)]

cat("Feature enrichment in heterogeneous vs non-heterogeneous loci:\n")
print(feature_comparison)
cat("\n")

# Save annotated results
fwrite(as.data.table(gws_annot_df),
       "results/annotation/gws_annotated.csv")
fwrite(feature_comparison,
       "results/annotation/feature_comparison.csv")
cat("Annotation results saved.\n\n")

# -----------------------------------------------------------------------------
# STEP 6: Extract nearest genes for pathway enrichment
# -----------------------------------------------------------------------------
# For pathway enrichment we need gene IDs (Entrez) for genes nearest
# to our GWS loci. annotatr provides gene_id in the annotation output
# for genic features.

cat("Extracting nearest genes...\n")

# Get gene IDs from annotated regions
# Filter to genic annotations (not intergenic or CpG)
genic_annot <- gws_annot_df[
  grepl("promoters|exons|introns|5UTR|3UTR", gws_annot_df$annot_simple), ]

# Extract unique Entrez gene IDs
gene_ids_all <- unique(na.omit(genic_annot$annot.gene_id))
cat("Unique genes at all GWS loci:          ", length(gene_ids_all), "\n")

# Genes at heterogeneous loci
gene_ids_het <- unique(na.omit(
  genic_annot$annot.gene_id[genic_annot$is_het == TRUE]
))
cat("Unique genes at heterogeneous loci:    ", length(gene_ids_het), "\n\n")

# Background gene list -- all genes in genome
all_genes <- keys(org.Hs.eg.db, keytype = "ENTREZID")

# -----------------------------------------------------------------------------
# STEP 7: GO enrichment analysis
# -----------------------------------------------------------------------------

cat("Running GO enrichment analysis...\n\n")

# GO enrichment for all GWS loci
go_all <- tryCatch(
  enrichGO(
    gene          = as.character(gene_ids_all),
    universe      = all_genes,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",          # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    readable      = TRUE           # convert IDs to gene symbols
  ),
  error = function(e) {
    cat("GO all error:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if(!is.null(go_all)) {
  cat("GO terms enriched (all GWS loci):", nrow(go_all@result), "\n")
  cat("Top 10 GO terms:\n")
  print(go_all@result[1:min(10, nrow(go_all@result)),
                      c("Description", "GeneRatio", "p.adjust")])
  fwrite(as.data.table(go_all@result),
         "results/annotation/go_enrichment_all.csv")
}

cat("\n")

# GO enrichment for heterogeneous loci
go_het <- tryCatch(
  enrichGO(
    gene          = as.character(gene_ids_het),
    universe      = as.character(gene_ids_all),  # use GWS genes as background
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    readable      = TRUE
  ),
  error = function(e) {
    cat("GO het error:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if(!is.null(go_het)) {
  cat("GO terms enriched (heterogeneous loci):", nrow(go_het@result), "\n")
  if(nrow(go_het@result) > 0) {
    cat("Top 10 GO terms:\n")
    print(go_het@result[1:min(10, nrow(go_het@result)),
                        c("Description", "GeneRatio", "p.adjust")])
  }
  fwrite(as.data.table(go_het@result),
         "results/annotation/go_enrichment_heterogeneous.csv")
}

cat("\n")

# KEGG pathway enrichment
cat("Running KEGG pathway enrichment...\n")

kegg_all <- tryCatch(
  enrichKEGG(
    gene          = as.character(gene_ids_all),
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  ),
  error = function(e) {
    cat("KEGG error:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if(!is.null(kegg_all)) {
  cat("KEGG pathways enriched:", nrow(kegg_all@result), "\n")
  if(nrow(kegg_all@result) > 0) {
    cat("Top 10 KEGG pathways:\n")
    print(kegg_all@result[1:min(10, nrow(kegg_all@result)),
                          c("Description", "GeneRatio", "p.adjust")])
  }
  fwrite(as.data.table(kegg_all@result),
         "results/annotation/kegg_enrichment.csv")
}

# -----------------------------------------------------------------------------
# STEP 8: Figure 9 -- Functional annotation summary
# -----------------------------------------------------------------------------

cat("\nGenerating Figure 9: Functional annotation...\n")

# Panel A: Overall annotation distribution
annot_plot_df <- data.table(
  Feature = annot_counts$Feature,
  Count   = annot_counts$Count,
  Pct     = annot_counts$Count / sum(annot_counts$Count) * 100
)

# Clean feature labels for plotting
annot_plot_df[, Feature_label := fcase(
  Feature == "promoters",   "Promoter",
  Feature == "exons",       "Exon",
  Feature == "introns",     "Intron",
  Feature == "5UTRs",       "5' UTR",
  Feature == "3UTRs",       "3' UTR",
  Feature == "intergenic",  "Intergenic",
  Feature == "CpG_islands", "CpG Island",
  Feature == "CpG_shores",  "CpG Shore",
  default = Feature
)]

p9a <- ggplot(annot_plot_df,
              aes(x = reorder(Feature_label, Pct),
                  y = Pct, fill = Feature_label)) +
  geom_col(alpha = 0.85) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "A. Genomic Feature Distribution",
    subtitle = "All EUR GWS T2D loci",
    x = NULL,
    y = "% of Annotations"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position  = "none")

# Panel B: Heterogeneous vs non-heterogeneous feature comparison
feat_plot_df <- melt(
  feature_comparison,
  id.vars      = "Feature",
  measure.vars = c("Het_pct", "NonHet_pct"),
  variable.name = "Group",
  value.name    = "Pct"
)
feat_plot_df[, Group := ifelse(Group == "Het_pct",
                               "Heterogeneous\n(Q p<0.05)",
                               "Non-heterogeneous")]

# Clean feature labels
feat_plot_df[, Feature_label := fcase(
  Feature == "promoters",   "Promoter",
  Feature == "exons",       "Exon",
  Feature == "introns",     "Intron",
  Feature == "5UTRs",       "5' UTR",
  Feature == "3UTRs",       "3' UTR",
  Feature == "intergenic",  "Intergenic",
  Feature == "CpG_islands", "CpG Island",
  Feature == "CpG_shores",  "CpG Shore",
  default = Feature
)]

p9b <- ggplot(feat_plot_df,
              aes(x = reorder(Feature_label, Pct),
                  y = Pct, fill = Group)) +
  geom_col(position = "dodge", alpha = 0.85) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Heterogeneous\n(Q p<0.05)" = "#e07b39",
    "Non-heterogeneous"         = "#4a90d9"
  )) +
  labs(
    title    = "B. Feature Distribution by Heterogeneity",
    subtitle = "Heterogeneous vs non-heterogeneous loci",
    x        = NULL,
    y        = "% of Annotations",
    fill     = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom")

fig9 <- p9a | p9b
fig9 <- fig9 + patchwork::plot_annotation(
  title = "Functional Annotation of EUR GWS T2D Loci",
  subtitle = paste0(
    "n = ", nrow(het_clean), " loci | ",
    "Heterogeneous: n = ", sum(het_clean$is_heterogeneous),
    " (Q p < 0.05)"
  )
)

ggsave("results/figures/fig9_functional_annotation.png",
       plot = fig9, width = 12, height = 6, dpi = 300)
cat("Saved: results/figures/fig9_functional_annotation.png\n")

# -----------------------------------------------------------------------------
# STEP 9: Figure 10 -- Pathway enrichment dotplot
# -----------------------------------------------------------------------------

cat("Generating Figure 10: Pathway enrichment...\n")

if(!is.null(go_all) && nrow(go_all@result) > 0) {
  
  # Select top 20 GO terms by adjusted p-value
  top_go <- as.data.table(go_all@result)[
    order(p.adjust)][1:min(20, .N)]
  
  # Parse GeneRatio to numeric
  top_go[, GR_num := sapply(GeneRatio, function(x) {
    parts <- as.numeric(strsplit(x, "/")[[1]])
    parts[1] / parts[2]
  })]
  
  fig10 <- ggplot(top_go,
                  aes(x     = GR_num,
                      y     = reorder(Description, -p.adjust),
                      color = p.adjust,
                      size  = Count)) +
    geom_point(alpha = 0.8) +
    scale_color_viridis_c(
      name      = "Adj. p-value",
      direction = -1,
      option    = "plasma"
    ) +
    scale_size_continuous(name = "Gene count",
                          range = c(2, 8)) +
    labs(
      title    = "GO Biological Process Enrichment at EUR GWS T2D Loci",
      subtitle = paste0("Top ", min(20, nrow(top_go)),
                        " terms | BH-adjusted p < 0.05"),
      x        = "Gene Ratio",
      y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(size = 9)
    )
  
  ggsave("results/figures/fig10_pathway_enrichment.png",
         plot = fig10, width = 10, height = 8, dpi = 300)
  cat("Saved: results/figures/fig10_pathway_enrichment.png\n\n")
  
} else {
  cat("No significant GO terms -- skipping Figure 10\n\n")
}

# -----------------------------------------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------------------------------------
cat("============================================================\n")
cat("SCRIPT 09 COMPLETE\n")
cat("============================================================\n")
cat("GWS loci annotated:              ", nrow(het_clean),                "\n")
cat("Heterogeneous loci (Q p<0.05):   ",
    sum(het_clean$is_heterogeneous),                                      "\n")
cat("Unique genes (all loci):         ", length(gene_ids_all),           "\n")
cat("Unique genes (het loci):         ", length(gene_ids_het),           "\n")
if(!is.null(go_all)) {
  cat("GO BP terms enriched (all):      ", nrow(go_all@result),          "\n")
}
if(!is.null(go_het)) {
  cat("GO BP terms enriched (het):      ", nrow(go_het@result),          "\n")
}
if(!is.null(kegg_all)) {
  cat("KEGG pathways enriched:          ", nrow(kegg_all@result),        "\n")
}
cat("Outputs saved to results/annotation/\n")
cat("Next step: 10_final_figures.R\n")
cat("============================================================\n")