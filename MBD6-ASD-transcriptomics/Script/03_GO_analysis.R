# =============================================================================
# Script 03: Gene Ontology (GO) Enrichment Analysis
# Project:   MBD6 Haploinsufficiency and ASD - Transcriptomic Analysis
# Course:    BMI 5332: Statistical Analysis of Genomic Data
#            UT Health Houston, Spring 2026
# Input:     data/processed/DEG_results.rds  (from 02_differential_expression.R)
# Output:    results/tables/Table5_go_BP.csv
#            results/figures/Fig_GO_BP_dotplot.png
#            results/figures/Fig_GO_BP_cnetplot.png
#            results/figures/Fig_GO_BP_emapplot.png
# =============================================================================

# ----- Load Libraries --------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(HGNChelper)
library(ggplot2)
library(ggnewscale)

# ----- Set Working Directory -------------------------------------------------
setwd("~/Desktop/BMI 5332_R/Project")

# ----- Load DEG Results ------------------------------------------------------
deg_data <- readRDS("Data/Processed/DEG_results.rds")
results_unique <- deg_data$results_unique
sig_DEGs <- deg_data$sig_DEGs
upregulated <- deg_data$upregulated
downregulated <- deg_data$downregulated


# =============================================================================
# STEP 1: Symbol Rescue with HGNChelper + Entrez ID Mapping
# =============================================================================
rescue_symbols <- function(gene_vec) {
  checked <- checkGeneSymbols(gene_vec)
  updated <- ifelse(!is.na(checked$Suggested.Symbol) &
                      checked$Suggested.Symbol != "",
                    checked$Suggested.Symbol,
                    gene_vec)
  bitr(updated, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}

sig_entrez_updated        <- rescue_symbols(sig_DEGs$GeneSymbol)
up_entrez_updated         <- rescue_symbols(upregulated$GeneSymbol)
down_entrez_updated       <- rescue_symbols(downregulated$GeneSymbol)
background_entrez_updated <- rescue_symbols(results_unique$GeneSymbol)

cat("Background mapped: ", nrow(background_entrez_updated), "\n")
cat("Sig DEGs mapped: ", nrow(sig_entrez_updated), "\n")
cat("Upregulated mapped: ", nrow(up_entrez_updated), "\n")
cat("Downregulated mapped: ", nrow(down_entrez_updated), "\n")
cat("Sig DEG failure rate: ", round((nrow(sig_DEGs) - nrow(sig_entrez_updated)) / nrow(sig_DEGs) * 100, 1), "%\n")


# =============================================================================
# STEP 2: GO Enrichment
# =============================================================================
run_enrichGO <- function(entrez_ids, background_ids, ontology) {
  enrichGO(
    gene = entrez_ids,
    universe = background_ids,
    OrgDb = org.Hs.eg.db,
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
}

go_BP <- run_enrichGO(sig_entrez_updated$ENTREZID, background_entrez_updated$ENTREZID, "BP")
go_MF <- run_enrichGO(sig_entrez_updated$ENTREZID, background_entrez_updated$ENTREZID, "MF")
go_CC <- run_enrichGO(sig_entrez_updated$ENTREZID, background_entrez_updated$ENTREZID, "CC")

cat("GO BP enriched terms: ", nrow(as.data.frame(go_BP)), "\n")
cat("GO MF enriched terms: ", nrow(as.data.frame(go_MF)), "\n")
cat("GO CC enriched terms: ", nrow(as.data.frame(go_CC)), "\n")

go_BP_up <- run_enrichGO(up_entrez_updated$ENTREZID, background_entrez_updated$ENTREZID, "BP")
go_BP_down <- run_enrichGO(down_entrez_updated$ENTREZID, background_entrez_updated$ENTREZID, "BP")

cat("GO BP upregulated enriched terms: ", nrow(as.data.frame(go_BP_up)), "\n")
cat("GO BP downregulated enriched terms: ", nrow(as.data.frame(go_BP_down)), "\n")

# Save GO BP table
go_BP_df <- as.data.frame(go_BP)[, c("Description", "GeneRatio", "BgRatio", "p.adjust", "geneID")]
write.csv(go_BP_df, "Results/Table5_go_BP.csv", row.names = FALSE)


# =============================================================================
# STEP 3: GO Visualizations
# =============================================================================

# --- Dotplot -----------------------------------------------------------------
png("Results/Fig7_GO_BP_dotplot.png", width = 10, height = 8, units = "in", res = 300)

dotplot(go_BP, 
        showCategory = 20,
        title = "GO Biological Process Enrichment", 
        font.size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

dev.off()

# --- cnetplot (genes driving pathways) ---------------------------------------
fc_vector <- ifelse(sig_DEGs$logFC > 0, 1, -1)
names(fc_vector) <- sig_DEGs$GeneSymbol

png("Results/Fig8_GO_BP_cnetplot.png", width = 15, height = 10, units = "in", res = 300)

p <- cnetplot(go_BP, 
              showCategory = 10, 
              foldChange = fc_vector) +
  scale_color_gradient2(low = "blue", 
                        mid = "white", 
                        high = "red",
                        midpoint = 0, 
                        guide = "none")
p + new_scale_color() +
  geom_point(
    data = data.frame(
      x = c(Inf, Inf), 
      y = c(Inf, Inf),
      group = c("Downregulated", "Upregulated")),
    aes(x = x, y = y, color = group), size = 5) +
  scale_color_manual(
    name = "Sig DEGs",
    values = c("Downregulated" = "blue", "Upregulated" = "red")) +
  labs(title = "GO BP Gene-Concept Network") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.margin = margin(10, 25, 10, 10))

dev.off()

# --- emapplot (relationships between enriched terms) -------------------------
go_BP_sim <- pairwise_termsim(go_BP)

png("Results/Unused_GO_BP_emapplot.png", width = 12, height = 10, units = "in", res = 300)

emapplot(go_BP_sim, showCategory = 28) +
  labs(title = "GO BP Enrichment Map") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

dev.off()


# =============================================================================
# SAVE: GO results for reference
# =============================================================================
saveRDS(list(
  go_BP = go_BP,
  go_MF = go_MF,
  go_CC = go_CC,
  go_BP_up = go_BP_up,
  go_BP_down = go_BP_down
), file = "Data/Processed/GO_results.rds")