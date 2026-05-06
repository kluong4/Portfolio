# MBD6 Haploinsufficiency and Autism Spectrum Disorder: Transcriptomic Analysis

## Overview
This project investigates the genome-wide transcriptional consequences of MBD6 haploinsufficiency in a neuronal model of Autism Spectrum Disorder (ASD). Using publicly available Agilent microarray data from MBD6-heterozygous SH-SY5Y human neuroblastoma cells, we identified 913 differentially expressed genes and linked them to ASD-associated genetic variants through eQTL-GWAS integration.

---

## Background
MBD6 (Methyl-CpG Binding Domain Protein 6) is a chromatin regulator implicated in ASD through rare copy number variants and de novo mutations. This analysis uses an isogenic neuronal cell model to characterize the downstream transcriptional consequences of MBD6 loss, with the goal of identifying disrupted biological pathways and ASD-relevant genetic signals.

**GEO Accession:** [GSE314093](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314093)  
**Platform:** GPL13497 — Agilent-028004 SurePrint G3 Human GE 8x60K Microarray

---

## Methods

| Step | Description | Tools |
|------|-------------|-------|
| Preprocessing | Background correction (normexp), quantile normalization, probe filtering | `limma` |
| Quality Control | PCA, sample correlation heatmap, MA plots | `ggplot2`, `reshape2` |
| Differential Expression | Linear model with empirical Bayes moderation, BH correction | `limma` |
| GO Enrichment | Biological Process, Molecular Function, Cellular Component | `clusterProfiler`, `org.Hs.eg.db` |
| eQTL-GWAS Integration | GTEx v11 Brain Cortex eQTLs overlapped with ASD GWAS variants (±500kb window) | `arrow`, `dplyr` |

---

## Key Results
- **913 DEGs** identified between MBD6-heterozygous and isogenic control cells (637 upregulated, 276 downregulated; FDR < 0.05, |log2FC| > 1)
    ![Fig5_VolcanoPlot](/MBD6-ASD-transcriptomics/Results/Figures/Fig5_VolcanoPlot.png)

    ![Fig6_DEG_Heatmap](/MBD6-ASD-transcriptomics/Results/Figures/Fig6_DEG_Heatmap.png)
  
- **Top dysregulated pathways:** immune activation, ERK/MAPK signaling, and stress-activated kinase signaling (GO Biological Process enrichment)
    ![Fig8_cnetplot](/MBD6-ASD-transcriptomics/Results/Figures/Fig8_GO_BP_cnetplot.png)
  
- **22 DEGs** harbored eQTLs overlapping ASD GWAS variants within ±500kb, including mitochondrial-linked genes **NDUFB4**, **GFM1**, and **SRGAP3**, linking MBD6 haploinsufficiency to known ASD genetic architecture
    ![Fig9_eQTL_GWAS_overlap](/MBD6-ASD-transcriptomics/Results/Figures/Fig9_eQTL_GWAS_overlap.png)
  
---

## Repository Structure

```
MBD6-ASD-transcriptomics/
│
├── README.md
│
├── Scripts/
│   ├── 01_EDA.R                  # Data loading, normalization, QC plots
│   ├── 02_DE_analysis.R          # DEG analysis, volcano plot, heatmap
│   ├── 03_GO_analysis.R          # GO enrichment + visualizations
│   └── 04_eQTL_GWAS_analysis.R   # GTEx eQTL x ASD GWAS overlap
│
├── Results/
│   ├── Figures/                   # All output plots (.png)
│   └── Tables/                    # All output tables (.csv)
│
├── MBD6_ASD_Report.pdf            # Full Report (.pdf)
│
└── 
```

> **Note:** Raw data files are not included due to size. See the Data Access section below.

---

## Requirements

**R version:** 4.3+

```r
install.packages(c("ggplot2", "ggrepel", "reshape2", "dplyr", "arrow"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "limma", "GEOquery", "ComplexHeatmap",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot"
))

install.packages("HGNChelper")
install.packages("ggnewscale")
```

**Run scripts in order:**
```
01_preprocessing.R → 02_differential_expression.R → 03_GO_enrichment.R → 04_eQTL_GWAS_integration.R
```

---

## Data Access

| Dataset | Source | Link |
|---------|--------|------|
| Microarray expression data | GEO: GSE314093 | [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314093) |
| Platform annotation | GEO: GPL13497 | [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13497) |
| GTEx v11 Brain Cortex eQTLs | GTEx Portal | [Link](https://gtexportal.org/home/downloads/adult-gtex/qtl) |
| ASD GWAS associations | GWAS Catalog (MONDO_0005258) | [Link](https://www.ebi.ac.uk/gwas/efotraits/MONDO_0005258) |

---

## Course Context
**BMI 5332:** Statistical Analysis of Genomic Data  
The University of Texas Health Science Center at Houston  
Spring 2026
