# Kim Luong — Bioinformatics Portfolio

M.S. Biomedical Informatics student at UT Health Houston with a B.S. in Biology from UT Austin. My work sits at the intersection of computational biology, genomics, and health data science — applying statistical modeling, bioinformatics pipelines, and data visualization to answer biological and clinical questions.

📧 kim2002.kl@gmail.com &nbsp;|&nbsp; 📍 San Antonio, TX &nbsp;|&nbsp; 🔗 [LinkedIn](https://www.linkedin.com/in/kim-ngoc-luong)

---

## Technical Skills

| Category | Tools & Languages |
|----------|------------------|
| **Languages** | Python, R, Bash/Shell, SQL |
| **Bioinformatics** | DESeq2, TCGAbiolinks, Bioconductor, limma, clusterProfiler, bedtools |
| **Data Analysis** | pandas, NumPy, dplyr, tidyr, data.table |
| **Visualization** | ggplot2, seaborn, matplotlib, ggrepel, ComplexHeatmap |
| **Databases** | TCGA/GDC, 1000 Genomes Project, UCSC Genome Browser, GEUVADIS, dbSNP, GWAS Catalog, NCBI, GEO |

---

## Projects

### 1. [MBD6 Haploinsufficiency and Autism Spectrum Disorder: Transcriptomic Analysis](https://github.com/kluong4/Portfolio/tree/main/MBD6-ASD-transcriptomics)
**Microarray Differential Expression Analysis | R, Bioconductor**

Genome-wide transcriptional analysis of MBD6 haploinsufficiency in a neuronal model of Autism Spectrum Disorder using Agilent microarray data from GEO (GSE314093).

- **Analysis:** Background correction & quantile normalization → limma linear modeling → GO enrichment → GTEx eQTL × ASD GWAS variant integration (±500kb window)
- **Key result:** 913 DEGs identified (637 up, 276 down); immune activation, ERK/MAPK, and stress-activated kinase pathways enriched; 22 DEGs overlap ASD GWAS variants including mitochondrial-linked genes NDUFB4, GFM1, and SRGAP3
- **Tools:** R, `limma`, `clusterProfiler`, `ComplexHeatmap`, `ggplot2`, `ggrepel`, `arrow`, `dplyr`

---

### 2. [Differential Gene Expression Analysis of Liver Hepatocellular Carcinoma (LIHC)](https://github.com/kluong4/Portfolio/tree/main/TCGA-LIHC-DEGanalysis)
**RNA-seq Differential Expression | R, DESeq2, Bioconductor**

End-to-end differential expression pipeline comparing AJCC Stage II vs. Stage III hepatocellular carcinoma tumors using TCGA-LIHC RNA-seq data accessed via the GDC portal.

- **Analysis:** GDC data query → PCA-based QC & outlier removal → DESeq2 negative binomial GLM → BH multiple testing correction → volcano plot visualization
- **Key result:** 715 significant DEGs (FDR < 0.05); RAB25 (log2FC = 5.45) and GCK (log2FC = −4.47) identified as top candidates linking tumor progression to metabolic reprogramming
- **Tools:** R, `TCGAbiolinks`, `DESeq2`, `ggplot2`, `ggrepel`

---

### 3. [Population and Evolutionary Analysis of SLC30A8 in Type 2 Diabetes](https://github.com/kluong4/Portfolio/tree/main/SLC30A8-T2D-analysis)
**Population Genetics · Evolutionary Genomics · Transcriptomics | Bash/Shell, R**

Integrative three-part analysis of the T2D-associated zinc transporter gene SLC30A8 using public genomic databases, combining population genetics, evolutionary conservation, and RNA-seq expression analysis.

- **Analysis 1 (Bash/Shell + R):** VCF extraction from 1000 Genomes Phase 1, bedtools exon/intron classification, population-stratified allele frequency analysis (AFR, AMR, ASN, EUR)
- **Analysis 2 (R):** PhyloP 100-way conservation scoring, high-conservation variant flagging (PhyloP > 2.27), exon structure overlay visualization
- **Analysis 3 (R):** GEUVADIS RNA-seq FPKM processing, per-population mean expression, GWAS gene filtering, variance ranking across 5 populations
- **Key result:** SLC30A8 variants are globally rare and non-population-specific; ranked 101st in expression variance among T2D GWAS genes — consistent with strong purifying selection
- **Tools:** Bash/Shell (`wget`, `awk`, `gzip`, `sed`, `bedtools`), R (`stringr`, `dplyr`, `ggplot2`, `data.table`, `readxl`)

---

### 4. [Chronic Disease Prevalence & Medicare Utilization Analysis](https://github.com/kluong4/Portfolio/tree/main/chronic-disease-analysis)
**Exploratory Data Analysis & Visualization | Python**

Exploratory analysis of chronic disease prevalence and Medicare utilization patterns across 21 conditions using CMS public data, with a focus on demographic disparities and geographic variation.

- **Analysis:** Prevalence heatmaps stratified by race and age group → Medicare utilization benchmarking (spending, readmissions, ER visits) across TX, CA, NY, FL → state-wise diverging bar chart analysis by gender and race
- **Key result:** Non-Hispanic Black beneficiaries show disproportionate cardiovascular disease burden; Native American populations exhibit the greatest geographic variability across conditions; Drug/Substance Abuse drives ER utilization far above other conditions
- **Tools:** Python, `pandas`, `NumPy`, `matplotlib`, `seaborn`

---

### 5. [Possum Morphology Regression Analysis](https://github.com/kluong4/Portfolio/tree/main/PossumMorphology-RegressionAnalysis)
**Multiple Linear Regression · Model Diagnostics | R**

Regression analysis of morphometric data from 104 mountain brushtail possums across 7 Australian trapping sites to identify predictors of total body length.

- **Analysis:** Multiple linear regression → ANOVA model comparison → backward selection → interaction testing → VIF multicollinearity check → Cook's D influence analysis → LOO PRESS cross-validation
- **Key result:** Model explains 81.2% of variance in total body length (R² = 0.8118); geographic site is a significant predictor (p = 2.5×10⁻¹²); skull width is non-significant after controlling for other variables
- **Tools:** R, `DAAG`, `GGally`, `ggfortify`, `car`

---

### 6. [Time Budget of Feeding in Northern Shovelers](https://github.com/kluong4/Portfolio/tree/main/NorthernShovelers-BehaviorAnalysis)
**Behavioral Ecology · Nonparametric Statistics | R**

Field observation study and statistical analysis of foraging behavior in Northern Shovelers (*Anas clypeata*) at Mueller Lake Park, Austin, TX across 12 days of scan sampling.

- **Analysis:** Shapiro-Wilk normality testing → Wilcoxon-Mann Whitney test (sex effect) → Kruskal-Wallis test (time-of-day effect) → paired t-tests across 4 foraging behavior categories
- **Key result:** Foraging frequency increases significantly across the day (p < 2.2×10⁻¹⁶); sex shows no significant effect (p = 0.1), suggesting human park activity as a behavioral confound
- **Tools:** R (base R, `ggplot2`)

---

### 7. [Train & Engineer Workers - Regression Analysis](https://github.com/kluong4/Portfolio/tree/main/T%26Eworkers-RegressionAnalysis)
**Linear Regression · Interaction Modeling | R**

Statistical analysis of federal railroad workforce data (246 records, FRA 2018) examining how work schedule type and sick leave usage predict total years of service.

- **Analysis:** Simple linear regression → interaction modeling → interaction plot visualization
- **Key result:** Effect of sick days on retention differs significantly by schedule type, providing quantitative evidence for schedule-based workforce policy
- **Tools:** R (base R, `ggplot2`)

---

### 8. [Lipidomics of Copper Stress on *Arabidopsis thaliana*](https://github.com/kluong4/Portfolio/tree/main/Arabidopsis-Lipidomics)
**Wet Lab · Lipid Extraction · Hydroponics | Freshman Research Initiative, UT Austin**

Controlled hydroponic experiment investigating the morphological effects of copper toxicity on *Arabidopsis thaliana*, with a lipid extraction pipeline designed for downstream LC-MS profiling.

- **Analysis:** Hydroponic culture (control vs. 10x copper) → root length measurement → pigmentation and biomass assessment → MTBE/MeOH lipid extraction (LC-MS results not collected within project timeline)
- **Key result:** 57% reduction in average root length in copper-stressed plants (3 cm vs. 7 cm); reduced pigmentation and biomass consistent with photosynthetic disruption
- **Tools:** Wet lab (hydroponic culture, liquid nitrogen freeze-drying, MTBE/MeOH extraction)

---

## Education

**M.S. Biomedical Informatics** — UT Health Houston (2025–2027, Expected)

**B.S. Biology** — The University of Texas at Austin (2021–2025)
- Pre-Health Professions for Science Majors Certificate
- University Honors: Fall 2021, 2022, 2023
