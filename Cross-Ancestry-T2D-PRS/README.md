# **Isolating the Drivers of Poor European-to-South Asian Type 2 Diabetes Risk Score Transferability**

---

## Overview

Polygenic risk scores (PRS) for type 2 diabetes (T2D) are trained predominantly on European (EUR)-ancestry GWAS, despite Central and South Asian (CSA) populations bearing the highest global disease burden. This project evaluates **why** EUR-derived T2D PRS transfer poorly to South Asian (SAS) individuals, testing three candidate explanations:

1. **Divergent genetic architecture** between populations
2. **Technical artifacts** (mismatched LD structure and allele frequencies)
3. **Limited discovery power** in non-European GWAS

**Hypothesis:** Limited CSA GWAS discovery power — not divergent architecture or technical artifacts — is the dominant driver of poor cross-ancestry T2D PRS transferability.

---

## Reports

- **Oral Presentation:** [File](Cross-Ancestry-T2D-PRS/oral_presentation/T2D_PRS_Presentation.pdf); [Presentation](https://youtu.be/e-Qbfxc2Rp8)
- **Poster:** [File](Cross-Ancestry-T2D-PRS/poster/SBMI-Poster_T2D-PRS_Kim-Luong.pdf); [Presentation](https://youtu.be/ym2-KCxHgV0)
- **Manuscript:** [File](Cross-Ancestry-T2D-PRS/paper/BMI6313_FullPaper_KimLuong.pdf)

---

## Key Findings

- EUR-derived PRS shifted the SAS score distribution upward by **1.27–1.32 SD**, producing ~4-fold enrichment of SAS individuals in the EUR-defined top risk decile
- Cross-score concordance between EUR- and CSA-derived scores was **negligible** (r = -0.001 to 0.047)
- Genome-wide genetic correlation was **indistinguishable from unity** (rg = 1.27, p = 6.7×10⁻¹¹), and 0 of 8,789 tested loci showed significant EUR–CSA heterogeneity after Bonferroni correction
- Directional concordance of effect estimates rose from **63.7%** at low statistical power to **100%** at high power
- **Conclusion:** Miscalibration stems mainly from limited CSA discovery sample size, not fundamental differences in diabetes biology. Expanding CSA-specific GWAS — not building a more sophisticated scoring model — is the most direct path to equitable T2D risk prediction.

---

## Data Sources

| Source | Description |
|---|---|
| [Pan-UK Biobank](https://pan.ukbb.broadinstitute.org/) | EUR and CSA T2D (ICD-10: E11) GWAS summary statistics |
| [1000 Genomes Project Phase 3](https://www.internationalgenome.org/) | Individual-level EUR and SAS genotypes (ancestry-matched reference panels) |

All data used are publicly available and de-identified.

---

## Methods Summary

- **Quality control:** MAF > 1%, missingness < 5%, Hardy-Weinberg p > 1×10⁻⁶; ancestry verified via PCA (503 EUR, 478 SAS individuals retained)
- **PRS methods:** Clumping-and-thresholding (C+T; LD r² < 0.1) and LDpred2-inf (joint genome-wide LD modeling)
- **Score configurations:** 4 configurations isolating effect-size, LD-reference, and discovery-power contributions
- **Outcome measures:** Mean score shift, top-decile enrichment, cross-score concordance (Pearson r), locus-level heterogeneity (Cochran's Q, Bonferroni-corrected), genome-wide genetic correlation (rg, via LD score regression)

---

## Repository Structure

```
Cross-Ancestry-T2D-PRS/
├── README.md                # This file
├── scripts/                 # Analysis code (QC, PCA, PRS construction, scoring, statistics)
├── results/                 # Output tables, figures, and summary statistics
├── poster/                  # Conference/course poster
├── paper/                   # Full written manuscript
└── oral presentation/       # Slide deck for oral presentation
```

---

## Acknowledgements

Thanks to classmates in BMI 6313: Scientific Writing in Healthcare, McWilliams School of Biomedical Informatics, UTHealth Houston, for feedback and peer review. Thanks also to the Pan-UK Biobank and 1000 Genomes Project consortia for making their data publicly available.

---

## Course Context
**BMI 6313:** Scientific Writing in Healthcare

**Institution:** The University of Texas Health Science Center at Houston 

**Presented:** Oral Presentation and Research Poster, Summer 2026
