
# Chronic Disease Prevalence & Medicare Utilization Analysis

Exploratory data analysis and visualization of chronic disease burden and healthcare utilization patterns among Medicare fee-for-service beneficiaries, using publicly available CMS data and Python.

---

## Overview

This project examines patterns in chronic disease prevalence and Medicare utilization across 21 conditions, with a focus on demographic disparities and geographic variation. The analysis is organized into three objectives, each producing interpretable visualizations that highlight trends across race/ethnicity, age group, gender, and U.S. geography.

---

## Skills Demonstrated

- Exploratory data analysis (EDA) on large healthcare datasets
- Multi-panel data visualization (heatmaps, grouped bar charts, diverging bar charts)
- Demographic and geographic stratification analysis
- Interpretation of Medicare utilization metrics (spending, readmissions, ER visits)
- Python: `pandas`, `numpy`, `matplotlib`, `seaborn`

---

## Dataset

| Field | Detail |
|-------|--------|
| **Source** | Centers for Medicare & Medicaid Services (CMS) Chronic Conditions Data Warehouse |
| **Link** | [data.cms.gov](https://data.cms.gov/) |
| **File** | `chronic-disease-data.csv` |
| **Conditions** | 21 chronic conditions |
| **Population** | Medicare fee-for-service beneficiaries |
| **Stratification** | Age group (`<65`, `65+`), Race/Ethnicity (5 groups), Gender, Geography (National + all U.S. states) |

> The raw data file is included in this repository. A data dictionary (`CC_R20_20200630_Data_Dectionary_MCC.pdf`) is also provided in `data/`.

---

## Objectives & Key Findings

### Objective 1: Chronic Disease Prevalence
**Visualization:** Faceted heatmap — prevalence by race/ethnicity and age group, paneled by geography

![Fig1_Heatmap](/chronic-disease-analysis/Results/Fig1_FacetedHeatmap_Race-Age.png)

- Hypertension and Hyperlipidemia are consistently the highest-prevalence conditions in the 65+ group across all geographies, with non-Hispanic Black beneficiaries showing a disproportionately high cardiovascular disease burden.
- In the <65 population, Drug/Substance Abuse, HIV/AIDS, and Schizophrenia show relatively elevated prevalence compared to the 65+ panels, reflecting conditions that disproportionately affect younger Medicare beneficiaries.
- Missing values (white cells) appear most frequently for Native American and Asian/Pacific Islander subgroups in certain states, likely reflecting data suppression due to small population sizes.

---

### Objective 2: Medicare Utilization Analysis
**Visualization:** Grouped bar charts — four utilization metrics by condition and geography (National vs. TX, CA, NY, FL)

![Fig2_BarChart](/chronic-disease-analysis/Results/Fig2_BarChart_UtilizationMetrics.png)

**Spending:**
- Stroke and Heart Failure are the highest-cost conditions under standardized per capita spending; California and New York show elevated actual spending, reflecting regional cost differences that are adjusted away in the standardized metric.
- Autism Spectrum Disorders and Hyperlipidemia rank among the lowest-cost conditions in both spending metrics.

**Hospital Readmissions:**
- HIV/AIDS and Chronic Hepatitis B & C have the highest 30-day readmission rates nationally — notably different from the spending rankings.
- Florida consistently exceeds the national average for readmissions across most conditions; California trends at or below it.

**ER Utilization:**
- Drug/Substance Abuse and Alcohol Abuse drive the highest ER visit rates by a wide margin, far exceeding high-spending conditions like Stroke and Heart Failure.
- New York and Florida show elevated ER utilization relative to the national average for the highest-burden conditions.

---

### Objective 3: State-Wise Prevalence Analysis
**Visualization:** Diverging bar charts — per-condition deviation from national baseline, stratified by gender and race/ethnicity

![Fig3_Gender](/chronic-disease-analysis/Results/Fig3_DivergingBarChart_Gender.png)

![Fig3_Race](/chronic-disease-analysis/Results/Fig4_DivergingBarChart_Race.png)

- For most conditions, male and female beneficiaries show consistent geographic deviation patterns, suggesting that state-level differences in chronic disease burden are generally not driven by one sex.
  - Exceptions: Depression and Osteoporosis show larger female deviations in high-prevalence states; HIV/AIDS shows a substantially larger male deviation in Washington D.C.
- Native American beneficiaries show the greatest geographic variability across conditions, with extreme deviations in Diabetes, COPD, and Chronic Kidney Disease.
- Diabetes shows the strongest geographic concentration: Puerto Rico ranks highest (driven by Hispanic beneficiaries); Vermont ranks lowest.
- Asian/Pacific Islander beneficiaries show the smallest deviations overall, indicating relatively uniform prevalence across states relative to their national baseline.

---

## Repository Structure

```
chronic-disease-analysis/
│
├── README.md
│
├── data/
│   ├── chronic-disease-data.csv
│   └── CC_R20_20200630_Data_Dectionary_MCC.pdf
│
├── results/
│   ├── Fig1_prevalence_heatmap.png
│   ├── Fig2_utilization_barcharts.png
│   └── Fig3_statewise_diverging_bars.png
│
└── Chronic_Disease_Analysis.ipynb
```

---

## Notes

- Some missing values reflect CMS data suppression policies for small subgroup populations and should not be interpreted as true zero prevalence.
- Results should be interpreted in the context of the Medicare fee-for-service population, which may not be representative of the broader U.S. population.

---

## Course Context

**Course:** BMI 5007: Methods in Health Data Science  
The University of Texas Health Science Center at Houston, Spring 2026
