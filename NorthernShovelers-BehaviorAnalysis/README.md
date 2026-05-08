# Northern Shoveler Foraging Behavior — Field Study & Statistical Analysis
Scan-sampling field study of foraging time budgets in Northern Shovelers (*Anas clypeata*) at Mueller Lake Park in Austin, Texas, to assess whether time of day and sex predict foraging behavior frequency, and to evaluate the influence of human presence on duck activity patterns in an urban park setting.

---

## Overview
This project applies non-parametric hypothesis testing in R to behavioral scan-sampling data collected from Northern Shovelers across nine observation sessions spanning morning, afternoon, and evening time blocks. The analysis investigates whether foraging behavior frequency — across four distinct behaviors (pecking, dabbling, dipping, upending) — differs significantly by time of day or sex, and situates findings within the context of established literature on shoveler time budgets and captive waterfowl behavior.

---

## Methods

| Step | Description |
|------|-------------|
| Ethogram Design | Defined 4 state foraging behaviors with precise operational criteria |
| Data Collection | Scan sampling at 10-minute intervals; 6 scans per observation hour |
| Normality Testing | Shapiro-Wilk test applied to all four behavior frequency distributions |
| Sex Comparison | Wilcoxon-Mann Whitney test (two independent non-normal groups) |
| Time of Day Comparison | Kruskal-Wallis test (three independent non-normal groups) |
| Behavior Comparisons | Pairwise paired t-tests across all 6 foraging behavior combinations |
| Visualization | Box-and-whisker plots and histograms by sex and time of day |
| Environmental Controls | Hourly temperature, human presence, and species co-occurrence logged |

**Tools:** R, `ggplot2`, Google Sheets, Microsoft PowerPoint, Microsoft Word

---

## Key Findings
- All four foraging behavior distributions failed the Shapiro-Wilk normality test (all p < 2.2×10⁻¹⁶), with right-skewed histograms, justifying the use of non-parametric tests throughout.
- **Sex** was not a significant predictor of foraging frequency (Wilcoxon W = 134,036, p = 0.10), consistent with the hypothesis that human presence in an urban park offsets the natural tendency for females to forage more than males (Afton, 1979; Rose et al., 2022).
- **Time of day** was a highly significant predictor of foraging frequency (Kruskal-Wallis χ² = 150.26, df = 2, p < 2.2×10⁻¹⁶), with foraging increasing markedly from morning to evening across all four behaviors.
- **Dabbling and dipping** were statistically equivalent in frequency (paired t-test p = 0.987) and collectively formed the dominant foraging strategy; **pecking and upending** were both significantly less frequent (p < 0.001).
- Morning inactivity (sleeping behavior) was identified as an unanticipated **confounding variable** contributing to low morning foraging rates, independent of human presence effects.

---

## Repository Structure
```
NorthernShoveler-ForagingBehavior/
│
├── README.md
├── Data/
│   └── dataset.csv                          — raw scan-sampling field observations
│
├── Results/
│   ├── Fig1_BehaviorHistograms.png          — distribution plots for all four foraging behaviors
│   ├── Fig2_DistributionBehaviorSex.png     — foraging frequency by Sex (Male, Female)
│   ├── Fig3_DistributionBehaviorTime.png    - foraging frequency by Time of Day (morning, afternoon, evening)
│   ├── Table1_ShapiroTestResults.csv        — Shapiro-Wilk's Test Results for Foraging behaviors
│   ├── Table2_SexSummary.csv                — Wilcoxon Test Results for Sex-Foraging Behaviors
│   ├── Table3_TimeSummary.csv               - Kruskal-Wallis Test Results for Time of Day-Foraging Behaviors
│   └── Table4_PairedTests.csv               - Paired t-tests Results Between Foraging Behaviors
│
├── Presentation/
│   └── TimeBudgetOfFeedingInNorthernShovelers_Presentation.pdf
│
└── Script/
    └── NorthernShovelers_Analysis.R         — full statistical analaysis R script
```

---

## Requirements
**R version:** 4.3+
```r
install.packages(c("ggplot2", "dplyr"))
```

**Dataset:** Collected via original field observation at Mueller Lake Park, Austin, TX (Feb 23 – Mar 6).  
No external download needed — raw data is included in `Data/Datasheet.csv`.

```r
shoveler_data <- read.csv("Data/Datasheet.csv")
```

---

## Skills Demonstrated

**Technical**
- Non-parametric hypothesis testing (Wilcoxon-Mann Whitney, Kruskal-Wallis) with correct test selection justified by normality diagnostics
- Pairwise behavioral comparisons using paired t-tests with effect direction interpretation
- Ethogram design with operationalized, reproducible behavioral definitions
- Scan-sampling data collection with multi-observer coordination and standardized recording
- Data visualization in R (`ggplot2`) — histograms, box-and-whisker plots
- Environmental covariate tracking and integration alongside behavioral data
- Scientific report writing in IMRaD format with primary literature citation

**Soft**
- Coordinated a 3-person field research team across multiple weeks with variable environmental conditions and shoveler population counts ranging from 0 to ~100 individuals
- Adapted observation schedules in real time in response to weather disruptions and low-count sessions
- Transparently reported unexpected confounding variables (morning inactivity, cold-front effects) and revised interpretations accordingly rather than minimizing inconvenient findings
- Communicated statistical results and conservation implications to a non-specialist audience in a formal class presentation

---

## Course Context
**BIO 359K:** Animal Behavior 

The University of Texas at Austin  

Spring 2023
