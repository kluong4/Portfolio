# Possum Morphology Regression Analysis

Multiple linear regression analysis of morphometric data from wild mountain brushtail possums (*Trichosurus cunninghami*) to identify predictors of total body length and explore geographic influences on body size variation across seven trapping sites in Australia.

---

## Overview

This project applies regression modeling techniques in R to the `possum` dataset from the DAAG package, which contains measurements from 104 possums trapped at seven sites spanning Southern Victoria to Central Queensland. The analysis investigates how physical traits (skull width, head length, tail length) and geographic site predict total body length, and rigorously evaluates model fit, influential observations, and potential overfitting.

---

## Methods

| Step | Description |
|------|-------------|
| Data Exploration | Correlation matrix and pairwise scatterplots (`ggpairs`) |
| Model Fitting | Multiple linear regression with site (categorical), skull width, head length, tail length |
| Model Interpretation | Coefficient interpretation for numerical and categorical predictors |
| Hypothesis Testing | ANOVA-based nested model comparison, `drop1()` F-tests |
| Interaction Testing | Interaction terms for head length × skull width and head length × tail length |
| Diagnostics | Residuals vs. fitted, Normal Q-Q, Scale-Location, Residuals vs. Leverage plots |
| Influence Analysis | Cook's D via `influencePlot()`, dfbetas for key observations |
| Multicollinearity | Variance Inflation Factors (VIF) |
| Model Selection | Backward stepwise selection (AIC) |
| Overfitting Check | SSE vs. LOO PRESS statistic |

**Tools:** R, `DAAG`, `GGally`, `ggfortify`, `car`

---

## Key Findings

- A multiple linear regression model incorporating site, head length, tail length, and skull width explained **81.2% of variance** in total body length (R² = 0.8118, F = 45.07, p < 2.2×10⁻¹⁶).
- **Geographic site** is a significant predictor of body size (p = 2.5×10⁻¹²), with possums at Bulburin averaging 5.1 cm shorter than the Cambarville reference site after controlling for morphometric variables.
- **Skull width** does not independently contribute to predicting total body length after accounting for head length, tail length, and site (p = 0.21), and was removed by backward model selection.
- Possum **BR1** was identified as the most influential observation via Cook's D and dfbetas analysis, affecting coefficients for tail length and site.
- Evidence of mild **overfitting** was detected (SSE = 360.09 vs. PRESS = 445.29), suggesting the model captures some noise in addition to signal.

---

## Repository Structure

```
PossumMorphology-RegressionAnalysis/
│
├── README.md
├── Possum_RegressionAnalysis_Pipeline.R   — full analysis pipeline
└── RegressionAnalysis_Report.pdf          — written report with figures and interpretation
```

---

## Requirements

**R version:** 4.3+

```r
install.packages(c("GGally", "ggfortify", "car"))

# DAAG contains the possum dataset
install.packages("DAAG")
```

**Dataset:** Built into the `DAAG` package — no external download needed.

```r
library(DAAG)
data(possum)
```

---

## Skills Demonstrated

- Multiple linear regression modeling with numerical and categorical predictors
- Model comparison and variable selection (ANOVA, `drop1()`, backward selection)
- Regression diagnostics (linearity, homoscedasticity, normality, leverage)
- Interaction effect testing and interpretation
- Influential observation detection (Cook's D, dfbetas)
- Multicollinearity assessment (VIF)
- Cross-validation for overfitting detection (PRESS statistic)
- Data visualization in R (`ggplot2`, `GGally`, `ggfortify`)

---

## Course Context

**SDS 324E:** Elements of Regression Analysis

The University of Texas at Austin

Fall 2024
