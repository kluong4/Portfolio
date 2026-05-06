# =============================================================================
# Script:    Regression Analysis Pipeline
# Project:   Uncovering Patterns in Possum Morphology
# Course:    SDS 324E: Elements of Regression Analysis
#            UT Austin, Fall 2024
# Hypotheses:
#           1) Including site significantly improves the model’s ability to predict total body length in possums.
#           2) Skull width does not significantly contribute to the prediction of total body length.
#           3) Head length affects total length depending on the tail length, and vice versa.
# =============================================================================

# ----- Load Libraries --------------------------------------------------------
library(DAAG)
library(GGally)
library(ggfortify)
library(car)

# ----- Set Working Directory -------------------------------------------------
setwd("/Volumes/DONGHOA/Kim's UT 25 BIO Undergrad/4 Fall 2024/SDS 324E/Project")

# ----- Load Clean Data -------------------------------------------------------
data(possum) # help(possum)

possum_data <- possum[, c("site", "skullw", "hdlngth", "taill", "totlngth")]
possum_data$site <- as.character(possum_data$site)

# ----- Exploring Marginal Relationships in the Data --------------------------
cor_analysis <- ggpairs(possum_data, lower = list(combo = "box_no_facet"))
cor_analysis

# ----- Multiple Linear Regression model --------------------------------------
possum_lm <- lm(totlngth ~ skullw + hdlngth + taill + site, data = possum_data)

# mlr results
summary(possum_lm) # Fig2

# residual standard error
summary(possum_data$totlngth)
1.96/75*100
1.96/96.5*100

# ----- Improving the Model ---------------------------------------------------
# diagnostic plot
autoplot(possum_lm) 

# no need for transformations

# Interaction Term Significance
possum_interaction <- lm(totlngth ~ hdlngth * skullw + taill + site, data = possum_data)
summary(possum_interaction)

anova(possum_lm, possum_interaction)

# ----- Hypothesis 1 ----------------------------------------------------------
# Reduced model without the site variable
reduced_model <- lm(totlngth ~ hdlngth + taill + skullw, data = possum_data)
# Perform an ANOVA to compare the two models
anova(reduced_model, possum_lm)

# ----- Hypothesis 2 ----------------------------------------------------------
drop1(possum_lm, test="F")

# ----- Hypothesis 3 ----------------------------------------------------------
# Interaction Term Significance
possum_interaction_hdlngth_taill <- lm(totlngth ~ hdlngth * taill + skullw + site, data = possum_data)
summary(possum_interaction_hdlngth_taill)

# ----- Robustness of Results -------------------------------------------------
### Influential Observations:
autoplot(possum_lm, which = 5) # from diagnostic plot
# highest leverage: BR1, BSF12, WW1

influencePlot(possum_lm_squared) # also show Cook's distance

possum_dfbetas <- dfbetas(possum_lm)
threshold <- 3 / sqrt(nrow(possum_data))
subset(possum_data, abs(possum_dfbetas[,"taill"]) > threshold)

subset(possum_data, abs(possum_dfbetas[,"site4"]) > threshold)

### Multicolinearity
vif(possum_lm)

### LOO Prediction Error and Overfitting:
press(possum_lm) #PRESS
sum((possum_lm$residuals)^2) #SSE

# ----- Automatic Selection --------------------------------------------------
possum_step <- step(possum_lm, direction = "backward")