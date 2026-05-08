# =============================================================================
# Script:    Linear Regression Analysis Pipeline
# Project:   Linear Regression Analysis of Train & Engineer Workers
# Course:    SDS 320E: Elements of Statistics
#            UT Austin, Fall 2022
# Input:     Data/TE_Worker_Data.xls  (from https://catalog.data.gov)
# =============================================================================

# ----- Load Libraries --------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)

# ----- Set Working Directory -------------------------------------------------
setwd("~/Desktop/T&Eworkers-RegressionAnalysis")

# ----- Load Data -------------------------------------------------------------
TE_Worker_Data <- read_excel("Data/TE_Worker_Data.xlsx", sheet = "Background")

# Extract Data
clean_data <- na.omit(TE_Worker_Data[, c("Sick days", "Schedule", "Total years of service")])
clean_data$`Schedule` <- factor(clean_data$`Schedule`, labels = c('fixed', 'variable'))
names(clean_data)=c("sick_days", "schedule", "service")

# Select 80 observations for each schedule group (n = 160)
mydata <- clean_data %>% group_by(schedule) %>% slice_sample(n=80)

# ----- Descriptive Analysis of Response Variable -----------------------------
boxplot(
  mydata$service, 
  data = mydata, 
  ylab='Frequency', 
  main='Distribution of Total Years of Service', 
  col='light blue')

summary(mydata$service)
IQR(mydata$service)

# test transformations - log()
boxplot(
  log(mydata$service), 
  data = mydata, 
  ylab='Frequency', 
  main='Distribution of Total Years of log(Service)', 
  col='light blue')

# test transformations - sqrt()
boxplot(
  sqrt(mydata$service), 
  data = mydata, 
  ylab='Frequency', 
  main='Distribution of the Square Root of Total Years of Service', 
  col='light blue')

# test transformations - Cube Root
boxplot(
  '^'(mydata$service, 1/3), 
  data = mydata, 
  ylab='Frequency', 
  main='Distribution of the Cube Root of Total Years of Service', 
  col='light blue')

# ----- Investigation of Numeric Explanatory Variable -------------------------
boxplot(
  mydata$sick_days, 
  data = ydata, 
  ylab='Frequency', 
  main='Distribution of Sick Days', 
  col='light pink')

summary(mydata$sick_days)
IQR(mydata$sick_days)

# linearity assumption - Scatterplot of Sick Days and Service
plot(
  service~sick_days, 
  data = mydata, 
  main='Service vs Sick Days', 
  pch=20)

# test transformations - log()
boxplot(
  log(mydata$sick_days), 
  data = mydata, 
  ylab='Frequency', 
  main='Distribution of log(Sick Days)', 
  col='light pink')

# test transformations - sqrt()
boxplot(
  sqrt(mydata$sick_days), 
  data = mydata, 
  ylab='Frequency', 
  main='Distribution of the Square Root of Sick Days', 
  col='light pink')

# test transformations - Cube Root
boxplot(
  '^'(mydata$sick_days, 1/3), 
  data = mydata, 
  ylab='Frequency', 
  main='Distribution of the Cube Root of Sick Days', 
  col='light pink')

# ----- Investigation of Categorical Explanatory Variable --------------------
boxplot(
  service~schedule, 
  data = mydata, 
  main = 'Service vs Schedule', 
  col = c('light blue', 'light green'))

# ----- Linear Regression Analysis --------------------------------------------

# Center sick_days variable
mydata$sick_days_c <- mydata$sick_days - mean(mydata$sick_days)

# Fit linear regression model with interaction
model1 <- lm(
  service ~ schedule * sick_days_c,
  data = mydata
)

# View model summary
summary(model1)


# ----- Assumptions -----------------------------------------------------------
# Normality of the Residual Assumption - Histogram of Residuals
hist(
  residuals(model1),
  main = "Model Residuals",
  xlab = "Residuals",
  ylab = "Frequency",
  col = "light blue",
  border = "black"
)

# Equal Variance Assumption - Residual Plot
plot(
  fitted(model1),
  residuals(model1),
  xlab = "Fitted Values",
  ylab = "Residuals",
  main = "Residual Plot"
)
abline(h = 0, col = "red")

# ----- Interaction Plot ------------------------------------------------------
ggplot(
  mydata,
  aes(
    x = sick_days_c,
    y = service,
    color = schedule
  )
) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", aes(fill = schedule)) +
  labs(
    title = "Service by Sick Days and Schedule",
    x = "Sick Days (days)",
    y = "Total Years of Service (years)",
    color = "Schedule",
    fill = "Schedule"
  ) +
  theme_classic()
