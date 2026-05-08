# =============================================================================
# Script:    Northern Shoveler Foraging Behavior Analysis
# Project:   Behavioral Ecology of Northern Shovelers
# Course:    BIO 359K: Animal Behavior
#            UT Austin, Spring 2023
# Input:     Data/dataset.csv
# Output:    Results/Fig1_BehaviorHistograms.png
#            Results/Fig2_DistributionBehaviorSex.png
#            Results/Fig3_DistributionBehaviorTime.png
#            Results/Table1_ShapiroTestResults.csv
#            Results/Table2_SexSummary.csv
#            Results/Table3_TimeSummary.csv
#            Results/Table4_PairedTests.csv
# =============================================================================

# ----- Load Libraries --------------------------------------------------------
library(ggplot2)
library(reshape2)
library(patchwork)
library(tidyr)
library(dplyr)

# ----- Set Working Directory -------------------------------------------------
setwd("~/Desktop/NorthernShovelers-BehaviorAnalysis")

# ----- Load Data -------------------------------------------------------------
mydata <- read.csv("Data/dataset.csv")
p1 <- mydata

# Clean inconsistent time labeling (fix trailing space issue)
p1[p1$Time.of.Day. == "evening ", ]$Time.of.Day. <- "evening"


# =============================================================================
# ----- Behavior - Normality Results ------------------------------------------
# =============================================================================

# ----- Create Total Behavior Counts -----------------------------------------
pecking_tot  <- p1$Pecking..f. + p1$Pecking..m.
dabbling_tot <- p1$Dabbling..f. + p1$Dabbling..m.
dipping_tot  <- p1$Dipping..f. + p1$Dipping..m.
upending_tot <- p1$Upending..f. + p1$Upending..m.

p1$pecking_tot  <- pecking_tot
p1$dabbling_tot <- dabbling_tot
p1$dipping_tot  <- dipping_tot
p1$upending_tot <- upending_tot


# ----- Figure 1: Histograms of Behaviors -------------------------------------
dist_pecking <- ggplot(data.frame(pecking_tot), aes(x = pecking_tot)) +
  geom_histogram(binwidth = 1, boundary = 0, fill = "coral", color = "black") +
  labs(x = "Pecking Frequency", title = "Pecking Distribution") +
  xlim(0, 10) +
  theme_classic()

dist_dipping <- ggplot(data.frame(dipping_tot), aes(x = dipping_tot)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(x = "Dipping Frequency", title = "Dipping Distribution") +
  theme_classic()

dist_dabbling <- ggplot(data.frame(dabbling_tot), aes(x = dabbling_tot)) +
  geom_histogram(bins = 18, fill = "seagreen", color = "black") +
  labs(x = "Dabbling Frequency", title = "Dabbling Distribution") +
  xlim(0, 20) +
  theme_classic()

dist_upending <- ggplot(data.frame(upending_tot), aes(x = upending_tot)) +
  geom_histogram(bins = 12, fill = "lavender", color = "black") +
  labs(x = "Upending Frequency", title = "Upending Distribution") +
  theme_classic()

combined_distributions <- (dist_pecking + dist_dipping) / (dist_dabbling + dist_upending)

ggsave("Results/Fig1_BehaviorHistograms.png", plot = combined_distributions, width = 12, height = 10, dpi = 700)


# ----- Table 1: Shapiro-Wilk Normality Tests --------------------------------
shapiro_results <- data.frame(
  Behavior = c("Pecking", "Dipping", "Dabbling", "Upending"),
  W_Statistic = c(
    shapiro.test(pecking_tot)$statistic,
    shapiro.test(dipping_tot)$statistic,
    shapiro.test(dabbling_tot)$statistic,
    shapiro.test(upending_tot)$statistic
  ),
  P_Value = c(
    shapiro.test(pecking_tot)$p.value,
    shapiro.test(dipping_tot)$p.value,
    shapiro.test(dabbling_tot)$p.value,
    shapiro.test(upending_tot)$p.value
  )
)

write.csv(shapiro_results, "Results/Table1_ShapiroTestResults.csv", row.names = FALSE)



# =============================================================================
# ----- Sex - Behavior Analysis -----------------------------------------------
# =============================================================================

# ----- Create Sex DataFrame --------------------------------------------------
male_df <- data.frame(
  freq = c(p1$Pecking..m., p1$Dabbling..m., p1$Dipping..m., p1$Upending..m.),
  behavior = rep(c("pecking","dabbling","dipping","upending"),
                 each = length(p1$Pecking..m.)),
  sex = "male"
)

female_df <- data.frame(
  freq = c(p1$Pecking..f., p1$Dabbling..f., p1$Dipping..f., p1$Upending..f.),
  behavior = rep(c("pecking","dabbling","dipping","upending"),
                 each = length(p1$Pecking..f.)),
  sex = "female"
)

sex_df <- rbind(male_df, female_df)

# ----- Table 2: Sex-Based Summary Statistics ---------------------------------
sex_summary_stats <- sex_df %>%
  group_by(behavior, sex) %>%
  summarise(
    mean = mean(freq, na.rm = TRUE),
    sd   = sd(freq, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(sex_summary_stats, "Results/Table2_SexSummary.csv", row.names = FALSE)

# ----- Figure 2: Sex Distribution Boxplot ------------------------------------
sex_df_updated <- sex_df %>%
  group_by(behavior, sex) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = sex,
              values_from = freq,
              id_cols = c(row_id, behavior)) %>%
  select(-row_id)

sex_df.melt <- melt(sex_df_updated, id.var = "behavior")
sex_df.melt$behavior <- factor(
  sex_df.melt$behavior,
  levels = c("pecking","dabbling","dipping","upending")
)

boxplot_sex <- ggplot(sex_df.melt, aes(x = variable, y = value, color = behavior)) +
  geom_boxplot() +
  labs(x = "Sex",
       y = "Frequency",
       title = "Foraging Behavior by Sex",
       color = "Behavior") +
  theme_minimal()

ggsave("Results/Fig2_DistributionBehaviorSex.png", boxplot_sex, width = 8, height = 6, dpi = 700)



# =============================================================================
# ----- Time - Behavior Analysis ----------------------------------------------
# =============================================================================

# ----- Create Time of Day DataFrame ------------------------------------------
time_df <- data.frame(
  freq = c(pecking_tot, dabbling_tot, dipping_tot, upending_tot),
  behavior = rep(c("pecking","dabbling","dipping","upending"),
                 each = length(pecking_tot)),
  time = p1$Time.of.Day.
)

# ----- Table 3: Time-of-Day Summary ------------------------------------------
time_summary_stats <- time_df %>%
  group_by(behavior, time) %>%
  summarise(
    mean = mean(freq, na.rm = TRUE),
    sd   = sd(freq, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(time_summary_stats, "Results/Table3_TimeSummary.csv", row.names = FALSE)

# ----- Figure 3: Time-of-Day Boxplot -----------------------------------------
time_df_updated <- time_df %>%
  group_by(behavior, time) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = time,
              values_from = freq,
              id_cols = c(row_id, behavior)) %>%
  select(-row_id)

time_df.melt <- melt(time_df_updated, id.var = "behavior")
time_df.melt$behavior <- factor(time_df.melt$behavior,
                                levels = c("pecking","dabbling","dipping","upending"))
time_df.melt$variable <- factor(time_df.melt$variable,
                                levels = c("morning","afternoon","evening"))

boxplot_time <- ggplot(time_df.melt,
                       aes(x = variable, y = value, color = behavior)) +
  geom_boxplot() +
  labs(x = "Time of Day",
       y = "Frequency",
       title = "Behavior Across Time",
       color = "Behavior") +
  theme_minimal()

ggsave("Results/Fig3_DistributionBehaviorTime.png", boxplot_time, width = 8, height = 6, dpi = 700)



# =============================================================================
# ----- Table 4: Paired t-tests ----------------------------------------------
# =============================================================================
comparisons <- list(
  c("pecking","dabbling"),
  c("pecking","dipping"),
  c("pecking","upending"),
  c("dabbling","dipping"),
  c("dabbling","upending"),
  c("dipping","upending")
)

paired_test_df <- do.call(rbind, lapply(comparisons, function(x) {
  test <- t.test(
    get(paste0(x[1], "_tot")),
    get(paste0(x[2], "_tot")),
    paired = TRUE
  )
  data.frame(
    behavior1 = x[1],
    behavior2 = x[2],
    t = as.numeric(test$statistic),
    df = as.numeric(test$parameter),
    p_value = test$p.value,
    conf_low = test$conf.int[1],
    conf_high = test$conf.int[2],
    mean_difference = as.numeric(test$estimate)
  )
}))

write.csv(paired_test_df,"Results/Table4_PairedTests", row.names = FALSE)
