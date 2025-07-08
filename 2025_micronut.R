# ===============
# Micronutrients
# ===============

library(ggstatsplot)
library(ggplot2)
library(phyloseq)
library(phyloseqCompanion)
library(performance)
library(lme4)
library(glmmTMB)
library(nortest)
library(MASS)
library(dplyr)
library(purrr)
library(tibble)
library(ggpubr)


micronut_parrot <- read.csv("##",header=T)
micronut_parrot_sub <- micronut_parrot %>%
  dplyr::select(sample, As, Ca, K, Cr, Fe, Zn, Se, Mg, P, Mn, site_condition)

# transform data into mg/100g or ug/100g -- data are usually given as ug/g
# transform as wet weight
percentWater <- 77.5 

dminerals <- micronut_parrot_sub %>%
  mutate(Ca = (Ca * ((100 - percentWater)/100))/10, #expressed as mg/100g
          K = (K * ((100 - percentWater)/100))/10, #expressed as mg/100g
          Mg = (Mg * ((100 - percentWater)/100))/10, #expressed as mg/100g
          P = (P * ((100 - percentWater)/100))/10, #expressed as mg/100g
          Mn = (Mn * ((100 - percentWater)/100))/10, #expressed as mg/100g
          Se = (Se * ((100 - percentWater)/100))*100, #expressed as ug/100g
          Zn = (Zn * ((100 - percentWater)/100))/10, #expressed as mg/100g
          Fe = (Fe * ((100 - percentWater)/100))/10, #expressed as mg/100g
          Cr = (Cr * ((100 - percentWater)/100))/10) #expressed as mg/100g

popfish_rarefied_2024_fish <-readRDS("##")

meta_popfish <- sample.data.frame(popfish_rarefied_2024_fish)
merge_popfish_micronut <- merge(meta_popfish, dminerals, by.x="sample",by.y="sample", no.dups=T)
merge_popfish_micronut <- as_tibble(merge_popfish_micronut) # so it prints a little nicer
merge_popfish_micronut <- merge_popfish_micronut %>%
  rename(site_condition = site_condition.x)
merge_popfish_micronut$SL <- as.numeric(gsub(",", ".", merge_popfish_micronut$SL))

# ===============
# Potassium
# ===============
K_parrot_m1 <- glmer(K ~ site_condition + (1|site_name), data = merge_popfish_micronut, family=gaussian(link="log"))
K_parrot_m2 <- lmer(K ~ site_condition + (1|site_name), data = merge_popfish_micronut)

AIC(K_parrot_m1, K_parrot_m2)

testDispersion(K_parrot_m1)
simulationOutput <- simulateResiduals(fittedModel = K_parrot_m1, plot = F)
plot(simulationOutput)
check_model(K_parrot_m1)
check_normality(K_parrot_m2)
ad.test(residuals(K_parrot_m1)) ## ok
ks.test(residuals(K_parrot_m1), "pnorm", mean = mean(residuals(K_parrot_m1)), sd = sd(residuals(K_parrot_m1)))
leveneTest(K ~ site_condition, merge_popfish_micronut)

library(multcomp)
summary(glht(K_parrot_m1, linfct = mcp(site_condition = "Tukey")), test = adjusted("holm"))

# Site effect
K_parrot_site_m1 <- kruskal.test(K ~ site_name, data=merge_popfish_micronut)
dunn_K_site <- dunn.test(merge_popfish_micronut$K, merge_popfish_micronut$site_name, list=T, method="bh")

variability <- merge_popfish_micronut %>%
  group_by(site_condition) %>%
  summarise(
    Mean = mean(K),
    SD = sd(K),
    Variance = var(K),
    CV = (SD / Mean) * 100  
  )

diff_percent_table <- variability %>%
  select(site_condition, Mean) %>%
  rename(Mean1 = Mean, Condition1 = site_condition) %>%
  crossing(
    variability %>% select(site_condition, Mean) %>% 
      rename(Mean2 = Mean, Condition2 = site_condition)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((Mean1 - Mean2) / Mean2) * 100
  )

# Calculate variance (pairwise differences for a condition)
pairwise_diff <- function(data) {
  combs <- combn(data$K, 2)  
  abs_diff <- abs(combs[1, ] - combs[2, ])  
  return(abs_diff)
}

# apply to each condition
B_diff <- pairwise_diff(subset(merge_popfish_micronut, site_condition == "B"))
LC_diff <- pairwise_diff(subset(merge_popfish_micronut, site_condition == "LC"))
LR_diff <- pairwise_diff(subset(merge_popfish_micronut, site_condition == "LR"))

all_diffs <- data.frame(
  diff = c(B_diff, LC_diff, LR_diff),
  condition = factor(rep(c("B", "LC", "LR"), times = c(length(B_diff), length(LC_diff), length(LR_diff))))
)

ggplot(all_diffs, aes(x = condition, y = diff)) +
  geom_boxplot() +
  labs(title = "comparisons of pairwise differences", y = "pairwise Differences", x = "condition") +
  theme_minimal()

# calculate variance for each condition
variance_b <- mean(all_diffs$diff[all_diffs$condition == "B"])
variance_LC <- mean(all_diffs$diff[all_diffs$condition == "LC"])
variance_LR <- mean(all_diffs$diff[all_diffs$condition == "LR"])

variability_table <- tibble(
  condition = c("B", "LC", "LR"),
  mean_variability = c(variance_b, variance_LC, variance_LR)
)

var_diff_table <- variability_table %>%
  rename(Condition1 = condition, MeanVar1 = mean_variability) %>%
  crossing(
    variability_table %>% rename(Condition2 = condition, MeanVar2 = mean_variability)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((MeanVar1 - MeanVar2) / MeanVar2) * 100
  )

# Perform Kruskal-Wallis test on the pairwise differences
kruskal.test(diff ~ condition, data = all_diffs)
dunn.test(all_diffs$dif, all_diffs$condition, list=T, method="bh")

# ===========
# Calcium
# ===========
remove_outliers <- function(df, column) {
  Q1 <- quantile(df[[column]], 0.25)
  Q3 <- quantile(df[[column]], 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  df %>% filter(df[[column]] >= lower_bound & df[[column]] <= upper_bound)
}

Ca_woutliers <- merge_popfish_micronut %>%
  do(remove_outliers(., "Ca")) %>%
  ungroup()

b <- boxcox(lm(Ca_woutliers$Ca ~ 1))
lambda <- b$x[which.max(b$y)]
lambda
new_x_exact <- (Ca_woutliers$Ca ^ lambda - 1) / lambda
Ca_woutliers$Ca_box <- new_x_exact

Ca_parrot_v1 <- lmer(Ca_box ~ site_condition + (1|site_name), data =Ca_woutliers)
Ca_parrot <- glmer(Ca_box ~ site_condition + (1|site_name), data =Ca_woutliers, gaussian(link="log"))
Ca_parrot_v2 <- glmer(Ca_box ~ site_condition + (1|site_name), data = Ca_woutliers, Gamma(link = "inverse"))

Ca_parrot_v1 <- lmer(Ca ~ site_condition + (1|site_name), data =Ca_woutliers)
Ca_parrot <- glmer(Ca ~ site_condition + (1|site_name), data =Ca_woutliers, gaussian(link="log"))
Ca_parrot_v2 <- glmer(Ca ~ site_condition + (1|site_name), data = Ca_woutliers, Gamma(link = "inverse"))

sim_residuals <- simulateResiduals(Ca_parrot_v1)
plot(sim_residuals)
testUniformity(sim_residuals)    
testDispersion(sim_residuals)   
testZeroInflation(sim_residuals) 

check_model(Ca_parrot_v1)
check_heteroscedasticity(Ca_parrot)

testDispersion(Ca_parrot_v1)
simulationOutput <- simulateResiduals(fittedModel =Ca_parrot_v1, plot = F)
plot(simulationOutput)

shapiro.test(residuals(Ca_parrot_v1))
ad.test(residuals(Ca_parrot_v1)) 
ks.test(residuals(Ca_parrot_v1), "pnorm", mean = mean(residuals(Ca_parrot_v1)), sd = sd(residuals(Ca_parrot_v1)))
leveneTest(Ca_box ~ site_condition, Ca_woutliers)

variability_Ca <- Ca_woutliers %>%
  group_by(site_name) %>%
  summarise(
    Mean = mean(Ca),
    SD = sd(Ca),
    Variance = var(Ca),
    Median = median(Ca),
    CV = (SD / Mean) * 100  
  )

site_median <- variability_Ca %>%
  select(site_name,Median) %>%
  distinct()

diff_percent_table <- site_means %>%
  rename(Site1 = site_name, Mean1 = Mean) %>%
  crossing(
    site_means %>%
      rename(Site2 = site_name, Mean2 = Mean)
  ) %>%
  filter(Site1 < Site2) %>%
  mutate(
    Percent_Diff = ((Mean1 - Mean2) / Mean2) * 100
  )

# Results by Site
Ca_woutliers <- merge_popfish_micronut %>%
  group_by(site_name) %>%
  do(remove_outliers(., "Ca")) %>%
  ungroup()

b <- boxcox(lm(Ca_woutliers$Ca ~ 1))
lambda <- b$x[which.max(b$y)]
lambda
new_x_exact <- (Ca_woutliers$Ca ^ lambda - 1) / lambda
Ca_woutliers$Ca_box <- new_x_exact

Ca_parrot_site_m1 <- lm(Ca_box ~ site_name, data =Ca_woutliers)
Ca_parrot_site_m1 <- aov((Ca_box) ~ site_name, data =Ca_woutliers)
check_model(Ca_parrot_site_m1)
check_heteroscedasticity((Ca_parrot_site_m1))

TukeyHSD(Ca_parrot_site_m1, "site_name", ordered = TRUE)

# Create a function to calculate pairwise differences for a condition
pairwise_diff <- function(data) {
  combs <- combn(data$Ca, 2)  
  abs_diff <- abs(combs[1, ] - combs[2, ])  
  return(abs_diff)
}

# Apply the function to each condition (B, LC, LR)
B_diff <- pairwise_diff(subset(Ca_woutliers, site_condition == "B"))
LC_diff <- pairwise_diff(subset(Ca_woutliers, site_condition == "LC"))
LR_diff <- pairwise_diff(subset(Ca_woutliers, site_condition == "LR"))

all_diffs <- data.frame(
  diff = c(B_diff, LC_diff, LR_diff),
  condition = factor(rep(c("B", "LC", "LR"), times = c(length(B_diff), length(LC_diff), length(LR_diff))))
)

# calculate variance for each condition
variance_b <- mean(all_diffs$diff[all_diffs$condition == "B"])
variance_LC <- mean(all_diffs$diff[all_diffs$condition == "LC"])
variance_LR <- mean(all_diffs$diff[all_diffs$condition == "LR"])

variability_table <- tibble(
  condition = c("B", "LC", "LR"),
  mean_variability = c(variance_b, variance_LC, variance_LR)
)

var_diff_table <- variability_table %>%
  rename(Condition1 = condition, MeanVar1 = mean_variability) %>%
  crossing(
    variability_table %>% rename(Condition2 = condition, MeanVar2 = mean_variability)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((MeanVar1 - MeanVar2) / MeanVar2) * 100
  )

kruskal.test(diff ~ condition, data = all_diffs)
dunn.test(all_diffs$dif, all_diffs$condition, list=T, method="bh")

# Computing by site_name
all_diffs_list <- Ca_woutliers %>%
  split(.$site_name) %>%
  map(~ pairwise_diff(.x))  

# Combine into one dataframe
all_diffs <- tibble(
  site_name = rep(names(all_diffs_list), times = map_int(all_diffs_list, length)),
  diff = unlist(all_diffs_list)
)

variability_table <- all_diffs %>%
  group_by(site_name) %>%
  summarise(mean_variability = mean(diff))

var_diff_table <- variability_table %>%
  rename(Site1 = site_name, MeanVar1 = mean_variability) %>%
  crossing(
    variability_table %>% rename(Site2 = site_name, MeanVar2 = mean_variability)
  ) %>%
  filter(Site1 != Site2) %>%
  mutate(
    Percent_Diff = ((MeanVar1 - MeanVar2) / MeanVar2) * 100
  )

var_diff_table_clean <- var_diff_table %>%
  filter(Site1 < Site2)

kruskal.test(diff ~ site_name, data = all_diffs)
dunn.test(all_diffs$diff, all_diffs$site_name, method = "bh", list=T)

# ===========
# Chromium
# ===========
Cr_woutliers <- merge_popfish_micronut %>%
  group_by(site_name) %>%
  do(remove_outliers(., "Cr")) %>%
  ungroup()

Cr_woutliers_v2 <- Cr_woutliers[-33,]

Cr_parrot_v1 <- lmer(Cr ~ site_condition + (1|site_name), data = Cr_woutliers)
Cr_parrot <- glmer(Cr ~ site_condition + (1|site_name), data = Cr_woutliers, gaussian(link="log"))
Cr_parrot_v2 <- glmer(Cr ~ site_condition + (1|site_name), data = Cr_woutliers, Gamma(link = "log"))

AIC(Cr_parrot, Cr_parrot_v2, Cr_parrot_v1)

check_model(Cr_parrot_v2)
check_heteroscedasticity(Cr_parrot)

shapiro.test(residuals(Cr_parrot_v2))
ad.test(residuals(Cr_parrot_v2)) 
ks.test(residuals(Cr_parrot_v2), "pnorm", mean = mean(residuals(Cr_parrot_v2)), sd = sd(residuals(Cr_parrot_v2)))
leveneTest(Cr ~ site_condition, Cr_woutliers_v2)

variability_cr <- Cr_woutliers %>%
  group_by(site_condition) %>%
  summarise(
    Mean = mean(Cr),
    SD = sd(Cr),
    Variance = var(Cr),
    CV = (SD / Mean) * 100  
  )

diff_percent_table_cr <- variability_cr %>%
  select(site_condition, Mean) %>%
  rename(Mean1 = Mean, Condition1 = site_condition) %>%
  crossing(
    variability_cr %>% select(site_condition, Mean) %>% 
      rename(Mean2 = Mean, Condition2 = site_condition)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((Mean1 - Mean2) / Mean2) * 100
  )

summary(glht(Cr_parrot_v2, linfct = mcp(site_condition = "Tukey")), test = adjusted("fdr"))

# Effect of site
Cr_parrot_site_m1 <- lm(log10(Cr) ~ site_name, data =Cr_woutliers)
Cr_parrot_site_m1 <- aov((Cr) ~ site_name, data =Cr_woutliers)
check_model(Cr_parrot_site_m1)
check_heteroscedasticity((Cr_parrot_site_m1))
Cr_parrot_site_m1 <- kruskal.test(Cr ~ site_name, data=Cr_woutliers) #NS


pairwise_diff <- function(data) {
  combs <- combn(data$Cr, 2)  
  abs_diff <- abs(combs[1, ] - combs[2, ])  
  return(abs_diff)
}

# Apply the function to each condition (B, LC, LR)
B_diff <- pairwise_diff(subset(Cr_woutliers, site_condition == "B"))
LC_diff <- pairwise_diff(subset(Cr_woutliers, site_condition == "LC"))
LR_diff <- pairwise_diff(subset(Cr_woutliers, site_condition == "LR"))

all_diffs <- data.frame(
  diff = c(B_diff, LC_diff, LR_diff),
  condition = factor(rep(c("B", "LC", "LR"), times = c(length(B_diff), length(LC_diff), length(LR_diff))))
)

ggplot(all_diffs, aes(x = condition, y = diff)) +
  geom_boxplot() +
  labs(title = "Comparison of Pairwise Differences", y = "Pairwise Differences", x = "Condition") +
  theme_minimal()

# Calculate variance for each condition
variance_b <- mean(all_diffs$diff[all_diffs$condition == "B"])
variance_LC <- mean(all_diffs$diff[all_diffs$condition == "LC"])
variance_LR <- mean(all_diffs$diff[all_diffs$condition == "LR"])

variability_table <- tibble(
  condition = c("B", "LC", "LR"),
  mean_variability = c(variance_b, variance_LC, variance_LR)
)

var_diff_table <- variability_table %>%
  rename(Condition1 = condition, MeanVar1 = mean_variability) %>%
  crossing(
    variability_table %>% rename(Condition2 = condition, MeanVar2 = mean_variability)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((MeanVar1 - MeanVar2) / MeanVar2) * 100
  )

var_diff_table_clean <- var_diff_table %>%
  filter(Condition1 < Condition2)

kruskal.test(diff ~ condition, data = all_diffs)
dunn.test(all_diffs$dif, all_diffs$condition, list=T, method="bh")


# =========
# Iron
# =========
Fe_woutliers <- merge_popfish_micronut %>%
  do(remove_outliers(., "Fe")) %>%
  ungroup()

Fe_parrot <- glmer(Fe ~ site_condition + (1|site_name), data = Fe_woutliers, gaussian(link="log"))
Fe_parrot_v1 <- glmer(Fe ~ site_condition + (1|site_name), data = Fe_woutliers, Gamma(link="log"))
Fe_parrot_v2 <- lmer(Fe ~ site_condition + (1|site_name), data = Fe_woutliers)

AIC(Fe_parrot, Fe_parrot_v2, Fe_parrot_v1)
check_model(Fe_parrot_v2)
ks.test(residuals(Fe_parrot_v2), "pnorm", mean = mean(residuals(Fe_parrot_v2)), sd = sd(residuals(Fe_parrot_v2)))

testDispersion(Fe_parrot_v2)
simulationOutput <- simulateResiduals(fittedModel =Fe_parrot_v2, plot = F)
plot(simulationOutput)

sim_residuals <- simulateResiduals(Fe_parrot_v2)
plot(sim_residuals)
testUniformity(sim_residuals)    
testDispersion(sim_residuals)   
testZeroInflation(sim_residuals) 

shapiro.test(residuals(Fe_parrot_v2))
ad.test(residuals(Fe_parrot_v2)) 

ks.test(residuals(Fe_parrot_v2), "pnorm", mean = mean(residuals(Fe_parrot_v2)), sd = sd(residuals(Fe_parrot_v2)))
leveneTest(Fe ~ site_condition, Fe_woutliers)

# NS
summary(glht(Fe_parrot_v2, linfct = mcp(site_condition = "Tukey")), test = adjusted("fdr"))

Fe_woutliers <- merge_popfish_micronut %>%
  group_by(site_name) %>%
  do(remove_outliers(., "Fe")) %>%
  ungroup()

b <- boxcox(lm(Fe_woutliers$Fe ~ 1))
lambda <- b$x[which.max(b$y)]
lambda
new_x_exact <- (Fe_woutliers$Fe ^ lambda - 1) / lambda
Fe_woutliers$Fe_box <- new_x_exact

# Effect of site
Fe_parrot_site_m1 <- lm(Fe_box ~ site_name, data =Fe_woutliers)
Fe_parrot_site_m1 <- aov(Fe_box ~ site_name, data =Fe_woutliers)
Fe_parrot_site_m1 <- kruskal.test(Fe ~ site_name, data=Fe_woutliers)
check_heteroscedasticity(Fe_parrot_site_m1)
check_model(Fe_parrot_site_m1)

## not relevant
pairwise_diff <- function(data) {
  combs <- combn(data$Fe, 2)  # Get all pairwise combinations of N
  abs_diff <- abs(combs[1, ] - combs[2, ])  # Calculate the absolute difference between the pairs
  return(abs_diff)
}

# Apply the function to each condition (B, LC, LR)
B_diff <- pairwise_diff(subset(Fe_woutliers, site_condition == "B"))
LC_diff <- pairwise_diff(subset(Fe_woutliers, site_condition == "LC"))
LR_diff <- pairwise_diff(subset(Fe_woutliers, site_condition == "LR"))

all_diffs <- data.frame(
  diff = c(B_diff, LC_diff, LR_diff),
  condition = factor(rep(c("B", "LC", "LR"), times = c(length(B_diff), length(LC_diff), length(LR_diff))))
)

ggplot(all_diffs, aes(x = condition, y = diff)) +
  geom_boxplot() +
  labs(title = "Comparison of Pairwise Differences", y = "Pairwise Differences", x = "Condition") +
  theme_minimal()

# Calculate variance for each condition
variance_b <- mean(all_diffs$diff[all_diffs$condition == "B"])
variance_LC <- mean(all_diffs$diff[all_diffs$condition == "LC"])
variance_LR <- mean(all_diffs$diff[all_diffs$condition == "LR"])

variability_table <- tibble(
  condition = c("B", "LC", "LR"),
  mean_variability = c(variance_b, variance_LC, variance_LR)
)

var_diff_table <- variability_table %>%
  rename(Condition1 = condition, MeanVar1 = mean_variability) %>%
  crossing(
    variability_table %>% rename(Condition2 = condition, MeanVar2 = mean_variability)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((MeanVar1 - MeanVar2) / MeanVar2) * 100
  )

var_diff_table_clean <- var_diff_table %>%
  filter(Condition1 < Condition2)

kruskal.test(diff ~ condition, data = all_diffs)
dunn.test(all_diffs$dif, all_diffs$condition, list=T, method="bh")

# =========
# Selenium
# =========
Se_parrot <- glmer(Se ~ site_condition + (1|site_name), data = merge_popfish_micronut, family=gaussian(link="log"))

check_model(Se_parrot)
ad.test(residuals(Se_parrot))
check_heteroscedasticity(Se_parrot)
check_normality(Se_parrot)
shapiro.test(residuals(Se_parrot))
ks.test(residuals(Se_parrot), "pnorm", mean = mean(residuals(Se_parrot)), sd = sd(residuals(Se_parrot)))
testDispersion(Se_parrot)
simulationOutput <- simulateResiduals(fittedModel =Se_parrot, plot = F)
plot(simulationOutput)
summary(Se_parrot)

variability_Se <- merge_popfish_micronut %>%
  group_by(site_condition) %>%
  summarise(
    Mean = mean(Se),
    SD = sd(Se),
    Variance = var(Se),
    CV = (SD / Mean) * 100  
  )

diff_percent_table_Se <- variability_Se %>%
  select(site_condition, Mean) %>%
  rename(Mean1 = Mean, Condition1 = site_condition) %>%
  crossing(
    variability_Se %>% select(site_condition, Mean) %>% 
      rename(Mean2 = Mean, Condition2 = site_condition)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((Mean1 - Mean2) / Mean2) * 100
  )

summary(glht(Se_parrot, linfct = mcp(site_condition = "Tukey")), test = adjusted("holm"))

# Effect of site
Se_parrot_site_m1 <- kruskal.test(Se ~ site_name, data=merge_popfish_micronut) 
dunn.test(merge_popfish_micronut$Se, merge_popfish_micronut$site_name, list=T, method="bh")


pairwise_diff <- function(data) {
  combs <- combn(data$Se, 2)  
  abs_diff <- abs(combs[1, ] - combs[2, ])  
  return(abs_diff)
}

# Apply the function to each condition (B, LC, LR)
B_diff <- pairwise_diff(subset(merge_popfish_micronut, site_condition == "B"))
LC_diff <- pairwise_diff(subset(merge_popfish_micronut, site_condition == "LC"))
LR_diff <- pairwise_diff(subset(merge_popfish_micronut, site_condition == "LR"))

all_diffs <- data.frame(
  diff = c(B_diff, LC_diff, LR_diff),
  condition = factor(rep(c("B", "LC", "LR"), times = c(length(B_diff), length(LC_diff), length(LR_diff))))
)

ggplot(all_diffs, aes(x = condition, y = diff)) +
  geom_boxplot() +
  labs(title = "Comparison of Pairwise Differences", y = "Pairwise Differences", x = "Condition") +
  theme_minimal()

# Calculate variance for each condition
variance_b <- mean(all_diffs$diff[all_diffs$condition == "B"])
variance_LC <- mean(all_diffs$diff[all_diffs$condition == "LC"])
variance_LR <- mean(all_diffs$diff[all_diffs$condition == "LR"])

variability_table <- tibble(
  condition = c("B", "LC", "LR"),
  mean_variability = c(variance_b, variance_LC, variance_LR)
)

var_diff_table <- variability_table %>%
  rename(Condition1 = condition, MeanVar1 = mean_variability) %>%
  crossing(
    variability_table %>% rename(Condition2 = condition, MeanVar2 = mean_variability)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((MeanVar1 - MeanVar2) / MeanVar2) * 100
  )

var_diff_table_clean <- var_diff_table %>%
  filter(Condition1 > Condition2)

kruskal.test(diff ~ condition, data = all_diffs)
dunn.test(all_diffs$dif, all_diffs$condition, list=T, method="bh")

# By site

# Create a list of pairwise diffs for each site_name
all_diffs_list <- merge_popfish_micronut %>%
  split(.$site_name) %>%
  map(~ pairwise_diff(.x))  # Apply your custom function

# Combine into one dataframe
all_diffs <- tibble(
  site_name = rep(names(all_diffs_list), times = map_int(all_diffs_list, length)),
  diff = unlist(all_diffs_list)
)

variability_table <- all_diffs %>%
  group_by(site_name) %>%
  summarise(mean_variability = mean(diff))

var_diff_table <- variability_table %>%
  rename(Site1 = site_name, MeanVar1 = mean_variability) %>%
  crossing(
    variability_table %>% rename(Site2 = site_name, MeanVar2 = mean_variability)
  ) %>%
  filter(Site1 != Site2) %>%
  mutate(
    Percent_Diff = ((MeanVar1 - MeanVar2) / MeanVar2) * 100
  )

var_diff_table_clean <- var_diff_table %>%
  filter(Site1 < Site2)

kruskal.test(diff ~ site_name, data = all_diffs)
dunn.test(all_diffs$diff, all_diffs$site_name, method = "bh", list=T)

# =========
# Magnesium
# =========
Mg_parrot <- glmer(Mg ~ site_condition + (1|site_name), data = merge_popfish_micronut, family=gaussian(link="log"))
Mg_parrot_v2 <- glmer(Mg ~ site_condition + (1|site_name), data = merge_popfish_micronut, family=gaussian(link="identity"))

AIC(Mg_parrot, Mg_parrot_v2)

check_model(Mg_parrot)
ad.test(residuals(Mg_parrot))
check_heteroscedasticity(Mg_parrot)
check_normality(Mg_parrot)
shapiro.test(residuals(Mg_parrot))
ks.test(residuals(Mg_parrot), "pnorm", mean = mean(residuals(Mg_parrot)), sd = sd(residuals(Mg_parrot)))
testDispersion(Mg_parrot)
simulationOutput <- simulateResiduals(fittedModel =Mg_parrot, plot = F)
plot(simulationOutput)
summary(Mg_parrot)

variability_Mg <- merge_popfish_micronut %>%
  group_by(site_condition) %>%
  summarise(
    Mean = mean(Mg),
    SD = sd(Mg),
    Variance = var(Mg),
    CV = (SD / Mean) * 100 
  )

diff_percent_table_Mg <- variability_Mg %>%
  select(site_condition, Mean) %>%
  rename(Mean1 = Mean, Condition1 = site_condition) %>%
  crossing(
    variability_Mg %>% select(site_condition, Mean) %>% 
      rename(Mean2 = Mean, Condition2 = site_condition)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((Mean1 - Mean2) / Mean2) * 100
  )


# Effect of Site
Mg_parrot_site_m1 <- kruskal.test(Mg ~ site_name, data=merge_popfish_micronut) #NS

pairwise_diff <- function(data) {
  combs <- combn(data$Mg, 2)  
  abs_diff <- abs(combs[1, ] - combs[2, ])  
  return(abs_diff)
}

# Apply the function to each condition (B, LC, LR)
B_diff <- pairwise_diff(subset(merge_popfish_micronut, site_condition == "B"))
LC_diff <- pairwise_diff(subset(merge_popfish_micronut, site_condition == "LC"))
LR_diff <- pairwise_diff(subset(merge_popfish_micronut, site_condition == "LR"))

all_diffs <- data.frame(
  diff = c(B_diff, LC_diff, LR_diff),
  condition = factor(rep(c("B", "LC", "LR"), times = c(length(B_diff), length(LC_diff), length(LR_diff))))
)

ggplot(all_diffs, aes(x = condition, y = diff)) +
  geom_boxplot() +
  labs(title = "Comparison of Pairwise Differences", y = "Pairwise Differences", x = "Condition") +
  theme_minimal()

# Calculate variance for each condition
variance_b <- mean(all_diffs$diff[all_diffs$condition == "B"])
variance_LC <- mean(all_diffs$diff[all_diffs$condition == "LC"])
variance_LR <- mean(all_diffs$diff[all_diffs$condition == "LR"])

variability_table <- tibble(
  condition = c("B", "LC", "LR"),
  mean_variability = c(variance_b, variance_LC, variance_LR)
)

var_diff_table <- variability_table %>%
  rename(Condition1 = condition, MeanVar1 = mean_variability) %>%
  crossing(
    variability_table %>% rename(Condition2 = condition, MeanVar2 = mean_variability)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((MeanVar1 - MeanVar2) / MeanVar2) * 100
  )

var_diff_table_clean <- var_diff_table %>%
  filter(Condition1 < Condition2)

kruskal.test(diff ~ condition, data = all_diffs)
dunn.test(all_diffs$dif, all_diffs$condition, list=T, method="bh")


# ===========
# Phosphorus
# ===========
P_parrot <- glmer(P ~ site_condition + (1|site_name), data = merge_popfish_micronut, family=gaussian(link="log"))
P_parrot_v2 <- glmer(P ~ site_condition + (1|site_name), data = merge_popfish_micronut)

check_model(P_parrot)
ks.test(residuals(P_parrot), "pnorm", mean = mean(residuals(P_parrot)), sd = sd(residuals(P_parrot)))
testDispersion(P_parrot)
simulationOutput <- simulateResiduals(fittedModel =P_parrot, plot = F)
plot(simulationOutput)
summary(P_parrot)

summary(P_parrot)

variability_P <- merge_popfish_micronut %>%
  group_by(site_condition) %>%
  summarise(
    Mean = mean(P),
    SD = sd(P),
    Variance = var(P),
    CV = (SD / Mean) * 100  
  )

diff_percent_table_P <- variability_P %>%
  select(site_condition, Mean) %>%
  rename(Mean1 = Mean, Condition1 = site_condition) %>%
  crossing(
    variability_P %>% select(site_condition, Mean) %>% 
      rename(Mean2 = Mean, Condition2 = site_condition)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((Mean1 - Mean2) / Mean2) * 100
  )

# Effect of Site
P_parrot_site_m1 <- kruskal.test(P ~ site_name, data=merge_popfish_micronut)
dunn.test(merge_popfish_micronut$P, merge_popfish_micronut$site_name, list=T, method="bh")

pairwise_diff <- function(data) {
  combs <- combn(data$P, 2)  # Get all pairwise combinations of N
  abs_diff <- abs(combs[1, ] - combs[2, ])  # Calculate the absolute difference between the pairs
  return(abs_diff)
}

# Apply the function to each condition (B, LC, LR)
B_diff <- pairwise_diff(subset(merge_popfish_micronut, site_condition == "B"))
LC_diff <- pairwise_diff(subset(merge_popfish_micronut, site_condition == "LC"))
LR_diff <- pairwise_diff(subset(merge_popfish_micronut, site_condition == "LR"))

all_diffs <- data.frame(
  diff = c(B_diff, LC_diff, LR_diff),
  condition = factor(rep(c("B", "LC", "LR"), times = c(length(B_diff), length(LC_diff), length(LR_diff))))
)

ggplot(all_diffs, aes(x = condition, y = diff)) +
  geom_boxplot() +
  labs(title = "Comparison of Pairwise Differences", y = "Pairwise Differences", x = "Condition") +
  theme_minimal()

# Calculate variance for each condition
variance_b <- mean(all_diffs$diff[all_diffs$condition == "B"])
variance_LC <- mean(all_diffs$diff[all_diffs$condition == "LC"])
variance_LR <- mean(all_diffs$diff[all_diffs$condition == "LR"])

variability_table <- tibble(
  condition = c("B", "LC", "LR"),
  mean_variability = c(variance_b, variance_LC, variance_LR)
)

var_diff_table <- variability_table %>%
  rename(Condition1 = condition, MeanVar1 = mean_variability) %>%
  crossing(
    variability_table %>% rename(Condition2 = condition, MeanVar2 = mean_variability)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((MeanVar1 - MeanVar2) / MeanVar2) * 100
  )

var_diff_table_clean <- var_diff_table %>%
  filter(Condition1 < Condition2)

kruskal.test(diff ~ condition, data = all_diffs)
dunn.test(all_diffs$dif, all_diffs$condition, list=T, method="bh")

# ===========
# Zinc
# ===========
Zn_woutliers <- merge_popfish_micronut %>%
  do(remove_outliers(., "Zn")) %>%
  ungroup()

Zn_parrot_v1 <- glmer(Zn ~ site_condition + (1|site_name), data = Zn_woutliers, family=gaussian(link="log"))
Zn_parrot_v2 <- lmer((Zn) ~ site_condition + (1|site_name), data = Zn_woutliers)

AIC(Zn_parrot_v2,Zn_parrot_v1)

check_model(Zn_parrot_v1)
ks.test(residuals(Zn_parrot_v1), "pnorm", mean = mean(residuals(Zn_parrot_v1)), sd = sd(residuals(Zn_parrot_v1)))
testDispersion(Zn_parrot_v1)
simulationOutput <- simulateResiduals(fittedModel =Zn_parrot_v1, plot = F)
plot(simulationOutput)
sim_residuals <- simulateResiduals(Zn_parrot_v1)
plot(sim_residuals)
testUniformity(sim_residuals)    
testDispersion(sim_residuals)   
testZeroInflation(sim_residuals) 
shapiro.test(residuals(Zn_parrot_v1))
ad.test(residuals(Zn_parrot_v1)) 
ks.test(residuals(Zn_parrot_v1), "pnorm", mean = mean(residuals(Zn_parrot_v1)), sd = sd(residuals(Zn_parrot_v1)))
leveneTest(Zn ~ site_condition, Zn_woutliers)

variability_Zn <- Zn_woutliers %>%
  group_by(site_condition) %>%
  summarise(
    Mean = mean(Zn),
    SD = sd(Zn),
    Variance = var(Zn),
    CV = (SD / Mean) * 100  
  )

diff_percent_table_Zn <- variability_Zn %>%
  select(site_condition, Mean) %>%
  rename(Mean1 = Mean, Condition1 = site_condition) %>%
  crossing(
    variability_Zn %>% select(site_condition, Mean) %>% 
      rename(Mean2 = Mean, Condition2 = site_condition)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((Mean1 - Mean2) / Mean2) * 100
  )

summary(glht(Zn_parrot_v1, linfct = mcp(site_condition = "Tukey")), test = adjusted("holm"))

# Effect of Site
Zn_parrot_site_m1 <- kruskal.test(Zn ~ site_name, data=Zn_woutliers)
dunn.test(merge_popfish_micronut$Zn, merge_popfish_micronut$site_name, list=T, method="bh")

pairwise_diff <- function(data) {
  combs <- combn(data$Zn, 2)  
  abs_diff <- abs(combs[1, ] - combs[2, ])  
  return(abs_diff)
}

# Apply the function to each condition (B, LC, LR)
B_diff <- pairwise_diff(subset(Zn_woutliers, site_condition == "B"))
LC_diff <- pairwise_diff(subset(Zn_woutliers, site_condition == "LC"))
LR_diff <- pairwise_diff(subset(Zn_woutliers, site_condition == "LR"))

all_diffs <- data.frame(
  diff = c(B_diff, LC_diff, LR_diff),
  condition = factor(rep(c("B", "LC", "LR"), times = c(length(B_diff), length(LC_diff), length(LR_diff))))
)

ggplot(all_diffs, aes(x = condition, y = diff)) +
  geom_boxplot() +
  labs(title = "Comparison of Pairwise Differences", y = "Pairwise Differences", x = "Condition") +
  theme_minimal()

# Calculate variance for each condition
variance_b <- mean(all_diffs$diff[all_diffs$condition == "B"])
variance_LC <- mean(all_diffs$diff[all_diffs$condition == "LC"])
variance_LR <- mean(all_diffs$diff[all_diffs$condition == "LR"])

variability_table <- tibble(
  condition = c("B", "LC", "LR"),
  mean_variability = c(variance_b, variance_LC, variance_LR)
)

var_diff_table <- variability_table %>%
  rename(Condition1 = condition, MeanVar1 = mean_variability) %>%
  crossing(
    variability_table %>% rename(Condition2 = condition, MeanVar2 = mean_variability)
  ) %>%
  filter(Condition1 != Condition2) %>%
  mutate(
    Percent_Diff = ((MeanVar1 - MeanVar2) / MeanVar2) * 100
  )

var_diff_table_clean <- var_diff_table %>%
  filter(Condition1 < Condition2)

kruskal.test(diff ~ condition, data = all_diffs)
dunn.test(all_diffs$dif, all_diffs$condition, list=T, method="bh")

## Merge micronut 
merge_popfish_micronut <- merge(meta_popfish, micronut_parrot, by.x="sample",by.y="sample", no.dups=T)
merge_popfish_micronut <- as_tibble(merge_popfish_micronut)
merge_popfish_micronut <- merge_popfish_micronut %>%
  rename(site_condition = site_condition.x)
merge_popfish_micronut$SL <- as.numeric(gsub(",", ".", merge_popfish_micronut$SL))
#write.csv(merge_popfish_micronut, "POPAMA_parrot_micronut.csv")

# Replace site_condition names before computing plots
replacement_rules <- c("I" = "B", "H" = "LR", "P" = "LC")
merge_popfish_micronut$site_condition <- str_replace_all(merge_popfish_micronut$site_condition, replacement_rules)

# Potassium
desired_order <- c("B", "LC", "LR")
merge_popfish_micronut$site_condition <- factor(merge_popfish_micronut$site_condition, levels = desired_order)

K_fish <- ggplot(data=merge_popfish_micronut, aes(x=site_condition, y=K, fill=site_condition)) +
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=0.9, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.7) +  # Set alpha value for transparency
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2") +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("K_fish") +
  xlab("") 

nogrid_K_fish <- K_fish + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

nogrid_K_fish

mean_values <- merge_popfish_micronut %>%
  group_by(site_condition) %>%
  summarise(mean_nutrient = mean(K, na.rm = TRUE))

# Chromium 
desired_order <- c("barrier", "lagoon_city", "lagoon_remote")
Cr_woutliers$site_condition <- factor(Cr_woutliers$site_condition, levels = desired_order)

Cr_fish <- ggplot(data=Cr_woutliers, aes(x=site_condition, y=Cr, fill=site_condition)) +
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=0.9, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.7) +  # Set alpha value for transparency
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2") +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Cr_fish") +
  xlab("") 

nogrid_Cr_fish <- Cr_fish + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

nogrid_Cr_fish

mean_values <- Cr_woutliers %>%
  group_by(site_condition) %>%
  summarise(mean_nutrient = mean(Cr, na.rm = TRUE))

# Iron
desired_order <- c("barrier", "lagoon_city", "lagoon_remote")
Fe_woutliers$site_condition <- factor(Fe_woutliers$site_condition, levels = desired_order)

Fe_fish <- ggplot(data=Fe_woutliers, aes(x=site_name, y=Fe, fill=site_name)) +
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=0.9, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.7) +  # Set alpha value for transparency
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2") +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Fe_fish") +
  xlab("") 

nogrid_Fe_fish <- Fe_fish + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

nogrid_Fe_fish

# Selenium
Se_woutliers$site_condition <- factor(Se_woutliers$site_condition, levels = desired_order)

Se_fish <- ggplot(data=merge_popfish_micronut, aes(x=site_condition, y=Se, fill=site_condition)) +
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=0.9, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.7) +  # Set alpha value for transparency
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2") +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Se_fish") +
  xlab("") 

nogrid_Se_fish <- Se_fish + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

nogrid_Se_fish

mean_values <- merge_popfish_micronut %>%
  group_by(site_condition) %>%
  summarise(mean_nutrient = mean(Se, na.rm = TRUE))

# Magnesium
Mg_woutliers$site_condition <- factor(Mg_woutliers$site_condition, levels = desired_order)

Mg_fish <- ggplot(data=merge_popfish_micronut, aes(x=site_condition, y=Mg, fill=site_condition)) +
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=0.9, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.7) +  # Mgt alpha value for transparency
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2") +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Mg_fish") +
  xlab("") 

nogrid_Mg_fish <- Mg_fish + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

nogrid_Mg_fish

mean_values <- merge_popfish_micronut %>%
  group_by(site_condition) %>%
  summarise(mean_nutrient = mean(Mg, na.rm = TRUE))

# Phosphorus
P_fish <- ggplot(data=merge_popfish_micronut, aes(x=site_condition, y=P, fill=site_condition)) +
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=0.9, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.7) +  # Pt alpha value for transparency
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2") +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("P_fish") +
  xlab("") 

nogrid_P_fish <- P_fish + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

nogrid_P_fish

mean_values <- merge_popfish_micronut %>%
  group_by(site_condition) %>%
  summarise(mean_nutrient = mean(P, na.rm = TRUE))

# Calcium
desired_order <- c("passe_boueni", "passe_chouazil", "baie_boueni","sada","ilot_mbouini","ngouja")

Ca_woutliers$site_name <- factor(Ca_woutliers$site_name, levels = desired_order)

Ca_fish <- ggplot(data=Ca_woutliers, aes(x=site_name, y=Ca, fill=site_name)) +
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=0.9, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.7) +  # Pt alpha value for transparency
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2") +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Ca_fish") +
  xlab("") 

nogrid_Ca_fish <- Ca_fish + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

nogrid_Ca_fish

# Zinc
Zn_woutliers$site_condition <- factor(Zn_woutliers$site_condition, levels = desired_order)

Zn_fish <- ggplot(data=Zn_woutliers, aes(x=site_condition, y=Zn, fill=site_condition)) +
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=0.9, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.7) +  # Pt alpha value for transparency
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2") +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Zn_fish") +
  xlab("") 

nogrid_Zn_fish <- Zn_fish + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

nogrid_Zn_fish

mean_values <- Zn_woutliers %>%
  group_by(site_condition) %>%
  summarise(mean_nutrient = mean(Zn, na.rm = TRUE))

ggarrange(nogrid_K_fish, nogrid_Cr_fish, nogrid_Zn_fish,nogrid_Se_fish, nogrid_Mg_fish, nogrid_P_fish, ncol = 2, 
          nrow = 3)

