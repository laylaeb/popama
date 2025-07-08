# ==============================================================================================
# Association patterns linking nutrients to ASVs identified in fish gut - example with Potassium
# ==============================================================================================

library(tidyverse)
library(dplyr)
library(phyloseq)
library(phyloseqCompanion)
library(moments)
library(lme4)
library(rsq)
library(lmerTest)
library(ggplot2)
library(MuMIn)
library(nortest)
library(car)
library(RColorBrewer)  
library(scales)  


micronut_fish <- as.data.frame(read.csv("##",header=T))
micronut_fish <- micronut_fish %>% select(K,Cr,Fe,Se,Mg,P,Zn,Ca,Mn,sample) # nutrients are given as ug/g

## Median, Lower and upper bound values for the different conditions
IQR_nutrient_wn <- micronut_fish  %>%
  summarize_at(vars(K:Mn),~ median(.x, na.rm = TRUE))

## Transform values of micronut into portion of 100 g 
# transform as wet weight
percentWater <- 77.5 

dminerals <- micronut_fish %>%
  mutate(Cr = (Cr * ((100 - percentWater)/100))/10, #expressed as mg/100g
         K = (K * ((100 - percentWater)/100))/10, #expressed as mg/100g
         Ca = (Ca * ((100 - percentWater)/100))/10, #expressed as mg/100g
         Mg = (Mg * ((100 - percentWater)/100))/10, #expressed as mg/100g
         P = (P * ((100 - percentWater)/100))/10, #expressed as mg/100g
         Se = (Se * ((100 - percentWater)/100))*100, #expressed as ug/100g
         Fe = (Fe * ((100 - percentWater)/100))/10, #expressed as mg/100g
         Zn = (Zn * ((100 - percentWater)/100))/10,
         Mn = (Mn * ((100 - percentWater)/100))/10)#expressed as mg/100g

## Load rarefied data
popfish_rarefied_2024_fish <- readRDS("##")
select_fish <- subset_samples(popfish_rarefied_2024_fish, sample_type %in% c("fish"))
relab_fish <- phyloseq::transform_sample_counts(select_fish, function(x) x/sum(x))

asv_rarefied_full <- otu_table(relab_fish, taxa_are_rows=F)
asv_rarefied_full <- t(asv_rarefied_full[rowSums(asv_rarefied_full) > 0,])

melt(asv_rarefied_full)
tax_fish <- as.data.frame(tax_table(select_fish))
meta_fish <- sample.data.frame(select_fish)

## All micronutrients and metadata
covariate_fish <- subset(meta_fish, select = c("site_condition", "site_name"))
rownames(dminerals) <- dminerals$sample
dminerals$sample <-NULL
dminerals$sample <- rownames(dminerals)
covariate_fish$sample <- rownames(covariate_fish)
merge_cov_micro <- merge(dminerals, covariate_fish, by="sample")

## Median values for the different conditions
IQR_nutrient <- merge_cov_micro %>%
  group_by(site_condition) %>%
  summarize_at(vars(K:Mn),~ median(.x, na.rm = TRUE))

## Load relative abundance of taxa
micronut_asvs_cond <- read.csv("##", header=T)

## For Potassium (K)
data_subset_K <- micronut_asvs_cond  %>%
  select(-c(Cr, Mg, P, Se, Fe, Zn,Ca, Mn))

## Prepare data for the loop
metadata <- data_subset_K[,c(1,2,3,4)]
asvs <- data_subset_K[,-c(1,2,3,4)]

hist(log(metadata$K))
skewness(metadata$K)

## Looping!

# Metadata and ASV data
metadata <- data_subset_K[, c(1,2,3,4)]
asvs <- data_subset_K[, -c(1,2,3,4)]

# Scale ASVs (z-score standardization)
asvs_scaled <- asvs
asvs_scaled[] <- lapply(asvs, function(x) if(sd(x) != 0) as.numeric(scale(x)) else x)

# Output directory for plots
output_dir <- "202503_model_diagnostics_plots_potassium_v2"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Pre-allocate results dataframe
test_results <- data.frame(
  ASV = colnames(asvs),
  p_value_asv = NA_real_,
  estimate_asv = NA_real_,
  p_value_site_condition_LC = NA_real_,
  estimate_site_condition_LC = NA_real_,
  p_value_site_condition_LR = NA_real_,
  estimate_site_condition_LR = NA_real_,
  Shapiro_p_value = NA_real_,
  Anderson_Darling_p_value = NA_real_,
  Levene_p_value = NA_real_,
  R2_marginal = NA_real_,
  R2_conditional = NA_real_,
  sign_asv = NA_real_,
  random_effect_adjusted = NA,
  Gamma_AIC = NA_real_,
  Gaussian_AIC = NA_real_,
  dispersion_parameter = NA_real_,
  stringsAsFactors = FALSE
)

for (a in seq_len(ncol(asvs_scaled))) {
  this_data <- data.frame(
    sample = rownames(asvs_scaled),
    K = metadata$K,
    site_condition = factor(metadata$site_condition),
    site_name = factor(metadata$site_name),
    asv = asvs_scaled[, a]
  )
  
  # skip ASVs without variation
  if (length(unique(this_data$asv)) == 1) {
    warning(paste("ASV", colnames(asvs_scaled)[a], "has no variation, skipped"))
    next
  }
  
  # fit Gaussian mixed model with random factor for sites
  m_gaussian <- tryCatch(
    lmer(K ~ asv + site_condition + (1 | site_name), data = this_data),
    error = function(e) {
      message("Error fitting Gaussian model for ASV: ", colnames(asvs_scaled)[a])
      return(NULL)
    }
  )
  
  if (is.null(m_gaussian)) next
  
  coef_summary <- summary(m_gaussian)$coefficients
  
  # check if 'asv' is in the model coefficients
  if (!"asv" %in% rownames(coef_summary)) next
  
  p_value_asv <- coef_summary["asv", "Pr(>|t|)"]
  estimate_asv <- coef_summary["asv", "Estimate"]
  
  # extract site_condition effects for LC and LR (fringing reefs exposed or not to cities)
  site_levels <- levels(this_data$site_condition)
  
  get_coef <- function(level) {
    nm <- paste0("site_condition", level)
    if (nm %in% rownames(coef_summary)) coef_summary[nm, ] else c(Estimate=NA, `Pr(>|t|)`=NA)
  }
  
  LC_coef <- get_coef(site_levels[2])
  LR_coef <- get_coef(site_levels[3])
  
  # stats
  residuals_gaussian <- residuals(m_gaussian)
  fitted_gaussian <- fitted(m_gaussian)
  
  shapiro_p <- tryCatch(shapiro.test(residuals_gaussian)$p.value, error=function(e) NA)
  ad_p <- tryCatch(ad.test(residuals_gaussian)$p.value, error=function(e) NA)
  levene_p <- tryCatch(car::leveneTest(residuals_gaussian ~ this_data$site_condition)$`Pr(>F)`[1], error=function(e) NA)
  
  r2_vals <- suppressWarnings(r.squaredGLMM(m_gaussian))
  r2_marginal <- if (!is.null(r2_vals)) r2_vals[1] else NA
  r2_conditional <- if (!is.null(r2_vals)) r2_vals[2] else NA
  
  sign_asv <- sign(estimate_asv)
  dispersion_gaussian <- var(residuals_gaussian)
  
  # try to fit Gamma model (to compare)
  m_gamma <- tryCatch(
    glmer(K ~ asv + site_condition + (1 | site_name), data = this_data, family = Gamma(link = "log")),
    error = function(e) {
      message("Error fitting Gamma model for ASV: ", colnames(asvs_scaled)[a])
      return(NULL)
    }
  )
  
  AIC_gamma <- if (!is.null(m_gamma)) AIC(m_gamma) else NA
  AIC_gaussian <- AIC(m_gaussian)
  
  random_effect_adj <- isSingular(m_gaussian)
  if (random_effect_adj) message("Singularity detected for ASV: ", colnames(asvs_scaled)[a])
  
  # store all results
  test_results[a, ] <- list(
    ASV = colnames(asvs_scaled)[a],
    p_value_asv = p_value_asv,
    estimate_asv = estimate_asv,
    p_value_site_condition_LC = LC_coef["Pr(>|t|)"],
    estimate_site_condition_LC = LC_coef["Estimate"],
    p_value_site_condition_LR = LR_coef["Pr(>|t|)"],
    estimate_site_condition_LR = LR_coef["Estimate"],
    Shapiro_p_value = shapiro_p,
    Anderson_Darling_p_value = ad_p,
    Levene_p_value = levene_p,
    R2_marginal = r2_marginal,
    R2_conditional = r2_conditional,
    sign_asv = sign_asv,
    random_effect_adjusted = random_effect_adj,
    Gamma_AIC = AIC_gamma,
    Gaussian_AIC = AIC_gaussian,
    dispersion_parameter = dispersion_gaussian
  )
  
  # save everything 
  jpeg(filename = file.path(output_dir, paste0("model_diagnostics_", colnames(asvs_scaled)[a], ".jpeg")),
       width = 1200, height = 400)
  par(mfrow = c(1, 3), mar = c(4,4,2,1))
  
  plot(fitted_gaussian, residuals_gaussian,
       xlab = "Fitted values", ylab = "Residuals",
       main = "Residuals vs Fitted (Gaussian Model)")
  abline(h = 0, col = "red", lty = 2)
  
  hist(residuals_gaussian, breaks = 20, col = "lightblue",
       main = "Histogram of Residuals (Gaussian)",
       xlab = "Residuals")
  
  qqnorm(residuals_gaussian, main = "Q-Q Plot of Residuals (Gaussian Model)")
  qqline(residuals_gaussian, col = "red")
  
  dev.off()
}

# resuts as CSV
write.csv(test_results, file = "202503_Potassium_model_test_results_v2.csv", row.names = FALSE)
test_results <- read.csv("202503_Potassium_model_test_results_v2.csv")

res_K_subset <- test_results[rev(order(test_results$R2_marginal)),]
res_K_subset <- res_K_subset[order(res_K_subset$p_value_asv),]

## first derivative
plot(-1*diff(rev(sort(res_K_subset$R2_marginal)))[1:100], type='b')

## second derivative
plot(-1*diff(diff(rev(sort(res_K_subset$R2_marginal))))[1:100], type='b')
abline(h=0)

res_K_subset_table <- res_K_subset[1:60,]
top_20 <- res_K_subset_table[1:60, "ASV"]

## final taxa that are associated with Potassium content
tax_res_K <- merge(res_K_subset_table, tax_fish, by.x="ASV",by.y="row.names")

## Plot
selected_asvs <- top_20[13:24] 

r2_values <- data.frame(asv_name = selected_asvs, r2_marginal = numeric(length(selected_asvs)))  

# produce loop for the selected ASVs
for (asv in selected_asvs) {
  
  this_data <- data.frame(
    K = metadata$K,
    site_condition = metadata$site_condition,
    site_name = metadata$site_name,
    asv_value = asvs[[asv]]
  )
  
  # all good?
  if (all(is.na(this_data$asv_value)) || var(this_data$asv_value, na.rm = TRUE) == 0) {
    message("ASV ", asv, " has missing or constant data. Skipping.")
    next
  }
  
  # fit the chosen model
  model_test <- tryCatch({
    lmer(K ~ asv_value + site_condition + (1 | site_name), data = this_data)
  }, error = function(e) {
    message("Error for ASV: ", asv)
    return(NULL)
  })
  
  if (is.null(model_test)) next
  
  # calculate R2
  r_squared <- r.squaredGLMM(model_test)
  print(paste("R² marginal for ASV ", asv, ": ", round(r_squared[1], 3), sep = ""))
  
  # store R2
  r2_values[r2_values$asv_name == asv, "r2_marginal"] <- r_squared[1]
  
  # create prediction for the model
  prediction_grid <- expand.grid(
    asv_value = seq(min(this_data$asv_value, na.rm = TRUE), max(this_data$asv_value, na.rm = TRUE), length.out = 100),
    site_condition = unique(this_data$site_condition)
  )
  
  # we only focus on the fixed effects for the prediction
  pred_values <- predict(model_test, newdata = prediction_grid, re.form = ~0, se.fit = TRUE)
  
  # keep predictions
  prediction_grid$predicted <- pred_values$fit
  
  # 95% confidence intervals
  conf.low <- prediction_grid$predicted - 1.96 * pred_values$se.fit
  conf.high <- prediction_grid$predicted + 1.96 * pred_values$se.fit
  prediction_grid$conf.low <- conf.low
  prediction_grid$conf.high <- conf.high
  prediction_grid$asv_name <- asv
  
  all_preds <- rbind(all_preds, prediction_grid)
  
  # store real points for plotting
  this_data$asv_name <- asv
  all_points <- rbind(all_points, this_data)
}

# plot relationships of relative abundance of each ASV and micronutrient
ggplot(all_preds, aes(x = asv_value, y = predicted, color = site_condition)) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = site_condition), 
              alpha = 0.2, color = NA) +
  geom_point(data = all_points,
             aes(x = asv_value, y = K, color = site_condition),
             alpha = 0.5, inherit.aes = FALSE) +
  facet_wrap(~asv_name, scales = "free", ncol = 3) +
  labs(
    x = "ASV relative abundance",
    y = "K (mg/100g)",
    color = "Reef Type",
    fill = "Reef Type"
  ) +
  scale_color_manual(values = dark2_colors[1:length(site_condition_labels)], 
                     labels = site_condition_labels) +  
  scale_fill_manual(values = dark2_colors[1:length(site_condition_labels)], 
                    labels = site_condition_labels) +   
  scale_x_continuous(labels = function(x) {
    ifelse(x == 0, "0", label_number()(x))  
  }) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11)
  ) +
  # add R2 on each relationship
  geom_text(data = r2_values, aes(x = Inf, y = Inf, label = paste("R² = ", round(r2_marginal, 3))),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.1, size = 3, color = "black")


## Select ASVs

# define ASVs to be included in the graph
selected_asvs <- c("ASV_00047", "ASV_00279", "ASV_00130", "ASV_00997", 
                   "ASV_00758", "ASV_01247", "ASV_00809", "ASV_00110", 
                   "ASV_00177", "ASV_00154", "ASV_00183", "ASV_00175")

# initialise
all_preds <- data.frame()
all_points <- data.frame()
r2_values <- data.frame(asv_name = selected_asvs, r2_marginal = numeric(length(selected_asvs)))  # Placeholder for R² values

# looping on the chosen ASVs
for (asv in selected_asvs) {
  
  this_data <- data.frame(
    K = metadata$K,
    site_condition = metadata$site_condition,
    site_name = metadata$site_name,
    asv_value = asvs[[asv]]
  )
  
  # sanity check
  if (all(is.na(this_data$asv_value)) || var(this_data$asv_value, na.rm = TRUE) == 0) {
    message("ASV ", asv, " a des données manquantes ou constantes. Ignoré.")
    next
  }
  
  
  model_test <- tryCatch({
    lmer(K ~ asv_value + site_condition + (1 | site_name), data = this_data)
  }, error = function(e) {
    message("Erreur pour ASV: ", asv)
    return(NULL)
  })
  
  if (is.null(model_test)) next
  
  # R2 
  r_squared <- r.squaredGLMM(model_test)
  print(paste("R² marginal pour ASV ", asv, ": ", round(r_squared[1], 3), sep = ""))
  
  
  r2_values[r2_values$asv_name == asv, "r2_marginal"] <- r_squared[1]
  
  # prediction grid
  prediction_grid <- expand.grid(
    asv_value = seq(min(this_data$asv_value, na.rm = TRUE), max(this_data$asv_value, na.rm = TRUE), length.out = 100),
    site_condition = unique(this_data$site_condition)
  )
  
  # predictions with fixed effects only
  pred_values <- predict(model_test, newdata = prediction_grid, re.form = ~0, se.fit = TRUE)
  
  
  prediction_grid$predicted <- pred_values$fit
  
  # confidence intervals
  conf.low <- prediction_grid$predicted - 1.96 * pred_values$se.fit
  conf.high <- prediction_grid$predicted + 1.96 * pred_values$se.fit
  
  prediction_grid$conf.low <- conf.low
  prediction_grid$conf.high <- conf.high
  prediction_grid$asv_name <- asv
  all_preds <- rbind(all_preds, prediction_grid)
  
  
  this_data$asv_name <- asv
  all_points <- rbind(all_points, this_data)
}


site_condition_labels <- c("B" = "Barrier", "LC" = "Fringing City", "LR" = "Fringing Remote")

dark2_colors <- brewer.pal(8, "Dark2")

ggplot(all_preds, aes(x = asv_value, y = predicted, color = site_condition)) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = site_condition), 
              alpha = 0.2, color = NA) +
  geom_point(data = all_points,
             aes(x = asv_value, y = K, color = site_condition),
             alpha = 0.5, inherit.aes = FALSE) +
  facet_wrap(~asv_name, scales = "free", ncol = 3) +
  labs(
    x = "ASV relative abundance",
    y = "K (mg/100g)",
    color = "Reef Type",
    fill = "Reef Type"
  ) +
  scale_color_manual(values = dark2_colors[1:length(site_condition_labels)], 
                     labels = site_condition_labels) +  # Légende des conditions de site
  scale_fill_manual(values = dark2_colors[1:length(site_condition_labels)], 
                    labels = site_condition_labels) +   # Remplissage avec la palette et labels
  scale_x_continuous(labels = function(x) {
    ifelse(x == 0, "0", label_number()(x))  # Appliquer un format personnalisé pour les labels de l'axe x
  }) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11)
  ) +
  # Ajouter le R² marginal sur chaque facette
  geom_text(data = r2_values, aes(x = Inf, y = Inf, label = paste("R² = ", round(r2_marginal, 3))),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.1, size = 3, color = "black")  # Positionner R² sur chaque facette

