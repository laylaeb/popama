# =================================================
# Alpha and Beta Diversity Analyses for parrotfish
# =================================================

## Alpha diversity
fish_only_beta <- subset_samples(popfish_rarefied_2024_fish, sample_type %in% c("fish"))
meta_fish_beta <- sample.data.frame(fish_only_beta)
sort(sample_sums(fish_only_beta))
replacement_rules <- c("I" = "B", "H" = "LR", "P" = "LC")
meta_fish_beta$site_condition <- str_replace_all(meta_fish_beta$site_condition, replacement_rules)

diversity_popfish <- microbiome::alpha(
  popfish_rarefied_2024_fish,
  index = c("observed", 
            "diversity_shannon", "evenness_pielou","diversity_gini_simpson"
  )
)
diversity_popfish_indices <- phyloseq::sample_data(diversity_popfish)

# Alpha div - Observed Richness
fish_only_alpha <- phyloseq::merge_phyloseq(fish_only_beta, diversity_popfish_indices)
meta_diversity_fish <- sample.data.frame(fish_only_alpha)

hist(meta_diversity_fish$observed) #skewed to the left

remove_outliers <- function(df, column) {
  Q1 <- quantile(df[[column]], 0.25)
  Q3 <- quantile(df[[column]], 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  df %>% filter(df[[column]] >= lower_bound & df[[column]] <= upper_bound)
}

meta_diversity_fish_no_outliers <- meta_diversity_fish %>%
  group_by(site_name) %>%
  do(remove_outliers(., "observed")) %>%
  ungroup()

# Plot observed richness by site
desired_order <- c("passe_boueni", "passe_chouazil", "baie_boueni", "sada", "ilot_mbouini", "ngouja")
meta_diversity_fish_no_outliers$site_name<- factor(meta_diversity_fish_no_outliers$site_name, levels = desired_order)

OR_fish <- ggplot(data=meta_diversity_fish_no_outliers, aes(x=site_name, y=observed, fill=site_condition)) +
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=0.9, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.7) +  # Set alpha value for transparency
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2") +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("OR_fish") +
  xlab("") 

nogrid_OR_fish <- OR_fish + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

group_means <- meta_diversity_fish_no_outliers %>%
  group_by(site_name) %>%
  summarize(mean_value = mean(observed, na.rm = TRUE), std=sd(observed, na.rm=T))

meta_diversity_fish_no_outliers$site_name <- as.factor(meta_diversity_fish_no_outliers$site_name)

# Compute stats for site
kruskal_observed <- kruskal.test(meta_diversity_fish_no_outliers$observed ~ site_name, data=meta_diversity_fish_no_outliers)
OR_site <- dunn.test(meta_diversity_fish_no_outliers$observed, meta_diversity_fish_no_outliers$site_name, list=T, method="bh")

## Beta Diversity

merged_popfish_unrarefied <- readRDS("202408_merged_popfish_filtered_final_fish.RData")
merged_popfish_unrarefied <- subset_samples(merged_popfish_unrarefied, sample != "G12") # outlier from PB

# Aitchison distance 
meta_fish_ait <- sample.data.frame(merged_popfish_unrarefied)

physeq_clr <- microbiome::transform(merged_popfish_unrarefied, "clr")

fish_clr <- subset_samples(physeq_clr, sample_type %in% c("fish"))
asv_fish_clr <- otu_table(fish_clr, taxa_are_rows=T)

fish_aitchison <- phyloseq::distance(fish_clr, method = "euclidean")

aitchison_pcoa <- ecodist::pco(fish_aitchison)

aitchison_fish_poca_df<- data.frame(pcoa1 = aitchison_pcoa$vectors[,1],
                                   pcoa2 = aitchison_pcoa$vectors[,2],
                                   pcoa3 = aitchison_pcoa$vectors[,3] )

aitchison_fish_poca_df <- cbind(aitchison_fish_poca_df,
                               meta_fish = meta_fish_ait$site_name)

pcoa_aitch_fish_plot <- ggplot(data = aitchison_fish_poca_df, aes(x=pcoa1, y=pcoa3,color=meta_fish)) +
  geom_point() +
  labs(x = "PC1",
       y = "PC3",
       title = "Atichison PCoA") +
  theme(title = element_text(size = 10)) # makes titles smaller

pcoa_aitch_fish_plot

ellipse_aitchison <- ggplot(data = aitchison_fish_poca_df , aes(x = pcoa1, y = pcoa2, color = meta_fish, fill = meta_fish)) +
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  stat_ellipse(aes(group = meta_fish),
               color = "black",
               alpha = 0.2,
               level = 0.95,
               type = "norm",
               geom = "polygon") +
  geom_point(size = 4) +  # Move this layer after stat_ellipse
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

variance_proportion_ai <-aitchison_pcoa$values / sum(aitchison_pcoa$values)

meta_fish_beta <- meta_fish_beta[]

# Betadisper
bdisp_fish_a<- betadisper(fish_aitchison, meta_fish_ait$site_name, type=c("centroid"))
bdisp_fish_a
aov_fish_a <-anova(bdisp_fish_a)
permutest(bdisp_fish_a, pairwise=T, permutation=how(nperm=999))

bdisp_fish_df_a <- as.data.frame(bdisp_fish_a$distances)
names(bdisp_fish_df_a) <- "distances"  

merged_bdisp_fish_a <- merge(meta_fish_beta, bdisp_fish_df_a, by = "row.names")

# Réordonnez les niveaux de site_name selon l'ordre souhaité
merged_bdisp_fish_a$site_name <- factor(merged_bdisp_fish_a$site_name, 
                                        levels = c("passe_boueni", "passe_chouazil", "baie_boueni", 
                                                   "sada", "ilot_mbouini", "ngouja"))

betadisp_boxplot_a <- ggplot(merged_bdisp_fish_a, aes(x = site_name, y = distances)) +
  geom_violin(aes(fill = site_condition), color = "black", alpha = 0.5) +  
  geom_boxplot(aes(fill = site_condition), color = "black", width = 0.1) + 
  geom_jitter(position=position_jitter(0.2), alpha=0.7) +  # Set alpha value for transparency
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Aitchison",
       x = "Site Name",
       y = "Distance",
       color = "Site Condition") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))


# PERMANOVA Aitchison
fish_ado_ai_cond <- adonis2(fish_aitchison ~ site_condition, data = meta_fish_ait, permutations = 999, strata = meta_fish_ait$site_name) ## NS
fish_ado_ai <- adonis2(fish_aitchison ~ site_name, permutations = 999, data=meta_fish_ait)
pairwise_fish_aitchison <- pairwise.adonis2(fish_aitchison ~ site_name, p.adjust.m="fdr", data=meta_fish_ait)

# In case when PERMDISP results are significant we have to complement our results by providing the results from manylgm models
#ab_fish <- mvabund(t(asv_fish_nmds_f))
#model_fish_site <- manyglm(ab_fish ~ site_name,
#                  data = meta_fish_beta, family="negative.binomial")

#anova_fish_site <- anova(model_fish_site, p.uni= "adjusted", nBoot = 99, pairwise.comp=meta_fish_beta$site_name, show.time=T)
#anova_fish_site <- anova_fish_site$uni.p


## Subset table
# Filtre pour les ASVs (0.001% et 5 individus)
asv_RA <- transform_sample_counts(fish_only_beta, function(x) x / sum(x))

asv_level <- t(otu_table(asv_RA, taxa_are_rows=T))
tax_level_asv <- tax_table(asv_RA)

seuil_presence <- 5  

# Calcul de l'abondance moyenne par ASV
mean_abundance <- colMeans(asv_level)

# Nombre d'échantillons où chaque ASV est présent (abondance > 0)
present_in_n_samples <- colSums(asv_level > 0)

keep_asvs <- names(mean_abundance[mean_abundance > 0.00001 & present_in_n_samples >= seuil_presence])

asv_filtered <- asv_level[, keep_asvs]
write.csv(asv_filtered, "202506_filtered_asv_0.001_5ind.csv")
