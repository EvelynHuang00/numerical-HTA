# Alternative PCA trait/region experiments from HTA_PCA.Rmd.
# These analyses were tested during development but were not retained as the
# primary PCA result in the final report.
# Run from the repository root.

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(recipes)
library(scales)

source("R/hta_core.R")
source("R/pca_helpers.R")

# ----------------------------------------------------------------------------
# 1. Load and clean the numerically dominated cohort
# ----------------------------------------------------------------------------
health_data <- read.csv("data/dirty_v3_path.csv", header = TRUE)
health_data <- na.omit(health_data)
health_data <- health_data %>% select(-c(random_notes, noise_col))

health_homo <- make_pseudo_homogeneous_numeric(
  health_data,
  n = 5000,
  seed = 123
)

# ----------------------------------------------------------------------------
# 2. Early experiment: PCA-defined regions and PCA-defined traits
#    Regions: PC1 x PC2 quantile bins
#    Traits: PC3 x PC4 x PC5 quantile bins
# ----------------------------------------------------------------------------
set.seed(123)
rec <- recipe(~ ., data = health_data) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_predictors())

rec_prep <- prep(rec, training = health_data)
X <- bake(rec_prep, new_data = NULL)
pca_fit_all <- prcomp(X, center = FALSE, scale. = FALSE)
scores_all <- as.data.frame(pca_fit_all$x)

summary(pca_fit_all)

lambda <- pca_fit_all$sdev^2
pve <- lambda / sum(lambda)
cum_pve <- cumsum(pve)

dfp <- data.frame(
  component = seq_along(lambda),
  lambda = lambda,
  cum_pct = 100 * cum_pve
)

scale_factor <- max(dfp$lambda) / 100

ggplot(dfp, aes(x = component)) +
  geom_col(aes(y = lambda), fill = "skyblue", alpha = 0.8) +
  geom_line(aes(y = cum_pct * scale_factor), color = "darkgreen", linewidth = 1) +
  geom_point(aes(y = cum_pct * scale_factor), color = "darkgreen", size = 2) +
  geom_hline(yintercept = 80 * scale_factor, linetype = "dashed", color = "darkgreen") +
  scale_x_continuous(breaks = seq_along(lambda)) +
  scale_y_continuous(
    name = "Eigenvalue (lambda)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Cumulative Var (%)")
  ) +
  labs(x = "Component", title = "Scree plot with cumulative variance") +
  theme_minimal()

bin_q <- function(v, bins) {
  cut(
    v,
    breaks = quantile(
      v,
      probs = seq(0, 1, length.out = bins + 1),
      na.rm = TRUE
    ),
    include.lowest = TRUE,
    labels = FALSE
  )
}

pc1_bin2 <- bin_q(scores_all$PC1, bins = 2)
pc2_bin2 <- bin_q(scores_all$PC2, bins = 2)
pc3_bin3 <- bin_q(scores_all$PC3, bins = 3)
pc4_bin3 <- bin_q(scores_all$PC4, bins = 3)
pc5_bin3 <- bin_q(scores_all$PC5, bins = 3)

health_data$region <- paste0("R_", pc1_bin2, "_", pc2_bin2)
health_data$trait <- paste0("T_", pc3_bin3, "_", pc4_bin3, "_", pc5_bin3)

table(health_data$region)
length(unique(health_data$trait))
head(sort(table(health_data$trait), decreasing = TRUE), 10)

res_pca_tr1 <- compute_hta_pca(
  df = health_data,
  pca_fit = pca_fit_all,
  region_pcs = c(1, 2),
  region_k = 2,
  trait_pcs = c(3, 4, 5),
  trait_k = 3,
  min_cell_n = 5
)

res_pca_tr2 <- compute_hta_pca(
  df = health_data,
  pca_fit = pca_fit_all,
  region_pcs = c(1, 2),
  region_k = 3,
  trait_pcs = c(3, 4, 5),
  trait_k = 3,
  min_cell_n = 5
)

res_pca_tr1$HTA
res_pca_tr1$local_hti
res_pca_tr1$cell_n
res_pca_tr2$HTA
res_pca_tr2$local_hti
res_pca_tr2$cell_n

# Within-cell trait composition for the 2x2 PC1 x PC2 region design.
pc1_br <- qbreaks(scores_all$PC1, 2)
pc2_br <- qbreaks(scores_all$PC2, 2)

df_vis <- health_data %>%
  mutate(
    PC1Bin = cut(scores_all$PC1, breaks = pc1_br, include.lowest = TRUE, right = TRUE),
    PC2Bin = cut(scores_all$PC2, breaks = pc2_br, include.lowest = TRUE, right = TRUE),
    cell = interaction(PC1Bin, PC2Bin, drop = TRUE)
  )

comp <- df_vis %>%
  count(PC1Bin, PC2Bin, cell, trait, name = "n") %>%
  group_by(cell) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(comp, aes(x = "cell", y = prop, fill = trait)) +
  geom_col(position = "stack") +
  coord_flip() +
  facet_grid(PC2Bin ~ PC1Bin) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Trait composition within each PC1×PC2 cell (2×2 grid)",
    x = NULL,
    y = "Proportion within cell"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(angle = 0)
  )

# Same early PCA-defined region/trait design on the pseudo-homogeneous cohort.
rec_homo <- recipe(~ ., data = health_homo) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_predictors())

rec_homo_prep <- prep(rec_homo, training = health_homo)
X_homo_all <- bake(rec_homo_prep, new_data = NULL)
pca_homo_all <- prcomp(X_homo_all, center = FALSE, scale. = FALSE)
scores_h_all <- as.data.frame(pca_homo_all$x)

pc1_bin2 <- bin_q(scores_h_all$PC1, 2)
pc2_bin2 <- bin_q(scores_h_all$PC2, 2)
pc3_bin3 <- bin_q(scores_h_all$PC3, 3)
pc4_bin3 <- bin_q(scores_h_all$PC4, 3)
pc5_bin3 <- bin_q(scores_h_all$PC5, 3)

health_homo$region <- paste0("R_", pc1_bin2, "_", pc2_bin2)
health_homo$trait <- paste0("T_", pc3_bin3, "_", pc4_bin3, "_", pc5_bin3)

res_homo_pc_regions <- compute_hta_pca(
  df = health_homo,
  pca_fit = pca_homo_all,
  region_pcs = c(1, 2),
  region_k = 2,
  trait_pcs = c(3, 4, 5),
  trait_k = 3,
  min_cell_n = 5
)

res_homo_pc_regions$HTA
res_homo_pc_regions$local_hti

tab <- with(health_homo, table(region, trait))
tapply(health_homo$trait, health_homo$region, function(x) length(unique(x)))
round(tab / rowSums(tab), 3)[1:4, 1:10]

prop_mat <- prop.table(tab, 1)
rowSums(prop_mat)
rowSums(round(prop_mat, 3))

# ----------------------------------------------------------------------------
# 3. Alternative Age x BMI region design with PC1-PC3 extreme-bin traits
# ----------------------------------------------------------------------------
health_data <- make_age_bmi_2x2(health_data)
health_homo <- make_age_bmi_2x2(health_homo)

pca_fit <- fit_numeric_pca_excluding_age_bmi(health_data)
pca_fit_homo <- fit_numeric_pca_excluding_age_bmi(health_homo)

scores <- as.data.frame(pca_fit$x)
scores_homo <- as.data.frame(pca_fit_homo$x)

# Original extreme-bin definition: bottom 20%, middle 60%, top 20%.
bin_extreme <- function(v, p = 0.2) {
  ql <- quantile(v, p, na.rm = TRUE, type = 7)
  qh <- quantile(v, 1 - p, na.rm = TRUE, type = 7)
  out <- rep("M", length(v))
  out[v <= ql] <- "L"
  out[v >= qh] <- "H"
  out
}

health_data$trait_extreme <- factor(
  paste0(
    "PC1_", bin_extreme(scores$PC1),
    "__PC2_", bin_extreme(scores$PC2),
    "__PC3_", bin_extreme(scores$PC3)
  )
)

table(health_data$trait_extreme)

res_extreme <- compute_ht_fixed_levels(
  health_data,
  region_col = "region_2x2",
  trait_col = "trait_extreme"
)

res_extreme$HTA
res_extreme$local_hti

health_homo$trait_extreme <- factor(
  paste0(
    "PC1_", bin_extreme(scores_homo$PC1),
    "__PC2_", bin_extreme(scores_homo$PC2),
    "__PC3_", bin_extreme(scores_homo$PC3)
  )
)

table(health_homo$trait_extreme)

res_extreme_homo <- compute_ht_fixed_levels(
  health_homo,
  region_col = "region_2x2",
  trait_col = "trait_extreme"
)

res_extreme_homo$HTA
res_extreme_homo$local_hti
prop.table(table(health_homo$region_2x2, health_homo$trait_extreme), margin = 1)

# ----------------------------------------------------------------------------
# 4. Experimental HTA_ML adjustment using within-region PC1 variance
# ----------------------------------------------------------------------------
# The PC1 tertile trait is recreated here because the original notebook used
# this representation before applying the variance adjustment.
health_data$trait_pc1_tertile <- factor(
  cut(
    scores$PC1,
    breaks = quantile(
      scores$PC1,
      probs = c(0, 1 / 3, 2 / 3, 1),
      na.rm = TRUE,
      type = 7
    ),
    include.lowest = TRUE,
    labels = c("Low", "Mid", "High")
  ),
  levels = c("Low", "Mid", "High")
)

health_homo$trait_pc1_tertile <- factor(
  cut(
    scores_homo$PC1,
    breaks = quantile(
      scores_homo$PC1,
      probs = c(0, 1 / 3, 2 / 3, 1),
      na.rm = TRUE,
      type = 7
    ),
    include.lowest = TRUE,
    labels = c("Low", "Mid", "High")
  ),
  levels = c("Low", "Mid", "High")
)

compute_hta_ml <- function(df, region_col, trait_col, pc1_scores) {
  dat <- df
  dat$pc1_score_ml <- pc1_scores

  ht_res <- compute_ht_fixed_levels(
    dat,
    region_col = region_col,
    trait_col = trait_col
  )

  local_hti_df <- data.frame(
    region = names(ht_res$local_hti),
    local_hti = as.numeric(ht_res$local_hti),
    stringsAsFactors = FALSE
  )

  region_stats <- do.call(
    rbind,
    lapply(split(dat, dat[[region_col]]), function(x) {
      n_r <- sum(!is.na(x$pc1_score_ml))
      sigma2_r <- if (n_r >= 2) var(x$pc1_score_ml, na.rm = TRUE) else 0

      data.frame(
        region = as.character(x[[region_col]][1]),
        n = n_r,
        sigma2 = sigma2_r,
        stringsAsFactors = FALSE
      )
    })
  )

  out <- merge(local_hti_df, region_stats, by = "region", all.x = TRUE)
  out$hti_ml <- out$local_hti * (out$sigma2 / out$n)
  hta_ml <- mean(out$hti_ml, na.rm = TRUE)

  list(
    HTA = ht_res$HTA,
    local_hti = ht_res$local_hti,
    region_stats = out,
    HTA_ML = hta_ml
  )
}

res_pc1_tertile_ml <- compute_hta_ml(
  df = health_data,
  region_col = "region_2x2",
  trait_col = "trait_pc1_tertile",
  pc1_scores = scores$PC1
)

res_pc1_tertile_ml$HTA
res_pc1_tertile_ml$HTA_ML
res_pc1_tertile_ml$region_stats
res_pc1_tertile_ml$region_stats[, c("region", "local_hti", "sigma2", "n", "hti_ml")]

res_pc1_tertile_homo_ml <- compute_hta_ml(
  df = health_homo,
  region_col = "region_2x2",
  trait_col = "trait_pc1_tertile",
  pc1_scores = scores_homo$PC1
)

res_pc1_tertile_homo_ml$HTA
res_pc1_tertile_homo_ml$HTA_ML
res_pc1_tertile_homo_ml$region_stats
res_pc1_tertile_homo_ml$region_stats[, c("region", "local_hti", "sigma2", "n", "hti_ml")]

# ----------------------------------------------------------------------------
# 5. Early single-configuration k-means trait prototype
#    The later, fuller k/PC-set sweep is in 03_clustering_sensitivity.R.
# ----------------------------------------------------------------------------
set.seed(1)
pcs_use <- c("PC1", "PC2", "PC3", "PC4")
K <- 6

km_het <- kmeans(scores[, pcs_use], centers = K, nstart = 50)
health_data$trait_cluster <- factor(paste0("Cluster_", km_het$cluster))
table(health_data$trait_cluster)

res_clust_het <- compute_ht_fixed_levels(
  health_data,
  region_col = "region_2x2",
  trait_col = "trait_cluster"
)

res_clust_het$HTA
res_clust_het$local_hti

set.seed(1)
km_homo <- kmeans(scores_homo[, pcs_use], centers = K, nstart = 50)
health_homo$trait_cluster <- factor(paste0("Cluster_", km_homo$cluster))
table(health_homo$trait_cluster)

res_clust_homo <- compute_ht_fixed_levels(
  health_homo,
  region_col = "region_2x2",
  trait_col = "trait_cluster"
)

res_clust_homo$HTA
res_clust_homo$local_hti
