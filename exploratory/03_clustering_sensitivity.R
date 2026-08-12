# K-means trait-definition experiments from HTA_clustering.Rmd.
# These are preserved as exploratory sensitivity analyses and are not part of
# the main result set presented in the final report.
# Run from the repository root.

library(tidyverse)
library(dplyr)
library(ggplot2)
library(cluster)

source("R/hta_core.R")
source("R/pca_helpers.R")
source("R/clustering_helpers.R")

# ----------------------------------------------------------------------------
# 1. Prepare heterogeneous and pseudo-homogeneous datasets
# ----------------------------------------------------------------------------
health_data <- read.csv("data/dirty_v3_path.csv", header = TRUE)
health_data <- na.omit(health_data)
health_data <- health_data %>% select(-c(random_notes, noise_col))

health_homo <- make_pseudo_homogeneous_numeric(
  health_data,
  n = 5000,
  seed = 123
)

health_data <- make_age_bmi_2x2(health_data)
health_homo <- make_age_bmi_2x2(health_homo)

pca_fit <- fit_numeric_pca_excluding_age_bmi(health_data)
pca_fit_homo <- fit_numeric_pca_excluding_age_bmi(health_homo)

scores <- as.data.frame(pca_fit$x)
scores_homo <- as.data.frame(pca_fit_homo$x)

stopifnot(nrow(health_data) == nrow(scores))
stopifnot(nrow(health_homo) == nrow(scores_homo))

# ----------------------------------------------------------------------------
# 2. Sweep number of clusters across PC sets
# ----------------------------------------------------------------------------
pc_sets <- list(
  PC4 = paste0("PC", 1:4),
  PC6 = paste0("PC", 1:6),
  PC8 = paste0("PC", 1:8)
)

k_values <- 2:15

het_results <- sweep_pcsets_k_ht(
  df = health_data,
  scores = scores,
  pc_sets = pc_sets,
  k_values = k_values,
  region_col = "region_2x2",
  nstart = 100,
  max_iter = 1000,
  seed = 1
)

homo_results <- sweep_pcsets_k_ht(
  df = health_homo,
  scores = scores_homo,
  pc_sets = pc_sets,
  k_values = k_values,
  region_col = "region_2x2",
  nstart = 100,
  max_iter = 1000,
  seed = 1
)

print(het_results)
print(homo_results)

ggplot(het_results, aes(x = k, y = HTA, color = pc_set, group = pc_set)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(x = "k (number of clusters)", y = "HTA", color = "PC set") +
  theme_bw()

ggplot(homo_results, aes(x = k, y = HTA, color = pc_set, group = pc_set)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(x = "k (number of clusters)", y = "HTA", color = "PC set") +
  theme_bw()

# ----------------------------------------------------------------------------
# 3. Diagnose the unexpectedly high homogeneous-dataset HTA
# ----------------------------------------------------------------------------
pcs_use <- paste0("PC", 1:4)
k <- 3

set.seed(1)
km <- kmeans(
  x = scores_homo[, pcs_use, drop = FALSE],
  centers = k,
  nstart = 100,
  iter.max = 1000,
  algorithm = "Lloyd"
)

df2 <- health_homo
df2$trait_cluster <- factor(paste0("Cluster_", km$cluster))

tab <- with(df2, table(region_2x2, trait_cluster))
prop <- prop.table(tab, margin = 1)

print(prop)

region_mix <- t(apply(prop, 1, function(p) {
  p <- p[p > 0]
  c(
    entropy = -sum(p * log(p)) / log(length(levels(df2$trait_cluster))),
    max_p = max(p),
    min_p = min(p)
  )
}))
print(region_mix)

sizes <- km$size
print(c(
  min = min(sizes),
  max = max(sizes),
  ratio = max(sizes) / min(sizes),
  sd = sd(sizes)
))

Z <- scores_homo[, paste0("PC", 1:4), drop = FALSE]
print(mean(silhouette(km$cluster, dist(Z))[, 3]))

# ----------------------------------------------------------------------------
# 4. Random-region sensitivity check
# ----------------------------------------------------------------------------
# The original notebook created the random region label on df2 but then tried
# to run the sweep on health_homo, so that block could not run as written.
# Here the same random-label experiment is attached to the dataset passed to
# the sweep; no HTA or k-means calculation is changed.
set.seed(123)
health_homo$region_2x2_rand <- sample(health_homo$region_2x2)
health_homo$region_2x2_rand <- factor(
  health_homo$region_2x2_rand,
  levels = levels(health_homo$region_2x2)
)

homo_random_region_results <- sweep_pcsets_k_ht(
  df = health_homo,
  scores = scores_homo,
  pc_sets = pc_sets,
  k_values = k_values,
  region_col = "region_2x2_rand",
  nstart = 100,
  max_iter = 1000,
  seed = 1
)

print(homo_random_region_results)

# ----------------------------------------------------------------------------
# 5. Permutation null for the selected clustering configuration
# ----------------------------------------------------------------------------
set.seed(1)
B <- 500

obs <- compute_ht_cluster(
  df = df2,
  region_col = "region_2x2",
  trait_col = "trait_cluster"
)
hta_obs <- obs$HTA

hta_null <- replicate(B, {
  df_tmp <- df2
  df_tmp$region_perm <- sample(df_tmp$region_2x2)
  compute_ht_cluster(
    df = df_tmp,
    region_col = "region_perm",
    trait_col = "trait_cluster"
  )$HTA
})

mean_null <- mean(hta_null)
sd_null <- sd(hta_null)
z <- (hta_obs - mean_null) / sd_null
p_right <- mean(hta_null >= hta_obs)
p_two <- 2 * min(mean(hta_null >= hta_obs), mean(hta_null <= hta_obs))

print(list(
  hta_obs = hta_obs,
  mean_null = mean_null,
  sd_null = sd_null,
  z = z,
  p_right = p_right,
  p_two = p_two
))
