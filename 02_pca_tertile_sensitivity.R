# PCA-tertile sensitivity analysis retained in the final report.
# Age and BMI define regions; PC1 from the remaining standardized numeric
# variables is divided into Low/Mid/High tertiles and used as the trait label.
# Run from the repository root.

library(dplyr)

source("R/hta_core.R")
source("R/pca_helpers.R")

# ----------------------------------------------------------------------------
# 1. Load and clean the numerically dominated cohort
# ----------------------------------------------------------------------------
health_data <- read.csv("data/dirty_v3_path.csv", header = TRUE)
health_data <- na.omit(health_data)
health_data <- health_data %>% select(-c(random_notes, noise_col))

# The original notebook constructs this comparison cohort from the cleaned data.
health_homo <- make_pseudo_homogeneous_numeric(
  health_data,
  n = 5000,
  seed = 123
)

# ----------------------------------------------------------------------------
# 2. Heterogeneous dataset: Age x BMI regions + PCA excluding Age/BMI
# ----------------------------------------------------------------------------
health_data <- make_age_bmi_2x2(health_data)

pca_fit <- fit_numeric_pca_excluding_age_bmi(health_data)
scores <- as.data.frame(pca_fit$x)

pc1_tertile <- cut(
  scores$PC1,
  breaks = quantile(
    scores$PC1,
    probs = c(0, 1 / 3, 2 / 3, 1),
    na.rm = TRUE,
    type = 7
  ),
  include.lowest = TRUE,
  labels = c("Low", "Mid", "High")
)

health_data$trait_pc1_tertile <- factor(
  pc1_tertile,
  levels = c("Low", "Mid", "High")
)

res_pc1_tertile <- compute_ht_fixed_levels(
  health_data,
  region_col = "region_2x2",
  trait_col = "trait_pc1_tertile"
)

cat("Heterogeneous dataset HTA:\n")
print(res_pc1_tertile$HTA)
print(res_pc1_tertile$local_hti)

cat("\nPC1 tertile proportions within each Age x BMI region:\n")
print(
  prop.table(
    table(health_data$region_2x2, health_data$trait_pc1_tertile),
    margin = 1
  )
)

# ----------------------------------------------------------------------------
# 3. Pseudo-homogeneous dataset under the same analysis design
# ----------------------------------------------------------------------------
health_homo <- make_age_bmi_2x2(health_homo)

pca_fit_homo <- fit_numeric_pca_excluding_age_bmi(health_homo)
scores_homo <- as.data.frame(pca_fit_homo$x)

pc1_tertile_homo <- cut(
  scores_homo$PC1,
  breaks = quantile(
    scores_homo$PC1,
    probs = c(0, 1 / 3, 2 / 3, 1),
    na.rm = TRUE,
    type = 7
  ),
  include.lowest = TRUE,
  labels = c("Low", "Mid", "High")
)

health_homo$trait_pc1_tertile <- factor(
  pc1_tertile_homo,
  levels = c("Low", "Mid", "High")
)

res_pc1_tertile_homo <- compute_ht_fixed_levels(
  health_homo,
  region_col = "region_2x2",
  trait_col = "trait_pc1_tertile"
)

cat("\nPseudo-homogeneous dataset HTA:\n")
print(res_pc1_tertile_homo$HTA)
print(res_pc1_tertile_homo$local_hti)

cat("\nPC1 tertile proportions within each Age x BMI region:\n")
print(
  prop.table(
    table(health_homo$region_2x2, health_homo$trait_pc1_tertile),
    margin = 1
  )
)

pca_tertile_summary <- data.frame(
  dataset = c("Heterogeneous", "Pseudo-homogeneous"),
  HTA = c(res_pc1_tertile$HTA, res_pc1_tertile_homo$HTA)
)

cat("\nPCA-tertile sensitivity summary:\n")
print(pca_tertile_summary)
