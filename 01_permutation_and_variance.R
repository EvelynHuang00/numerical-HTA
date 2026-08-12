# Exploratory permutation and variance analyses from numerical_HTA.Rmd.
# These blocks were part of the research process but were not central to the
# descriptive result set retained in the final report.
# Run from the repository root.

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)

source("R/hta_core.R")

AD_data <- read.csv("data/alzheimers_disease_data.csv", header = TRUE)
psedo_data <- read.csv("data/psedo_homogeneous_dataset.csv", header = TRUE)

ad_trait_cols <- c("Gender", "Ethnicity", "EducationLevel", "Smoking")
psedo_trait_cols <- c("SBP_grp", "LDL_grp")

# ----------------------------------------------------------------------------
# 1. Early single-variable HTI checks
# ----------------------------------------------------------------------------
min(AD_data$Age)
max(AD_data$Age)

AD_data$age_group <- cut(
  AD_data$Age,
  breaks = c(60, 70, 80, 90),
  labels = c("1", "2", "3"),
  include.lowest = TRUE
)

p <- prop.table(table(AD_data$age_group))
C <- length(p)
HTI_age <- -sum(p * log(p)) / log(C)
HTI_age

min(AD_data$BMI)
max(AD_data$BMI)

AD_data$bmi_group <- cut(
  AD_data$BMI,
  breaks = c(0, 18.5, 25, 30, 100),
  labels = c("under", "normal", "over", "obese"),
  include.lowest = TRUE
)

p <- prop.table(table(AD_data$bmi_group))
C <- length(p)
HTI_bmi <- -sum(p * log(p)) / log(C)
HTI_bmi

# ----------------------------------------------------------------------------
# 2. Cell-z permutation analysis across grid resolutions
# ----------------------------------------------------------------------------
res_2x2_cellz <- compute_hta(
  AD_data,
  age_k = 2,
  bmi_k = 2,
  perm_B = 500,
  trait_cols = ad_trait_cols,
  p_method = "cell_z",
  debug = TRUE
)

res_3x3_cellz <- compute_hta(
  AD_data,
  age_k = 3,
  bmi_k = 3,
  perm_B = 500,
  trait_cols = ad_trait_cols,
  p_method = "cell_z",
  debug = TRUE
)

res_4x4_cellz <- compute_hta(
  AD_data,
  age_k = 4,
  bmi_k = 4,
  perm_B = 500,
  trait_cols = ad_trait_cols,
  p_method = "cell_z",
  debug = TRUE
)

res_5x5_cellz <- compute_hta(
  AD_data,
  age_k = 5,
  bmi_k = 5,
  perm_B = 500,
  trait_cols = ad_trait_cols,
  p_method = "cell_z",
  debug = TRUE
)

res_6x6_cellz <- compute_hta(
  AD_data,
  age_k = 6,
  bmi_k = 6,
  perm_B = 500,
  trait_cols = ad_trait_cols,
  p_method = "cell_z",
  debug = TRUE
)

list(
  `2x2` = list(HTA = res_2x2_cellz$HTA, p_val = res_2x2_cellz$perm$p_one),
  `3x3` = list(HTA = res_3x3_cellz$HTA, p_val = res_3x3_cellz$perm$p_one),
  `4x4` = list(HTA = res_4x4_cellz$HTA, p_val = res_4x4_cellz$perm$p_one),
  `5x5` = list(HTA = res_5x5_cellz$HTA, p_val = res_5x5_cellz$perm$p_one),
  `6x6` = list(HTA = res_6x6_cellz$HTA, p_val = res_6x6_cellz$perm$p_one)
)

# The original notebook contained the following block and explicitly noted
# that it no longer worked because cell_z does not store null_hta:
# hist(is.numeric(res_4x4$perm$null_hta), breaks = 40)
# abline(v = res_4x4$HTA, col = "red", lwd = 2)

# ----------------------------------------------------------------------------
# 3. Detailed permutation diagnostics for the 3x3 grid
# ----------------------------------------------------------------------------
out <- compute_hta(
  AD_data,
  age_k = 3,
  bmi_k = 3,
  perm_B = 500,
  trait_cols = ad_trait_cols,
  p_method = "cell_z",
  debug = TRUE
)

out$perm$p_one
out$diag$perm_B
out$diag$n_unique_traits_kept
summary(out$diag$w_perm)
summary(out$diag$Hcheck_sd)

out2 <- compute_hta(
  AD_data,
  age_k = 3,
  bmi_k = 3,
  perm_B = 500,
  trait_cols = ad_trait_cols,
  p_method = "hta_tail",
  debug = TRUE
)

out2$perm$p_one

hist(
  out2$perm$null_hta,
  breaks = 40,
  main = "HTA from observed data",
  xlab = "HTA under permutation"
)
abline(v = out$HTA, lwd = 2, col = "red")

# ----------------------------------------------------------------------------
# 4. Pseudo-homogeneous permutation experiments
# ----------------------------------------------------------------------------
psedo_res_2x2_cellz <- compute_hta(
  psedo_data,
  age_k = 2,
  bmi_k = 2,
  perm_B = 500,
  trait_cols = psedo_trait_cols,
  p_method = "cell_z",
  debug = TRUE
)

psedo_res_3x3_cellz <- compute_hta(
  psedo_data,
  age_k = 3,
  bmi_k = 3,
  perm_B = 500,
  trait_cols = psedo_trait_cols,
  p_method = "cell_z",
  debug = TRUE
)

psedo_res_4x4_cellz <- compute_hta(
  psedo_data,
  age_k = 4,
  bmi_k = 4,
  perm_B = 500,
  trait_cols = psedo_trait_cols,
  p_method = "cell_z",
  debug = TRUE
)

list(
  `2x2` = list(HTA = psedo_res_2x2_cellz$HTA, p_val = psedo_res_2x2_cellz$perm$p_one),
  `3x3` = list(HTA = psedo_res_3x3_cellz$HTA, p_val = psedo_res_3x3_cellz$perm$p_one),
  `4x4` = list(HTA = psedo_res_4x4_cellz$HTA, p_val = psedo_res_4x4_cellz$perm$p_one)
)

psedo_res_2x2_tail <- compute_hta(
  psedo_data,
  age_k = 2,
  bmi_k = 2,
  perm_B = 500,
  trait_cols = psedo_trait_cols,
  p_method = "hta_tail",
  debug = TRUE
)

psedo_res_3x3_tail <- compute_hta(
  psedo_data,
  age_k = 3,
  bmi_k = 3,
  perm_B = 500,
  trait_cols = psedo_trait_cols,
  p_method = "hta_tail",
  debug = TRUE
)

psedo_res_4x4_tail <- compute_hta(
  psedo_data,
  age_k = 4,
  bmi_k = 4,
  perm_B = 500,
  trait_cols = psedo_trait_cols,
  p_method = "hta_tail",
  debug = TRUE
)

list(
  `2x2` = list(HTA = psedo_res_2x2_tail$HTA, p_val = psedo_res_2x2_tail$perm$p_one),
  `3x3` = list(HTA = psedo_res_3x3_tail$HTA, p_val = psedo_res_3x3_tail$perm$p_one),
  `4x4` = list(HTA = psedo_res_4x4_tail$HTA, p_val = psedo_res_4x4_tail$perm$p_one)
)

# ----------------------------------------------------------------------------
# 5. Tau-squared / mixed-model exploration
# ----------------------------------------------------------------------------
df_2 <- make_cells_quantile(AD_data, K = 2)
df_3 <- make_cells_quantile(AD_data, K = 3)
df_4 <- make_cells_quantile(AD_data, K = 4)

# Smoking
m_smoke_2 <- glmer(
  Smoking ~ 1 + (1 | cell_id),
  data = df_2,
  family = binomial
)
VarCorr(m_smoke_2)

m_smoke_3 <- glmer(
  Smoking ~ 1 + (1 | cell_id),
  data = df_3,
  family = binomial
)
VarCorr(m_smoke_3)

m_smoke_4 <- glmer(
  Smoking ~ 1 + (1 | cell_id),
  data = df_4,
  family = binomial
)
VarCorr(m_smoke_4)

# Gender
m_gender_2 <- glmer(
  Gender ~ 1 + (1 | cell_id),
  data = df_2,
  family = binomial
)
VarCorr(m_gender_2)

m_gender_3 <- glmer(
  Gender ~ 1 + (1 | cell_id),
  data = df_3,
  family = binomial
)
VarCorr(m_gender_3)

m_gender_4 <- glmer(
  Gender ~ 1 + (1 | cell_id),
  data = df_4,
  family = binomial
)
VarCorr(m_gender_4)

# Ethnicity and EducationLevel
df_3 <- df_3 %>%
  mutate(
    Ethnicity_num = as.numeric(Ethnicity),
    EducationLevel_num = as.numeric(EducationLevel),
    Ethnicity_z = as.numeric(scale(Ethnicity_num)),
    EducationLevel_z = as.numeric(scale(EducationLevel_num))
  )

m_eth <- lmer(
  Ethnicity_z ~ 1 + (1 | cell_id),
  data = df_3,
  REML = TRUE
)
VarCorr(m_eth)

m_edu <- lmer(
  EducationLevel_z ~ 1 + (1 | cell_id),
  data = df_3,
  REML = TRUE
)
VarCorr(m_edu)

vc_eth <- VarCorr(m_eth)
tau2_eth <- as.numeric(vc_eth$cell_id)^2
sigma2_eth <- attr(vc_eth, "sc")^2
ICC_eth <- tau2_eth / (tau2_eth + sigma2_eth)

vc_edu <- VarCorr(m_edu)
tau2_edu <- as.numeric(vc_edu$cell_id)^2
sigma2_edu <- attr(vc_edu, "sc")^2
ICC_edu <- tau2_edu / (tau2_edu + sigma2_edu)

tau2_eth
ICC_eth
tau2_edu
ICC_edu

# Numerical variable: SleepQuality
df_2 <- df_2 %>% mutate(SleepQuality_z = as.numeric(scale(SleepQuality)))
df_3 <- df_3 %>% mutate(SleepQuality_z = as.numeric(scale(SleepQuality)))
df_4 <- df_4 %>% mutate(SleepQuality_z = as.numeric(scale(SleepQuality)))

# Preserve the original model calls, which use SleepQuality rather than
# the newly created SleepQuality_z.
m_sp_2 <- glmer(
  SleepQuality ~ 1 + (1 | cell_id),
  data = df_2
)
VarCorr(m_sp_2)

m_sp_3 <- glmer(
  SleepQuality ~ 1 + (1 | cell_id),
  data = df_3
)
VarCorr(m_sp_3)

m_sp_4 <- glmer(
  SleepQuality ~ 1 + (1 | cell_id),
  data = df_4
)
VarCorr(m_sp_4)
