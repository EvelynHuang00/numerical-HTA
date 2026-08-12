# Main analysis retained in the final report:
# categorical Numerical HTA, pseudo-homogeneous comparison, and Alpha Simpson benchmark.
# Run from the repository root.

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)

source("R/hta_core.R")

dir.create("figures", showWarnings = FALSE)

# ----------------------------------------------------------------------------
# 1. Load data
# ----------------------------------------------------------------------------
AD_data <- read.csv("data/alzheimers_disease_data.csv", header = TRUE)
psedo_data <- read.csv("data/psedo_homogeneous_dataset.csv", header = TRUE)

ad_trait_cols <- c("Gender", "Ethnicity", "EducationLevel", "Smoking")
psedo_trait_cols <- c("SBP_grp", "LDL_grp")

# ----------------------------------------------------------------------------
# 2. Numerical HTA across Age x BMI grid resolutions
# ----------------------------------------------------------------------------
# The final report is descriptive, so permutation testing is not part of this
# public-facing script. The same permutation code is retained under exploratory/.
ks <- 2:6

hta_results <- lapply(ks, function(k) {
  compute_hta(
    AD_data,
    age_k = k,
    bmi_k = k,
    trait_cols = ad_trait_cols,
    min_cell_n = 5,
    perm_B = 0
  )
})
names(hta_results) <- paste0(ks, "x", ks)

hta_resolution <- data.frame(
  k = ks,
  HTA = vapply(hta_results, function(x) x$HTA, numeric(1)),
  cells_used = vapply(hta_results, function(x) x$cells_used, numeric(1)),
  total_n = vapply(hta_results, function(x) x$total_n, numeric(1))
)

print(hta_resolution)

fig_a1 <- ggplot(hta_resolution, aes(x = k, y = HTA)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = ks) +
  labs(
    x = "k (Age×BMI bins per axis)",
    y = "HTA"
  ) +
  theme_minimal()

print(fig_a1)
ggsave("figures/figure_A1_hta_resolution.png", fig_a1, width = 6.5, height = 4.5, dpi = 300)

# ----------------------------------------------------------------------------
# 3. Within-cell trait composition: heterogeneous Alzheimer cohort
# ----------------------------------------------------------------------------
res_2x2 <- hta_results[["2x2"]]
comp_2x2 <- make_comp_df(AD_data, res_2x2, ad_trait_cols, "2x2")

fig_a2 <- ggplot(comp_2x2, aes(x = "cell", y = prop, fill = trait)) +
  geom_col(position = "stack") +
  coord_flip() +
  facet_grid(BMIBin ~ AgeBin) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Trait composition within each Age×BMI cell (2×2 grid)",
    x = NULL,
    y = "Proportion within cell"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(angle = 0)
  )

print(fig_a2)
ggsave("figures/figure_A2_heterogeneous_composition.png", fig_a2, width = 10, height = 6.5, dpi = 300)

# ----------------------------------------------------------------------------
# 4. Pseudo-homogeneous comparison
# ----------------------------------------------------------------------------
psedo_ks <- 2:4

psedo_hta_results <- lapply(psedo_ks, function(k) {
  compute_hta(
    psedo_data,
    age_k = k,
    bmi_k = k,
    trait_cols = psedo_trait_cols,
    min_cell_n = 5,
    perm_B = 0
  )
})
names(psedo_hta_results) <- paste0(psedo_ks, "x", psedo_ks)

psedo_hta <- data.frame(
  k = psedo_ks,
  HTA = vapply(psedo_hta_results, function(x) x$HTA, numeric(1)),
  cells_used = vapply(psedo_hta_results, function(x) x$cells_used, numeric(1)),
  total_n = vapply(psedo_hta_results, function(x) x$total_n, numeric(1))
)

print(psedo_hta)

psedo_res_2x2 <- psedo_hta_results[["2x2"]]
psedo_comp_2x2 <- make_comp_df(
  psedo_data,
  psedo_res_2x2,
  trait_cols = psedo_trait_cols,
  "2x2"
)

fig_a3 <- ggplot(psedo_comp_2x2, aes(x = "cell", y = prop, fill = trait)) +
  geom_col(position = "stack") +
  coord_flip() +
  facet_grid(BMIBin ~ AgeBin) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Trait composition within each Age×BMI cell with synthetic homogeneous data",
    x = NULL,
    y = "Proportion within cell"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

print(fig_a3)
ggsave("figures/figure_A3_homogeneous_composition.png", fig_a3, width = 9, height = 6, dpi = 300)

# ----------------------------------------------------------------------------
# 5. Alpha Simpson benchmark under the same Age x BMI grids
# ----------------------------------------------------------------------------
benchmark_ks <- 2:4

ad_gs_results <- lapply(benchmark_ks, function(k) {
  compute_gs(
    AD_data,
    age_k = k,
    bmi_k = k,
    trait_cols = ad_trait_cols,
    min_cell_n = 5
  )
})

psedo_gs_results <- lapply(benchmark_ks, function(k) {
  compute_gs(
    psedo_data,
    age_k = k,
    bmi_k = k,
    trait_cols = psedo_trait_cols,
    min_cell_n = 5
  )
})

benchmark_df <- bind_rows(
  data.frame(
    dataset = "Heterogeneous",
    k = benchmark_ks,
    HTA = vapply(hta_results[paste0(benchmark_ks, "x", benchmark_ks)], function(x) x$HTA, numeric(1)),
    Alpha_Simpson = vapply(ad_gs_results, function(x) x$GS_HTA, numeric(1))
  ),
  data.frame(
    dataset = "Homogeneous",
    k = benchmark_ks,
    HTA = vapply(psedo_hta_results, function(x) x$HTA, numeric(1)),
    Alpha_Simpson = vapply(psedo_gs_results, function(x) x$GS_HTA, numeric(1))
  )
) %>%
  mutate(
    grid = factor(paste0(k, "×", k), levels = c("2×2", "3×3", "4×4")),
    dataset = factor(dataset, levels = c("Heterogeneous", "Homogeneous"))
  )

print(benchmark_df)

plot_long <- benchmark_df %>%
  select(dataset, grid, HTA, Alpha_Simpson) %>%
  pivot_longer(
    cols = c(HTA, Alpha_Simpson),
    names_to = "Metric",
    values_to = "Score"
  )

fig_a4 <- ggplot(plot_long, aes(x = grid, y = Score, group = Metric, shape = Metric)) +
  geom_line(linewidth = 0.8, color = "black") +
  geom_point(size = 3, color = "black") +
  facet_wrap(~dataset, nrow = 1) +
  labs(
    x = "Age × BMI grid",
    y = "Score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 16),
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

print(fig_a4)
ggsave("figures/figure_A4_hta_vs_alpha_simpson.png", fig_a4, width = 9, height = 5.5, dpi = 300)
