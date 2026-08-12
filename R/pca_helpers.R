# Helper functions for the PCA-based HTA extension.
# These preserve the calculations used in HTA_PCA.Rmd while removing repeated setup code.

# Construct the 2x2 Age x BMI regions used in the PCA experiments.
make_age_bmi_2x2 <- function(df) {
  age_q <- quantile(df$Age, probs = 0.5, na.rm = TRUE, type = 7)
  bmi_q <- quantile(df$BMI, probs = 0.5, na.rm = TRUE, type = 7)

  df$age_bin <- ifelse(df$Age <= age_q, "Age_Low", "Age_High")
  df$bmi_bin <- ifelse(df$BMI <= bmi_q, "BMI_Low", "BMI_High")

  df$region_2x2 <- paste(df$age_bin, df$bmi_bin, sep = "__")
  df$region_2x2 <- factor(
    df$region_2x2,
    levels = c(
      "Age_Low__BMI_Low",
      "Age_Low__BMI_High",
      "Age_High__BMI_Low",
      "Age_High__BMI_High"
    )
  )

  df
}

# PCA pipeline used in the retained extension:
# exclude Age/BMI, keep numeric variables, standardize, then run prcomp().
fit_numeric_pca_excluding_age_bmi <- function(df) {
  X <- subset(df, select = -c(Age, BMI))
  X_num <- X[, sapply(X, is.numeric), drop = FALSE]
  X_scaled <- scale(X_num, center = TRUE, scale = TRUE)
  prcomp(X_scaled, center = FALSE, scale. = FALSE)
}

# HTA implementation used in the PCA notebook.
# Important: C is the number of trait levels in the full analyzed dataset,
# matching the original compute_ht() implementation in HTA_PCA.Rmd.
compute_ht_fixed_levels <- function(df,
                                    region_col = "region_2x2",
                                    trait_col = "trait") {
  ok <- !is.na(df[[region_col]]) & !is.na(df[[trait_col]])
  d <- df[ok, , drop = FALSE]

  trait_levels <- levels(d[[trait_col]])
  if (is.null(trait_levels) || anyNA(trait_levels)) {
    trait_levels <- sort(unique(d[[trait_col]]))
  }
  C <- length(trait_levels)
  if (C <= 1) stop("Need >= 2 trait categories to compute HTI/HTA.")

  hti_region <- function(x) {
    tab <- table(factor(x, levels = trait_levels))
    p <- as.numeric(tab) / sum(tab)
    p <- p[p > 0]
    -sum(p * log(p)) / log(C)
  }

  regions <- levels(d[[region_col]])
  if (is.null(regions) || anyNA(regions)) {
    regions <- sort(unique(d[[region_col]]))
  }

  n_by_region <- sapply(regions, function(r) sum(d[[region_col]] == r))
  local_hti <- sapply(
    regions,
    function(r) hti_region(d[[trait_col]][d[[region_col]] == r])
  )
  weights <- n_by_region / sum(n_by_region)
  HTA <- sum(weights * local_hti)

  list(
    HTA = HTA,
    local_hti = local_hti,
    n_by_region = n_by_region,
    weights = weights
  )
}

# Earlier PCA experiment: PC-defined regions and PC-defined traits.
# Preserved for exploratory_analysis/02_pca_alternative_traits.R.
compute_hta_pca <- function(df, pca_fit,
                            region_pcs = c(1, 2), region_k = 2,
                            trait_pcs = c(3, 4, 5), trait_k = 3,
                            min_cell_n = 5) {
  scores <- as.data.frame(pca_fit$x)

  if (nrow(df) != nrow(scores)) {
    stop("nrow(df) != nrow(pca_fit$x). Fit PCA on the same df rows you pass here.")
  }

  for (j in seq_len(ncol(scores))) {
    df[[paste0("PC", j)]] <- scores[[j]]
  }

  pcA <- paste0("PC", region_pcs[1])
  pcB <- paste0("PC", region_pcs[2])

  brA <- qbreaks(df[[pcA]], region_k)
  brB <- qbreaks(df[[pcB]], region_k)

  df$Reg1Bin <- cut(df[[pcA]], breaks = brA, include.lowest = TRUE, right = TRUE)
  df$Reg2Bin <- cut(df[[pcB]], breaks = brB, include.lowest = TRUE, right = TRUE)
  df$cell <- interaction(df$Reg1Bin, df$Reg2Bin, drop = TRUE)

  trait_bin_cols <- character(length(trait_pcs))
  trait_breaks <- vector("list", length(trait_pcs))

  for (i in seq_along(trait_pcs)) {
    pcname <- paste0("PC", trait_pcs[i])
    br <- qbreaks(df[[pcname]], trait_k)
    trait_breaks[[i]] <- br

    colnm <- paste0(pcname, "Bin")
    df[[colnm]] <- cut(df[[pcname]], breaks = br, include.lowest = TRUE, right = TRUE)
    df[[colnm]] <- factor(df[[colnm]])
    trait_bin_cols[i] <- colnm
  }
  names(trait_breaks) <- paste0("PC", trait_pcs)

  complete_traits <- stats::complete.cases(df[, trait_bin_cols, drop = FALSE])
  traits_all <- interaction(df[, trait_bin_cols, drop = FALSE], drop = TRUE)
  traits <- as.character(traits_all[complete_traits])
  cells <- droplevels(df$cell[complete_traits])

  by_cell <- split(traits, cells)
  cell_n <- sapply(by_cell, length)
  keep <- cell_n >= min_cell_n
  by_cell <- by_cell[keep]
  cell_n <- cell_n[keep]

  if (length(by_cell) == 0) {
    stop("No cells left after filtering; use coarser region_k or lower min_cell_n.")
  }

  local_hti <- sapply(by_cell, hti)
  w <- cell_n
  HTA <- sum(w * local_hti) / sum(w)

  list(
    HTA = HTA,
    cells_used = length(by_cell),
    total_n = sum(w),
    region_breaks = list(Reg1 = brA, Reg2 = brB),
    trait_breaks = trait_breaks,
    local_hti = local_hti,
    cell_n = w,
    cells_used_names = names(by_cell)
  )
}

# Synthetic more-homogeneous cohort used in HTA_PCA.Rmd and HTA_clustering.Rmd.
# The generation sequence and parameter choices are unchanged.
make_pseudo_homogeneous_numeric <- function(health_data, n = 5000, seed = 123) {
  set.seed(seed)

  health_homo <- health_data[rep(1, n), , drop = FALSE]

  tab <- sort(table(health_data$Gender, useNA = "no"), decreasing = TRUE)
  lv <- names(tab)
  p <- rep(0, length(lv))
  p[1] <- 0.8
  if (length(lv) >= 2) p[2] <- 0.1
  if (length(lv) > 2) p[3:length(lv)] <- 0
  health_homo$Gender <- sample(lv, n, replace = TRUE, prob = p / sum(p))

  tab <- sort(table(health_data$Medical.Condition, useNA = "no"), decreasing = TRUE)
  lv <- names(tab)
  p <- rep(0, length(lv))
  p[1] <- 0.85
  p[2] <- 0.10
  if (length(lv) > 2) p[3:length(lv)] <- 0.05 / (length(lv) - 2)
  health_homo$Medical.Condition <- sample(lv, n, replace = TRUE, prob = p)

  mu <- mean(health_data$Age, na.rm = TRUE); s <- sd(health_data$Age, na.rm = TRUE)
  health_homo$Age <- rnorm(n, mu, 0.5 * s)

  mu <- mean(health_data$Glucose, na.rm = TRUE); s <- sd(health_data$Glucose, na.rm = TRUE)
  health_homo$Glucose <- rnorm(n, mu, 0.5 * s)

  mu <- mean(health_data$Blood.Pressure, na.rm = TRUE); s <- sd(health_data$Blood.Pressure, na.rm = TRUE)
  health_homo$Blood.Pressure <- rnorm(n, mu, 0.5 * s)

  mu <- mean(health_data$BMI, na.rm = TRUE); s <- sd(health_data$BMI, na.rm = TRUE)
  health_homo$BMI <- rnorm(n, mu, 0.5 * s)

  mu <- mean(health_data$Oxygen.Saturation, na.rm = TRUE); s <- sd(health_data$Oxygen.Saturation, na.rm = TRUE)
  health_homo$Oxygen.Saturation <- rnorm(n, mu, 0.5 * s)

  mu <- mean(health_data$LengthOfStay, na.rm = TRUE); s <- sd(health_data$LengthOfStay, na.rm = TRUE)
  health_homo$LengthOfStay <- as.integer(round(rnorm(n, mu, 0.5 * s)))

  mu <- mean(health_data$Cholesterol, na.rm = TRUE); s <- sd(health_data$Cholesterol, na.rm = TRUE)
  health_homo$Cholesterol <- rnorm(n, mu, 0.5 * s)

  mu <- mean(health_data$Triglycerides, na.rm = TRUE); s <- sd(health_data$Triglycerides, na.rm = TRUE)
  health_homo$Triglycerides <- rnorm(n, mu, 0.5 * s)

  mu <- mean(health_data$HbA1c, na.rm = TRUE); s <- sd(health_data$HbA1c, na.rm = TRUE)
  health_homo$HbA1c <- rnorm(n, mu, 0.5 * s)

  p1 <- mean(health_data$Smoking == 1, na.rm = TRUE)
  health_homo$Smoking <- rbinom(n, 1, prob = min(max(p1, 0.2), 0.8))

  p1 <- mean(health_data$Alcohol == 1, na.rm = TRUE)
  health_homo$Alcohol <- rbinom(n, 1, prob = min(max(p1, 0.2), 0.8))

  mu <- mean(health_data$Physical.Activity, na.rm = TRUE); s <- sd(health_data$Physical.Activity, na.rm = TRUE)
  health_homo$Physical.Activity <- rnorm(n, mu, 0.5 * s)

  mu <- mean(health_data$Diet.Score, na.rm = TRUE); s <- sd(health_data$Diet.Score, na.rm = TRUE)
  health_homo$Diet.Score <- rnorm(n, mu, 0.5 * s)

  p1 <- mean(health_data$Family.History == 1, na.rm = TRUE)
  health_homo$Family.History <- rbinom(n, 1, prob = min(max(p1, 0.2), 0.8))

  mu <- mean(health_data$Stress.Level, na.rm = TRUE); s <- sd(health_data$Stress.Level, na.rm = TRUE)
  health_homo$Stress.Level <- rnorm(n, mu, 0.5 * s)

  mu <- mean(health_data$Sleep.Hours, na.rm = TRUE); s <- sd(health_data$Sleep.Hours, na.rm = TRUE)
  health_homo$Sleep.Hours <- rnorm(n, mu, 0.5 * s)

  # Same cleanup used in the original notebook if these exploratory columns exist.
  health_homo <- health_homo[, !(names(health_homo) %in% c("region", "trait")), drop = FALSE]

  health_homo
}
