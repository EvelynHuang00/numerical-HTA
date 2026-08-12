# Numerical HTA helper functions
# Extracted from numerical_HTA.Rmd without changing the underlying calculations.

# ---- HTI / HTA utilities ----
hti <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(0)
  if (is.factor(x)) x <- droplevels(x)
  tab <- table(x)
  p <- as.numeric(tab) / sum(tab)
  p <- p[p > 0]
  C <- length(p)
  if (C <= 1) return(0)
  -sum(p * log(p)) / log(C)
}

qbreaks <- function(v, k) {
  br <- unique(quantile(v, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE))
  if (length(br) < k + 1) br <- unique(pretty(v, n = k))
  br
}

compute_hta <- function(df, age_k, bmi_k,
                        trait_cols,
                        min_cell_n = 5, perm_B = 0, seed = 1,
                        p_method = c("hta_tail", "cell_z"),
                        debug = FALSE) {

  p_method <- match.arg(p_method)

  # Coerce traits to factors
  for (nm in trait_cols) {
    if (!is.factor(df[[nm]])) df[[nm]] <- factor(df[[nm]])
  }

  # Build Age/BMI bins
  age_br <- qbreaks(df$Age, age_k)
  bmi_br <- qbreaks(df$BMI, bmi_k)
  df$AgeBin <- cut(df$Age, breaks = age_br, include.lowest = TRUE, right = TRUE)
  df$BMIBin <- cut(df$BMI, breaks = bmi_br, include.lowest = TRUE, right = TRUE)
  df$cell   <- interaction(df$AgeBin, df$BMIBin, drop = TRUE)

  # Trait combinations on complete rows only
  complete_traits <- stats::complete.cases(df[, trait_cols])
  traits_all <- interaction(df[, trait_cols], drop = TRUE)
  traits     <- as.character(traits_all[complete_traits])
  cells      <- droplevels(df$cell[complete_traits])

  # Split by cell and enforce minimum cell size
  by_cell <- split(traits, cells)
  cell_n  <- sapply(by_cell, length)
  keep    <- cell_n >= min_cell_n
  by_cell <- by_cell[keep]
  cell_n  <- cell_n[keep]

  if (length(by_cell) == 0) {
    stop("No cells left after filtering; use a coarser grid or lower min_cell_n.")
  }

  # Local HTIs and weights
  local_hti <- sapply(by_cell, hti)
  w <- cell_n
  HTA <- sum(w * local_hti) / sum(w)

  res <- list(
    HTA        = HTA,
    cells_used = length(by_cell),
    total_n    = sum(w),
    age_breaks = age_br,
    bmi_breaks = bmi_br,
    local_hti  = local_hti,
    cell_n     = w,
    cells_used_names = names(by_cell)
  )

  # Permutation section (kept-only pool)
  if (perm_B > 0) {
    set.seed(seed)

    sp_idx_all <- split(seq_along(traits), cells)
    sp_idx      <- sp_idx_all[names(by_cell)]
    idx_kept    <- unlist(sp_idx, use.names = FALSE)
    traits_kept <- traits[idx_kept]

    w_perm <- vapply(sp_idx, length, integer(1))
    ends   <- cumsum(w_perm)
    starts <- c(1, head(ends, -1) + 1)

    if (p_method == "hta_tail") {
      hta_perm <- numeric(perm_B)
      for (b in seq_len(perm_B)) {
        perm <- sample(traits_kept)
        local_perm <- mapply(function(s, e) hti(perm[s:e]), starts, ends)
        hta_perm[b] <- sum(w_perm * local_perm) / sum(w_perm)
      }
      p_one <- mean(hta_perm >= HTA)
      res$perm <- list(B = perm_B, null_hta = hta_perm, p_one = p_one)

    } else { # "cell_z"
      H <- matrix(NA_real_, nrow = perm_B, ncol = length(w_perm))
      for (b in seq_len(perm_B)) {
        perm <- sample(traits_kept)
        H[b, ] <- mapply(function(s, e) hti(perm[s:e]), starts, ends)
      }
      mu_r <- colMeans(H)
      sd_r <- apply(H, 2, stats::sd)

      keep_r <- which(sd_r > 0 & is.finite(sd_r))
      if (length(keep_r) == 0L) {
        hta_perm <- as.vector(H %*% w_perm) / sum(w_perm)
        Z <- NA_real_
        p_one <- mean(hta_perm >= HTA)
      } else {
        num <- sum(local_hti[keep_r] - mu_r[keep_r])
        den <- sqrt(sum(sd_r[keep_r]^2))
        Z <- num / den
        p_one <- 1 - pnorm(Z)
      }
      res$perm <- list(
        B = perm_B,
        Z = Z,
        p_one = p_one,
        mu_r = mu_r,
        sd_r = sd_r,
        used_regions = keep_r
      )
    }

    if (debug) {
      Bcheck <- min(50L, perm_B)
      Hcheck <- matrix(NA_real_, nrow = Bcheck, ncol = length(w_perm))
      for (b in seq_len(Bcheck)) {
        perm <- sample(traits_kept)
        Hcheck[b, ] <- mapply(function(s, e) hti(perm[s:e]), starts, ends)
      }
      res$diag <- list(
        perm_B = perm_B,
        n_unique_traits_kept = length(unique(traits_kept)),
        w_perm = w_perm,
        Hcheck_sd = apply(Hcheck, 2, sd)
      )
    }
  }

  res
}

# ---- Trait-composition utilities ----
make_comp_df <- function(df, res, trait_cols, label) {
  df2 <- df %>%
    mutate(
      AgeBin = cut(
        Age,
        breaks = res$age_breaks,
        include.lowest = TRUE,
        right = TRUE
      ),
      BMIBin = cut(
        BMI,
        breaks = res$bmi_breaks,
        include.lowest = TRUE,
        right = TRUE
      )
    )

  df2$trait <- interaction(df2[, trait_cols], drop = TRUE)

  df2 %>%
    filter(!is.na(AgeBin), !is.na(BMIBin), !is.na(trait)) %>%
    group_by(AgeBin, BMIBin, trait) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(AgeBin, BMIBin) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    mutate(grid = label)
}

# ---- Quantile-cell utility used in exploratory tau-squared analysis ----
make_cells_quantile <- function(df, K) {
  probs <- seq(0, 1, length.out = K + 1)

  age_breaks <- quantile(df$Age, probs = probs, na.rm = TRUE)
  age_breaks <- unique(age_breaks)

  bmi_breaks <- quantile(df$BMI, probs = probs, na.rm = TRUE)
  bmi_breaks <- unique(bmi_breaks)

  df %>%
    mutate(
      age_bin = cut(Age, breaks = age_breaks, include.lowest = TRUE),
      bmi_bin = cut(BMI, breaks = bmi_breaks, include.lowest = TRUE),
      cell_id = interaction(age_bin, bmi_bin, drop = TRUE)
    )
}

# ---- Gini-Simpson utilities ----
gs_local <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  if (is.factor(x)) x <- droplevels(x)
  tab <- table(x)
  p <- as.numeric(tab) / sum(tab)
  if (length(p) <= 1) return(0)
  1 - sum(p^2)
}

compute_gs <- function(df, age_k, bmi_k,
                       trait_cols,
                       min_cell_n = 5) {

  # Coerce traits to factors
  for (nm in trait_cols) {
    if (!is.factor(df[[nm]])) df[[nm]] <- factor(df[[nm]])
  }

  # Build Age/BMI bins using the same logic as compute_hta
  age_br <- qbreaks(df$Age, age_k)
  bmi_br <- qbreaks(df$BMI, bmi_k)
  df$AgeBin <- cut(df$Age, breaks = age_br, include.lowest = TRUE, right = TRUE)
  df$BMIBin <- cut(df$BMI, breaks = bmi_br, include.lowest = TRUE, right = TRUE)
  df$cell   <- interaction(df$AgeBin, df$BMIBin, drop = TRUE)

  # Trait combinations on complete rows only
  complete_traits <- stats::complete.cases(df[, trait_cols])
  traits_all <- interaction(df[, trait_cols], drop = TRUE)
  traits <- as.character(traits_all[complete_traits])
  cells  <- droplevels(df$cell[complete_traits])

  # Split by cell and enforce minimum size
  by_cell <- split(traits, cells)
  cell_n  <- sapply(by_cell, length)
  keep    <- cell_n >= min_cell_n
  by_cell <- by_cell[keep]
  cell_n  <- cell_n[keep]

  if (length(by_cell) == 0) {
    stop("No cells left after filtering; use a coarser grid or lower min_cell_n.")
  }

  # Local Gini-Simpson and HTA-style weighted average
  local_gs <- sapply(by_cell, gs_local)
  w <- cell_n
  GS_HTA <- sum(w * local_gs) / sum(w)

  list(
    GS_HTA      = GS_HTA,
    cells_used  = length(by_cell),
    total_n     = sum(w),
    age_breaks  = age_br,
    bmi_breaks  = bmi_br,
    local_gs    = local_gs,
    cell_n      = w,
    cells_used_names = names(by_cell)
  )
}
