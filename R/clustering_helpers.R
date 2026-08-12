# Helper functions for the k-means HTA experiments.
# Extracted from HTA_clustering.Rmd without changing the calculation logic.

# HTA implementation used in the clustering notebook.
# This version keeps the notebook's fixed global trait-level normalization and
# atomic-vector coercion used to avoid list-column issues.
compute_ht_cluster <- function(df, region_col = "region_2x2", trait_col = "trait") {

  ok <- !is.na(df[[region_col]]) & !is.na(df[[trait_col]])
  d <- df[ok, , drop = FALSE]

  to_atomic_chr <- function(x) {
    if (is.list(x)) {
      return(vapply(x, function(z) as.character(z)[1], character(1)))
    }
    as.character(x)
  }

  d[[region_col]] <- factor(to_atomic_chr(d[[region_col]]))
  d[[trait_col]]  <- factor(to_atomic_chr(d[[trait_col]]))

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

# Run one k-means fit and compute HTA from the resulting cluster labels.
run_kmeans_ht <- function(df, scores, pcs_use, k,
                          region_col = "region_2x2",
                          nstart = 100, max_iter = 1000, seed = 1) {
  stopifnot(all(pcs_use %in% colnames(scores)))
  stopifnot(nrow(df) == nrow(scores))
  stopifnot(region_col %in% names(df))

  set.seed(seed)
  km <- kmeans(
    x = scores[, pcs_use, drop = FALSE],
    centers = k,
    nstart = nstart,
    iter.max = max_iter
  )

  df2 <- df
  df2$trait_cluster <- factor(paste0("Cluster_", km$cluster))

  res <- tryCatch(
    compute_ht_cluster(
      df2,
      region_col = region_col,
      trait_col = "trait_cluster"
    ),
    error = function(e) {
      list(HTA = NA_real_, local_hti = NA, .err = e$message)
    }
  )

  list(
    HTA = res$HTA,
    local_hti = res$local_hti,
    iter = km$iter,
    tot_withinss = km$tot.withinss,
    size = km$size,
    err_msg = res$.err
  )
}

# Sweep k for one fixed PC set.
sweep_k_ht <- function(df, scores, pcs_use, k_values,
                       region_col = "region_2x2",
                       nstart = 100, max_iter = 1000, seed = 1) {
  out <- lapply(k_values, function(k) {
    r <- run_kmeans_ht(
      df = df,
      scores = scores,
      pcs_use = pcs_use,
      k = k,
      region_col = region_col,
      nstart = nstart,
      max_iter = max_iter,
      seed = seed
    )

    data.frame(
      pcs = paste(pcs_use, collapse = ","),
      d = length(pcs_use),
      k = k,
      HTA = r$HTA,
      err_msg = ifelse(is.null(r$err_msg), NA_character_, r$err_msg),
      iter = r$iter,
      tot_withinss = r$tot_withinss,
      size_min = min(r$size),
      size_max = max(r$size)
    )
  })

  do.call(rbind, out)
}

# Sweep both PC sets and k values.
sweep_pcsets_k_ht <- function(df, scores, pc_sets, k_values,
                              region_col = "region_2x2",
                              nstart = 100, max_iter = 1000, seed = 1) {
  all_res <- do.call(
    rbind,
    lapply(names(pc_sets), function(set_name) {
      pcs_use <- pc_sets[[set_name]]

      tmp <- sweep_k_ht(
        df = df,
        scores = scores,
        pcs_use = pcs_use,
        k_values = k_values,
        region_col = region_col,
        nstart = nstart,
        max_iter = max_iter,
        seed = seed
      )

      tmp$pc_set <- set_name
      tmp
    })
  )

  all_res <- all_res[, c("pc_set", setdiff(names(all_res), "pc_set"))]
  all_res[order(all_res$pc_set, all_res$k), , drop = FALSE]
}
