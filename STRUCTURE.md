# Code structure

This repository reorganizes three exploratory R Markdown notebooks into a small set of reusable functions, two report-aligned analysis scripts, and separate exploratory scripts.

## Public-facing analysis

### `analysis/01_categorical_hta.R`
Contains the analysis that forms the main result in the final report:

1. Numerical HTA on the Alzheimer cohort.
2. Sensitivity to Age × BMI grid resolution.
3. Within-cell trait composition.
4. Comparison with the pseudo-homogeneous dataset.
5. Alpha Simpson benchmark.
6. Generation of Figures A1–A4.

### `analysis/02_pca_tertile_sensitivity.R`
Contains the PCA extension retained in the final report:

1. Age and BMI define a 2×2 regional grid.
2. Age and BMI are excluded from PCA.
3. Remaining numeric variables are standardized.
4. PC1 is split into Low/Mid/High tertiles.
5. PC1 tertiles define the trait labels used for HTA.
6. The same analysis is repeated on the pseudo-homogeneous numerical cohort.

## Reusable functions

### `R/hta_core.R`
Functions from `numerical_HTA.Rmd`:

- `hti()`
- `qbreaks()`
- `compute_hta()`
- `make_comp_df()`
- `make_cells_quantile()`
- `gs_local()`
- `compute_gs()`

The permutation branches inside `compute_hta()` are retained unchanged even though the main report-facing script runs HTA descriptively with `perm_B = 0`.

### `R/pca_helpers.R`
Functions extracted from the PCA notebook:

- construction of the 2×2 Age × BMI regions;
- PCA on standardized numeric variables after excluding Age and BMI;
- the fixed-trait-level HTA calculation used in the PCA notebook;
- the earlier PC-defined region/trait HTA function;
- generation of the pseudo-homogeneous numerical cohort.

### `R/clustering_helpers.R`
K-means + HTA sweep functions from the clustering notebook.

## Exploratory analyses

### `exploratory/01_permutation_and_variance.R`
Early single-variable HTI checks, permutation diagnostics, tau-squared/ICC, and mixed-model experiments from `numerical_HTA.Rmd`.

### `exploratory/02_pca_alternative_traits.R`
PCA formulations that were tested but not retained as the final PCA result, including PC-defined regions/traits, extreme-bin PC traits, and the experimental `HTA_ML` adjustment.

### `exploratory/03_clustering_sensitivity.R`
K-means trait definitions, PC-set and cluster-count sweeps, diagnostics for unexpectedly high HTA in the homogeneous cohort, random-region checks, and a permutation null.

## Logic-preservation notes

The goal of this reorganization is structural rather than methodological. Metric definitions, binning choices, PCA standardization, synthetic-data generation parameters, k-means settings, and permutation settings were preserved from the original notebooks.

A few non-methodological fixes were necessary so the reorganized scripts are internally consistent:

- Absolute local file paths were replaced with repository-relative `data/...` paths.
- The original grid-resolution convenience loop omitted `trait_cols`; the cleaned script passes the same categorical trait definition used everywhere else in that analysis.
- The clustering notebook created a randomized region variable on one object and then tried to use it on another. The cleaned exploratory script attaches the same randomized labels to the object actually passed to the sweep.
- Repeated definitions and duplicate plots were consolidated into helper files.

Two HTI normalization implementations from the original work are intentionally not collapsed into one. The main categorical pipeline uses the `hti()` convention in `numerical_HTA.Rmd`, while the later PCA/clustering notebooks use a fixed number of global trait levels in their `compute_ht()` calculation. Keeping them separate avoids silently changing the original experimental logic.

## Archive

`archive/` contains the three original Rmd files for provenance. It is excluded by `.gitignore` so the public GitHub repository can stay concise. Remove the `archive/` ignore rule only if you explicitly want the raw research notebooks committed.
