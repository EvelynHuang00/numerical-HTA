# Numerical HTA for Tabular Health Data

This project adapts the HeTerogeneity Average (HTA), originally designed around local spatial regions, to tabular health data. The main analysis defines local regions in Age × BMI covariate space and measures within-region heterogeneity using joint demographic and behavioral trait combinations.

The repository is organized around the analyses retained in the final report, while earlier methodological experiments are kept separately under `exploratory/`.

## Main analyses

### 1. Categorical Numerical HTA

Age and BMI define quantile-based regions. Joint combinations of Gender, Ethnicity, EducationLevel, and Smoking define patient traits. Local normalized-entropy scores are aggregated by region size to obtain a global HTA score.

The analysis also evaluates sensitivity to grid resolution, compares the Alzheimer cohort with a pseudo-homogeneous dataset, and benchmarks Numerical HTA against Alpha Simpson diversity.

### 2. PCA tertile sensitivity analysis

For a numerically dominated health cohort, Age and BMI remain the region-defining variables and are excluded from PCA. PCA is fit to the remaining standardized numeric variables, and PC1 is divided into Low/Mid/High tertiles to define an alternative trait representation.

This branch is presented as a sensitivity analysis because the tertile construction can itself create balanced trait categories and therefore high apparent heterogeneity.

## Repository structure

```text
.
├── README.md
├── STRUCTURE.md
├── R/
│   ├── hta_core.R
│   ├── pca_helpers.R
│   └── clustering_helpers.R
├── analysis/
│   ├── 01_categorical_hta.R
│   └── 02_pca_tertile_sensitivity.R
├── exploratory/
│   ├── 01_permutation_and_variance.R
│   ├── 02_pca_alternative_traits.R
│   └── 03_clustering_sensitivity.R
├── figures/
├── data/
│   └── dirty_v3_path.csv
│   └── psedo_homogeneous_dataset.csv
│   └── alzheimers_disease_data.csv
├── report/
│   └── HTA_final_report.pdf
└── archive/                 # local provenance; ignored by Git by default
```

## Running the analysis

Run scripts from the repository root after placing the required CSV files in `data/`:

```r
source("analysis/01_categorical_hta.R")
source("analysis/02_pca_tertile_sensitivity.R")
```

The first script saves the report-aligned figures to `figures/`.

See `STRUCTURE.md` for a detailed mapping from the original Rmd notebooks to the cleaned project structure and for notes on the few non-methodological fixes made during reorganization.
