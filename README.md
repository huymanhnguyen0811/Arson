# Arson & Wildfire Accelerant Classification — Computational Workflow

> **TL;DR**: This repository contains an end‑to‑end R workflow to preprocess GC×GC‑MS data, engineer a robust set of features, test imputation + normalization strategies, and train/benchmark classifiers to distinguish sources of accelerants used in wildfire arson.

---

## 📌 About this project

This repo accompanies the manuscript:

> **“An open-access computational fingerprinting workflow for source classifications of neat gasoline using GC×GC-TOFMS and Machine Learning”**

---

## 🧰 System requirements

- **OS**: Developed on Windows 11 (11th Gen Intel® Core™ i7‑11800H, 16 GB RAM)
- **R**: 4.3.1 (or newer recommended)
- **RStudio**: 2023.06.0+421 “Mountain Hydrangea” (or newer)

> Other platforms should work, but you may need to adapt installation steps (e.g., for `lightgbm`).

---

## 📦 R packages

Core packages used throughout the workflow:

- **Data I/O & wrangling**: `readxl`, `writexl`, `tidyverse` (`dplyr`, `tidyr`, `purrr`), `data.table`
- **Visualization**: `ggplot2`, `grid`, `gridExtra`, `ggpubr`, `ggsignif`, `viridis`, `ggforce`
- **Multivariate / PCA**: `FactoMineR`, `factoextra`, `vegan`
- **Modeling**: `caret`, `randomForest`, `xgboost`, `lightgbm`, `pROC`
- **Imputation**: `missForest`, `VIM` (kNN), `softImpute` (if used)
- **Heatmaps**: `pheatmap`

Install (first time only):

```r
pkgs <- c(
  "readxl","writexl","tidyverse","data.table","ggplot2","grid","gridExtra",
  "ggpubr","ggsignif","viridis","ggforce","FactoMineR","factoextra","vegan",
  "caret","randomForest","xgboost","pROC","missForest","VIM","pheatmap"
  # optional: "lightgbm","softImpute"
)
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
```

---

## 🗂️ Repository layout (recommended)

```
.
├─ data/                       # input spreadsheets (not committed)
├─ scripts/                    # optional: utility R scripts
├─ figures/                    # figures saved by the workflow
├─ results/                    # exported tables, metrics
├─ Arson workflow instruction.Rmd
└─ README.md
```

> **Data**: The workflow expects an input file named **`Gasolines_BOP_110424.xlsx`** and a target/compound list such as **`Shortened ILR Compound List PF001A 07-06-2024.xlsx`**. Place these in `data/` (or update paths accordingly).

---

## 🚀 Getting started

1. **Clone and open** this repository in RStudio.
2. **Set the project root**:
   ```r
   # in the Rmd, prefer a project root like this:
   knitr::opts_knit$set(root.dir = getwd())
   ```
3. **Put input files** into `./data/`.
4. **Knit** the analysis:
   ```r
   rmarkdown::render("Arson workflow instruction.Rmd",
                     output_format = rmarkdown::md_document(variant = "gfm"))
   ```
   Artifacts (figures, tables) are saved to `figures/` and `results/` as configured in the Rmd.

---

## 🔄 End‑to‑end workflow (high‑level)

The Rmd is organized into numbered steps. Below is a faithful summary to help readers follow along in plain Markdown.

### STEP 1.1 — Data import
- Read **`Gasolines_BOP_110424.xlsx`** (multiple sheets are bound into one data frame with a `Sample_name` column).
- Map codes to **octane ratings** (`Gas_87`,`Gas_89`,`Gas_91`,`Gas_94`) and **seasons** (`blue`, `purple`, `orange`).  
- Map station IDs (e.g., `F001`→`Station_1`, …, `F010`→`Station_10`).
- Rename GC×GC‑MS columns to consistent names: `RT1`, `RT2`, `Ion1`, `Ion2`.

### STEP 1.2 — Remove solvent, column bleed, BTEX, and `MF = 0`
- Filter out common interferents/bleed (e.g., **Carbon disulfide**, cyclic siloxanes, etc.).
- Optional filters for BTEX can be toggled as needed.
- Produce **QA histograms** for peak distributions prior to normalization.

#### STEP 1.2B — Coverage vs. intensity threshold
- Sweep a peak‑area threshold; for each threshold, compute:
  - % total peak area retained
  - number of peaks retained
- Plot coverage curves to justify the selected cutoff.

### STEP 2A — Grouping by retention & ions
- Group features by (`RT1`,`RT2`,`Ion1`,`Ion2`) windows; use benchmark compound(s) (e.g., **toluene**) to check RT stability across batches.

### STEP 2B — Data compression
- Collapse features that share `Ion1` (m/z) ranges and matched `RT1/RT2` windows to reduce redundancy and stabilize alignment.

### STEP 3 — Feature elimination / prioritization
- Remove uninformative or unstable features based on grouping logic, QC consistency, and domain knowledge.

### STEP 4 — Missingness filtering
- Drop features with **>90% missing** overall.

### STEP 4 — Imputation & normalization grid
Imputation methods implemented include (examples shown in code):
- Constant fill (e.g., **0.001**)
- Mean / Median
- **kNN** (`VIM::kNN`)
- **Random Forest** (`missForest`)
- (Optionally) **SoftImpute/SVD**

Normalization choices include:
- **None**
- **Log** transform
- **Min–Max** scaling
- **Z‑score** scaling

Each (imputation × normalization) combination is evaluated downstream.

### Evaluation & model training
- **Cross‑validation** / OOB fallback when data are sparse (adaptive inner CV size).
- **Random Forest** hyper‑parameter search over `mtry` and `ntree` with OOB or `caret` CV depending on sample size.
- Metrics recorded per combo:
  - **Accuracy**, **AUC**, **F1‑weighted**, **Kappa**, **MCC (multiclass)**
  - Unsupervised **cluster resolution**
  - **Correlation score** vs. raw distribution
  - **K–S p‑value** (distributional safeguard)
- Aggregate results table is ranked primarily by **Accuracy** and **AUC**, then by structure‑preserving metrics (cluster resolution, correlation), and statistical diagnostics.

### Statistical testing
- Pairwise significance testing per feature across classes:
  - Shapiro–Wilk normality check → **t‑test** (parametric) or **Wilcoxon** (non‑parametric), paired where appropriate.
  - Multiple‑comparison control applied to p‑values (see Rmd for details).

### Visualization (selected)
- Coverage curves, RT histograms (per benchmark), PCA/ordination plots
- Heatmaps (`pheatmap`), ROC curves (`pROC`), confusion matrices
- (Optional) feature importance / explainability visualizations

---

## 📁 Inputs

- `data/Gasolines_BOP_110424.xlsx` — main input (multiple sheets; one per sample group)
- `data/Shortened ILR Compound List PF001A 07-06-2024.xlsx` — compound/target list
- (If you use different filenames or locations, update the paths in the Rmd or add a config block.)

---

## 📤 Outputs

Typical outputs include:
- `results/` — tables summarizing imputation × normalization performance, rankings, and statistical tests
- `figures/` — coverage plots, RT histograms, PCA, ROC, confusion matrices, and heatmaps

> Exact filenames and locations are controlled within the Rmd; adjust to your preferred structure.

---

## 🔁 Reproducibility tips

- **Root directory**: avoid machine‑specific paths; rely on project root (`getwd()`), or use `here::here()`.
- **Seeds**: set `set.seed(123)` (done throughout) for reproducible splits/imputation.
- **Package versions**: record with `sessionInfo()` when publishing results.

---

## 📝 Citation

If you build on this work, please cite the manuscript (title below; year/journal TBD):

> Nguyen et al. **“The Use of Computational Fingerprinting Techniques to Distinguish Sources of Accelerants Used in Wildfire Arson.”**

---

## 📄 License

Add a license (e.g., MIT, Apache‑2.0) to clarify reuse. If data are proprietary or sensitive, **do not commit raw spreadsheets**; share de‑identified or synthetic examples instead.

---

## 🙏 Acknowledgments

Thanks to collaborators and co‑authors who contributed sample collection, analytical chemistry, and domain expertise.
