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

## 🔄 End‑to‑end workflow

The Rmd is organized into numbered steps. Below is a faithful summary to help readers follow along in plain Markdown.

### STEP 1.1 — Data import
- Read **`Gasolines_BOP_110424.xlsx`** (multiple sheets are bound into one data frame with a `Sample_name` column).
- Map codes to **octane ratings** (`Gas_87`,`Gas_89`,`Gas_91`,`Gas_94`) and **seasons** (`blue`, `purple`, `orange`).  
- Map station IDs (e.g., `F001`→`Station_1`, …, `F010`→`Station_10`).
- Rename GC×GC‑MS columns to consistent names: `RT1`, `RT2`, `Ion1`, `Ion2`.

```r
target_comp <- read_xlsx(path = "Shortened ILR Compound List PF001A 07-06-2024.xlsx")
# ASTM <- read_xlsx(path = "ILR Compound List 05-15-2024_Without DieselASTM.xlsx")
file_path <- "Gasolines_BOP_110424.xlsx"

dfs <- excel_sheets(file_path) %>%
  set_names() %>%
  map(~ read_excel(file_path, sheet = .x) %>% mutate(Sample_name = .x))

df_step1.1 <- bind_rows(dfs) %>%
  dplyr::select(-c("RMF", "Area %")) %>%
  mutate(Octane_rating = ifelse(str_detect(Sample_name, "A"), "Gas_87", 
                                ifelse(str_detect(Sample_name, "B"), "Gas_89",
                                       ifelse(str_detect(Sample_name, "C"), "Gas_91", "Gas_94")))) %>%
  mutate(sampling_season = ifelse(str_detect(Sample_name, "b"), "blue",
                                  ifelse(str_detect(Sample_name, "p"), "purple", "orange"))) %>%
  mutate(gas_station = ifelse(str_detect(Sample_name, "F001"), "Station_1",
                              ifelse(str_detect(Sample_name, "F002"), "Station_2",
                                     ifelse(str_detect(Sample_name, "F003"), "Station_3",
                                            ifelse(str_detect(Sample_name, "F004"), "Station_4",
                                                   ifelse(str_detect(Sample_name, "F005"), "Station_5",
                                                          ifelse(str_detect(Sample_name, "F006"), "Station_6",
                                                                 ifelse(str_detect(Sample_name, "F007"), "Station_7",
```

### STEP 1.2 — Remove solvent, column bleed, BTEX, and `MF = 0`
- Filter out common interferents/bleed (e.g., **Carbon disulfide**, cyclic siloxanes, etc.).
- Optional filters for BTEX can be toggled as needed.
- Produce **QA histograms** for peak distributions prior to normalization.

#### Quality Assurace — Coverage vs. intensity threshold
- Sweep a peak‑area threshold; for each threshold, compute:
  - % total peak area retained
  - number of peaks retained
- Plot coverage curves to justify the selected cutoff.

```r
# Calculate the mean/median percentage of all sample in df_step1.2
coverage_list <- c()
for (threshold in c(seq(from = 0, to = 1000000, by = 50000))) {
  df_filter_area <- df_step1.2 %>%
    dplyr::filter(Area > threshold) %>%
    dplyr::group_by(Sample_name) %>%
    dplyr::summarise(across(Area, base::sum))
  
  coverage <- c()
  for (sample in df_filter_area$Sample_name) {
    coverage <- c(coverage, df_filter_area[which(df_filter_area$Sample_name == sample),]$Area*100/sum(df_step1.2[which(df_step1.2$Sample_name == sample),]$Area))
  }
  coverage_list <- c(coverage_list, mean(coverage))
}

num_comp_list <- c()

# Precompute sample counts for all thresholds using vectorized operations
for (threshold in seq(from = 0, to = 1000000, by = 50000)) {
  
  # Filter the data based on threshold
  df_filter_area <- df_step1.2 %>%
    dplyr::filter(Area > threshold) %>%
    dplyr::group_by(Sample_name) %>%
    dplyr::summarize(num_comp = n())   # Count occurrences per Sample_name
  
  # Compute the mean of num_comp and store it in num_comp_list
  num_comp_list <- c(num_comp_list, mean(df_filter_area$num_comp))
}

df <- data.frame(list(thres = c(seq(from = 0, to = 1000000, by = 50000)), remain = num_comp_list))
plot <- ggplot(data = df,
               aes(x = thres, y = remain)) +
  geom_col(fill = "skyblue") +
  geom_text(aes(label = round(remain, digits = 0)), color = "black", 
            angle = 90, hjust = 1.5, size = 9) +
  geom_text(aes(label = paste0(round(coverage_list, digits =1),"%"), 
                y = remain + 5),  # Adjust 'y = remain + 5' for positioning above bars
            color = "black", size = 9, vjust = 0) +     # Adjust vjust for fine-tuning
  scale_x_continuous(breaks = seq(from = 0, to = 1000000, by = 50000),
                     # remove space between plotted data and xy-axes
                     expand = c(0,0)) +
  theme_classic(base_size = 25) +
  theme(axis.text.x = element_text(size = 25, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 25),
        strip.text = element_text(size = 30, face = "bold"),   # Make facet label text bold
        strip.background = element_blank()) +
  labs(x = "Stepwise thresholds of peak area for removal of compounds with low signal to noise", 
       y = "Number of peak remains after applying thresholds of peak area") 
plot
```

```r
filter_list <- c("^Carbon disulfide$", 
                "Cyclotrisiloxane..hexamethyl",
                "Cyclotetrasiloxane..octamethyl"
                # "^Benzene$",
                # "^Toluene$",
                # "^Ethylbenzene$",
                # "Xylene"
                )

# ^Carbon disulfide$ 75.890 - 75.959, 77.881 - 77.948

df_step1.2 <- copy(df_step1.1) # %>%
  # filter(MF > 0)

for (filter_comp in filter_list) {
  df_step1.2 <- df_step1.2 %>%
      filter(!grepl(filter_comp, Compound))
}

df_step2 <- df_step1.2 %>%
  filter(Area > 300000) %>%
  arrange(RT1, RT2)
```

### STEP 2A — Grouping by retention & ions
- Group features by (`RT1`,`RT2`,`Ion1`,`Ion2`) windows; use benchmark compound(s) (e.g., **toluene**) to check RT stability across batches.

#### Quality assurance 1. Select optimal alignment window by confirming targeted compounds at each step of Data processing
```r
# For each target compounds,  Apply different alignment windows and record number of samples that have the =target compounds
combi <- tidyr::crossing(
  # RT1
  c(0.1 ,0.15, 0.2, 0.25, 0.3, 0.35, 0.4), 
  # RT2
  c(0.12,0.13,0.14,0.15, 0.16, 0.17, 0.18, 0.19, 0.2))

df <- data.frame(RT1_window=integer(), RT2_window=integer(), target=character(), proportion=integer())
  
for (i in 1:nrow(combi)) {
  for (j in 1:nrow(target_comp)) {
    # Catch all peaks in dataframe that falling into the window with target compounds as center of the window
    idx1 <- which((df_step2$Ion1 <= (target_comp[j,]$Ion1 + 0.1) & df_step2$Ion1 >= (target_comp[j,]$Ion1 - 0.1) & 
                   df_step2$Ion2 <= (target_comp[j,]$Ion2 + 0.1) & df_step2$Ion2 >= (target_comp[j,]$Ion2 - 0.1)) |
                    
                   (df_step2$Ion1 <= (target_comp[j,]$Ion2 + 0.1) & df_step2$Ion1 >= (target_comp[j,]$Ion2 - 0.1) &
                   df_step2$Ion2 <= (target_comp[j,]$Ion1 + 0.1) & df_step2$Ion2 >= (target_comp[j,]$Ion1 - 0.1)))
    
    temp <- df_step2[idx1,]
    minrt1 <- max(temp$RT1)
    maxrt1 <- min(temp$RT1)
    minrt2 <- max(temp$RT2)
    maxrt2 <- max(temp$RT2)
    
    
    idx2 <- which(temp$RT1 <= (target_comp[j,]$RT1 + as.numeric(combi[i,][1])) & 
                    temp$RT1 >= (target_comp[j,]$RT1 - as.numeric(combi[i,][1])) &
                    temp$RT2 <= (target_comp[j,]$RT2 + as.numeric(combi[i,][2])) & 
                    temp$RT2 >= (target_comp[j,]$RT2 - as.numeric(combi[i,][2])))
    
    df[nrow(df) + 1,] <- c(as.numeric(combi[i,][1]),
                           as.numeric(combi[i,][2]),
                           target_comp[j,]$Compound,
                           # paste0(length(unique(temp[idx2,]$Sample_name)), "/", 71))
                           as.numeric(100*(length(unique(temp[idx2,]$Sample_name)) / length(unique(temp$Sample_name)))))
  }
}

df$proportion <- as.numeric(df$proportion)

summary_df <- df %>% pivot_wider(names_from = target, values_from = proportion)

# Which window have the highest proportion of samples that have the target compounds
summary_df[, 3:ncol(summary_df)] <- lapply(summary_df[, 3:ncol(summary_df)], as.numeric)
max_row <- summary_df[which.max(rowSums(summary_df[, 3:ncol(summary_df)])), ]
print(max_row)

# Adding the number of samples where the compounds were found with matching Ion1 and Ion2 to colnames of each compound
i <- 3
for (j in 1:nrow(target_comp)) {
  idx1 <- which((df_step2$Ion1 <= (target_comp[j,]$Ion1 + 0.1) & df_step2$Ion1 >= (target_comp[j,]$Ion1 - 0.4) & 
                   df_step2$Ion2 <= (target_comp[j,]$Ion2 + 0.1) & df_step2$Ion2 >= (target_comp[j,]$Ion2 - 0.1)) |
                  
                  (df_step2$Ion1 <= (target_comp[j,]$Ion2 + 0.1) & df_step2$Ion1 >= (target_comp[j,]$Ion2 - 0.1) &
                     df_step2$Ion2 <= (target_comp[j,]$Ion1 + 0.1) & df_step2$Ion2 >= (target_comp[j,]$Ion1 - 0.1)))
  
  temp <- df_step2[idx1,]
  idx2 <- which(temp$RT1 <= (target_comp[j,]$RT1 + as.numeric(max_row$RT1_window)) & 
                  temp$RT1 >= (target_comp[j,]$RT1 - as.numeric(max_row$RT1_window)) &
                  temp$RT2 <= (target_comp[j,]$RT2 + as.numeric(max_row$RT2_window)) & 
                  temp$RT2 >= (target_comp[j,]$RT2 - as.numeric(max_row$RT2_window)))
  
  colnames(summary_df)[i] <- paste0(target_comp[j,]$Compound, " (n = ", length(unique(temp$Sample_name)), ")")
  i <- i + 1 
}

View(summary_df)
```
#### Quality assurance 2. Examine the distribution of RT of target compounds after alignment with pF001A
```r
# Get the column names from the existing data frame
column_names <- c( "Target compound", "RT", "Retention time of pF001A", "Min.", "1st Qu.", "Median",  "Mean", "3rd Qu.", "Max.", "Max. RT - Min. RT")
# Create an empty data frame with the same column names
summary_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(summary_df) <- column_names

summary_target_compounds <- list()
for (j in 1:nrow(target_comp)) {
    # First for each target compound, matching Major and minor ion 
    idx1 <- which((df_step2$Ion1 <= (target_comp[j,]$Ion1 + 0.2) & df_step2$Ion1 >= (target_comp[j,]$Ion1 - 0.2) & 
                   df_step2$Ion2 <= (target_comp[j,]$Ion2 + 0.2) & df_step2$Ion2 >= (target_comp[j,]$Ion2 - 0.2)) |
                    
                   (df_step2$Ion1 <= (target_comp[j,]$Ion2 + 0.2) & df_step2$Ion1 >= (target_comp[j,]$Ion2 - 0.2) &
                   df_step2$Ion2 <= (target_comp[j,]$Ion1 + 0.2) & df_step2$Ion2 >= (target_comp[j,]$Ion1 - 0.2)))
    
    temp <- df_step2[idx1,]
    
    # Then, for each target compound, matching within RT +- 0.1 wrt the pF001A retention time. 
    idx2 <- which(temp$RT1 <= (target_comp[j,]$RT1 + 0.1) & 
                  temp$RT1 >= (target_comp[j,]$RT1 - 0.1) &
                  temp$RT2 <= (target_comp[j,]$RT2 + 0.1) & 
                  temp$RT2 >= (target_comp[j,]$RT2 - 0.1))
    
    
    
    if (nrow(temp[idx2, ]) == 0) {
      summary_df[nrow(summary_df) + 1,] <- c(paste0(target_comp[j,]$Compound, " was not found with matching ions and within RT1/RT2 windows of 0.1 of pF001A"), NA, NA, NA, NA, NA, NA, NA, NA, NA)
    } else {
      
      # Make descriptive stats summary of RT1 of all peaks that was aligned to target compound in pF001A
      summary_df[nrow(summary_df) + 1,] <- c(target_comp[j,]$Compound,
                                             "RT1",
                                             as.numeric(target_comp[j,]$RT1),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT1))[[1]]),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT1))[[2]]),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT1))[[3]]),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT1))[[4]]),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT1))[[5]]),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT1))[[6]]),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT1))[[6]] - summary(as.numeric(temp[idx2,]$RT1))[[1]]))
      
      # Make descriptive stats summary of RT2 of all peaks that was aligned to target compound in pF001A
      summary_df[nrow(summary_df) + 1,] <- c(target_comp[j,]$Compound,
                                             "RT2",
                                             as.numeric(target_comp[j,]$RT2),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT2))[[1]]),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT2))[[2]]),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT2))[[3]]),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT2))[[4]]),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT2))[[5]]),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT2))[[6]]),
                                             as.numeric(summary(as.numeric(temp[idx2,]$RT2))[[6]] - summary(as.numeric(temp[idx2,]$RT2))[[1]]))
    }
}

View(summary_df)
```

### STEP 2B — Data compression
- Collapse features that share `Ion1` (m/z) ranges and matched `RT1/RT2` windows to reduce redundancy and stabilize alignment.

```r
# Define tolerances
tolerances <- list(RT1 = 0.1, RT2 = 0.1, Ion1 = 0.5, Ion2 = 0.5)

# Use pF001A as a base
df_all <- df_step2 %>% 
  filter(Sample_name %in% "pF001A") %>% 
  filter(`Signal to Noise` > 10)

df_all$Feature <- 1:nrow(df_all)

# Loop through the samples apart from pF001A
for (sample in setdiff(unique(df_step2$Sample_name), c("pF001A", "bF001A", "bF007B"))) {
  # print(sample)
  df <- df_step2 %>% 
    filter(Sample_name %in% sample) %>% 
    filter(`Signal to Noise` > 10)
  df$Feature <- NA
  
  # Go through each row
  for (i in 1:nrow(df)) {
    row <- df[i, ]
    mask <- (
      abs(df_all$RT1 - row$RT1) <= tolerances$RT1 &
      abs(df_all$RT2 - row$RT2) <= tolerances$RT2 & 
      abs(df_all$Ion1 - row$Ion1) <= tolerances$Ion1 &
      abs(df_all$Ion2 - row$Ion2) <= tolerances$Ion2
    )
    
    idx <- which(mask)
    # If there is a match between a peak and the existing peak list, then assign the same Feature number to that peak 
    if (any(mask)) {
      row$Feature <- unique(df_all[idx, ]$Feature)[1]
    } else { # If not a match, then create new identity for the new Feature 
      row$Feature <- max(df_all$Feature) + 1
    }
    
    # adding the peak that have matchs in df_all
    df_all <- bind_rows(df_all, row)
  }
}

df_all <- df_all %>%
  # remove all peak with RT1 > 36 (which does not belong Gasoline)
  filter(RT1 < 30)

# Create metadata for next data analysis
metadata <- df_all %>%
  dplyr::select(
    Feature,
    Sample_name,
    Area) %>%
  group_by(Feature,
           Sample_name) %>%
  dplyr::summarise(across(Area, base::mean)) %>% # If feature appear in 
  tidyr::pivot_wider(names_from = Sample_name,
                     values_from = Area) 

#### Adding avg RT and Ions to the master df for Gwen
mean_rt1 <- c()
mean_rt2 <- c()
mean_ion1 <- c()
mean_ion2 <- c()

for (i in 1:nrow(metadata)) {
  mean_rt1 <- c(mean_rt1, mean(df_all[which(df_all$Feature %in% metadata[i, ]$Feature),]$RT1))
  mean_rt2 <- c(mean_rt2, mean(df_all[which(df_all$Feature  %in% metadata[i, ]$Feature),]$RT2))
  mean_ion1 <- c(mean_ion1, max(df_all[which(df_all$Feature  %in% metadata[i, ]$Feature),]$Ion1)) 
  mean_ion2 <- c(mean_ion2, max(df_all[which(df_all$Feature  %in% metadata[i, ]$Feature),]$Ion2))
}

metadata$RT1 <- mean_rt1
metadata$RT2 <- mean_rt2
metadata$Ion1 <- mean_ion1
metadata$Ion2 <- mean_ion2

metadata <- metadata %>% relocate(RT1, RT2, Ion1, Ion2, .after = 1)

core_metadata <- metadata %>% 
  # relocate(Chemical_group, .after = 1) %>% 
  arrange(RT1, RT2)
```

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

> Nguyen et al. **“An open-access computational fingerprinting workflow for source classifications of neat gasoline using GC×GC-TOFMS and Machine Learning”**

---

## 📄 License

Add a license (e.g., MIT, Apache‑2.0) to clarify reuse. If data are proprietary or sensitive, **do not commit raw spreadsheets**; share de‑identified or synthetic examples instead.

---

## 🙏 Acknowledgments

Thanks to collaborators and co‑authors who contributed sample collection, analytical chemistry, and domain expertise.
