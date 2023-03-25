# Instruction to Workflow - Arson

```{r, warning = FALSE, message = FALSE, echo = FALSE}
library(DiagrammeR)
library(DiagrammeRsvg)
library(ggplot2)
library(readxl)
library(ggpubr)
library(gtable)
library(gridExtra)
library(viridis)
library(wesanderson)
library(tidyverse)
library(lubridate)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(data.table)
library(purrr)
library(missMDA)
library(tsne)
library(Rtsne)
library(rgl)
library(collapse)
library(plotly)
library(umap)
library(rsvg)
library(knitr)
library(kableExtra)
library(grid)
```

# Arson Data Workflow

```{r, warning = FALSE, message = FALSE, echo = FALSE, fig.height=10, fig.width=18}
knitr::opts_chunk$set(dpi = 1000)
DiagrammeR::grViz("digraph my_flowchart  {
      graph[splines = ortho] # to get 90 degree angles and straight lines.
      node [layout = dot, fontname = Helvetica, shape = box, fixedsize = false, width = 4, height = 1]

      node1[label = <<FONT COLOR='blue' POINT-SIZE='25'><b>Data import </b></FONT>>];
      node2[label = <<b><font color='blue' POINT-SIZE='25'>Filtering out column bleed / solvent / BTEX and minimum area limit observations (associated with GCxGC system)     </font></b>>]
      node3[label = <<b><font color='blue' POINT-SIZE='25'>Sorting data in ascending order of RT1 and RT2      </font></b>>]
      node4[label = <<b><font color='blue' POINT-SIZE='25'>Grouping identical compounds <br/>based on RT1, RT2, Ion1 threshold      </font></b>>]
      node5[label = <<b><font color='blue' POINT-SIZE='25'>Identify shared and unique compound    <br/>groups across all 31 IL samples        </font></b>>]
      node6[label = <<b><font color='blue' POINT-SIZE='25'>Data Normalization (TSN as 1st layer)   </font></b>>]
      node6a[label = <<b><font color='red' POINT-SIZE='25'>Total Sum Normalization   </font></b>>]
      node7[label = <<b><font color='blue' POINT-SIZE='25'>Statistical test to identify potential markers      </font></b>>]
      node8[label = <<b><font color='blue' POINT-SIZE='25'>Test potential markers with clustering analysis      </font></b>>]
      node8a[label = <<b><font color='red' POINT-SIZE='25'>Principle Component Analysis    </font></b>>]
      node8b[label = <<b><font color='red' POINT-SIZE='25'>Hierarchical Clustering on   <br/>Principle Components     </font></b>>]
      node8c[label = <<b><font color='red' POINT-SIZE='25'>T-distributed stochastic <br/>neighbor embedding (t-SNE)     </font></b>>]
      node8d[label = <<b><font color='red' POINT-SIZE='25'>Uniform Manifold Approximation    <br/>and Projection (UMAP)</font></b>>]
      node9[label = <<b><font color='blue' POINT-SIZE='25'>Regression Classification   </font></b>>]
      node10[label = <<b><font color='blue' POINT-SIZE='25'>Quality Control (cross-validation with ASTM) </font></b>>]

      blank1[label = '', width = 0.01, height = 0.01]
      blank1a[label = '', width = 0.01, height = 0.01]
      blank1b[label = '', width = 0.01, height = 0.01]
      blank1c[label = '', width = 0.01, height = 0.01]
      m1 [label = <<b><font color='red' POINT-SIZE='25'>Pair-wise comparison</font></b>>]
      m2 [label = <<b><font color='red' POINT-SIZE='25'>Multiple comparison</font></b>>]

      node7 -> blank1[dir = none];
      blank1 -> blank1a[dir = none, minlen = 12];
      {rank = same; blank1 blank1a}
      blank1a -> blank1b[dir = none, maxlen = 1];
      blank1a -> blank1c[dir = none, maxlen = 1];
      {rank = same; blank1b blank1c}
      blank1b -> m1[maxlen = 1];
      {rank = same; blank1b m1}
      blank1c -> m2[maxlen = 1];
      {rank = same; blank1c m2}
      blank1 -> node8;

      m1a [label = <<b><font color='darkgreen' POINT-SIZE='25'>Kolmogorov-Smirnov test    </font></b>>]
      m1b [label = <<b><font color='darkgreen' POINT-SIZE='25'>Mann-Whitney <br/>Rank Sum test</font></b>>]
      m1c [label = <<b><font color='darkgreen' POINT-SIZE='25'>Student t-Test</font></b>>]

      blank2[label = '', width = 0.01, height = 0.01]
      m1 -> blank2[dir = none, maxlen = 1];
      blank2 -> m1a[maxlen = 1];
      blank2 -> m1b[maxlen = 1];
      blank2 -> m1c[maxlen = 1];
      {rank = same; m1a m1b m1c}

      m2a [label = <<b><font color='darkgreen' POINT-SIZE='25'>ANOVA</font></b>>]
      m2b [label = <<b><font color='darkgreen' POINT-SIZE='25'>ANCOVA</font></b>>]

      blank3[label = '', width = 0.01, height = 0.01]
      m2 -> blank3[dir = none, maxlen = 1];
      blank3 -> m2a[maxlen = 1];
      blank3 -> m2b[maxlen = 1];
      {rank = same; m2a m2b}

      blank4[label = '', width = 0.01, height = 0.01]
      node6 -> blank4[dir = none];
      blank4 -> node6a[minlen = 6];
      {rank = same; blank4 node6a}
      blank4 -> node7


      node8 -> node8a;

      node8a -> node8b;

      node8 -> node8c;

      node8 -> node8d;


     # create undirected edge from source to dummy node
      edge [dir=normal]
      node1 -> node2 -> node3 -> node4 -> node5 -> node10 -> node6;
      node8 -> node9[minlen = 6];
    }")

```

<br>

# Summary of findings

```{r, warning = FALSE, message = FALSE, echo = FALSE}
# Functions -------------------------------------------------------------------------------------------------------
# Filtering matched compound names
filtering <- function(df, filter_list) {
  clean_data <-  copy(df)
  for (ele in filter_list) {
    clean_data <- clean_data %>%
      filter(!grepl(ele, Compound))
  }
  return(clean_data)
}

# Filtering limit of observations
limit_obser <- function(df_list, file_list, cap) {
  df_list_filter_area <- list()
  df_list_removed_area <- list()
  for (i in 1:length(df_list)) {
    df_list_filter_area[[i]] <- df_list[[i]] %>%
      filter(., Area > cap) %>%
      mutate(sample_name = file_list[[i]]) %>%
      mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                                ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                       ifelse(str_detect(sample_name, "D"), "Diesel", "Gas"))))
    
    df_list_removed_area[[i]] <- df_list[[i]] %>%
      filter(., Area <= cap) %>%
      mutate(sample_name = file_list[[i]]) %>%
      mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                                ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                       ifelse(str_detect(sample_name, "D"), "Diesel", "Gas"))))
  }
  return(list(df_list_filter_area, df_list_removed_area))
}
  
# Notin function
`%notin%` <- Negate(`%in%`)

# Grouping compounds based on RT1, RT2, and Ion1
grouping_comp <- function(data, rt1thres, rt2thres, ion1thres, ion2thres) {

  # create empty list, each sub-list is a compound group with following criteria:
  # RT1 threshold +- 0.4
  # RT2 threshold +- 0.1
  # Ion1 threshold +-0.1
  
  # Initialize the compound column filled with NA values
  data$compound <- NA
  i <- 1
  for (row in 1:nrow(data)) {
    # filter data by index, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    rt1 <- data[row,]$RT1
    rt2 <- data[row,]$RT2
    ion1 <- data[row,]$Ion1
    ion2 <- data[row,]$Ion2
    idx <- which(data$RT1 <= (rt1 + rt1thres) & data$RT1 >= (rt1 - rt1thres) & 
                 data$RT2 <= (rt2 + rt2thres) & data$RT2 >= (rt2 - rt2thres) & 
                 data$Ion1 <= (ion1 + ion1thres) & data$Ion1 >= (ion1 - ion1thres) & 
                 data$Ion2 <= (ion2 + ion2thres) & data$Ion2 >= (ion2 - ion2thres) &
                 is.na(data$compound))
    if (identical(idx, integer(0))) {
      next
    }
    else {
      data[idx, "compound"] <- paste0("Compound_", i, ".")
      i <- i + 1
    }
  }
  return(data)
}

# Filtering similar and unique compound
comp_filter <- function(data, n) {
  all_similar_compounds_idx <- c()
  all_other_compounds_idx <- c()
  all_unique_compounds_idx <- c()
  
  for (comp_grp in unique(data$compound)) {
    # filter data by indexing, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    idx <- which(grepl(comp_grp, data$compound, fixed = TRUE))
    
    if (length(unique(data[idx,]$sample_name)) > (n - 1)) {
      all_similar_compounds_idx <- c(all_similar_compounds_idx, idx)
    }
    else if (length(unique(data[idx,]$sample_name)) < 2) {
      all_unique_compounds_idx <- c(all_unique_compounds_idx, idx)
    }
    else {
      all_other_compounds_idx <- c(all_other_compounds_idx, idx)
    }
  }
  return(list(all_similar_compounds_idx, all_other_compounds_idx, all_unique_compounds_idx))
}

# Pair wise comparison
pairwise_test <- function(df, p_val_threshold, test_choice = "t.test") {
  pairwise_test <- list()
  i <- 1
  for (com_grp in unique(df$compound)) {
    templist <- list()
    # iterates through every fuel_type
    for (fuel1 in unique(df$fuel_type)) { #  sample1 in unique(df$sample_name) / fuel1 in unique(df$fuel_type) 
      idx1 <- which(df$fuel_type == fuel1 & df$compound == com_grp) #  df$sample_name == sample1 / df$fuel_type == fuel1
      if (length(idx1) < 2) {
        next
      }
      else {
        for (fuel2 in unique(df$fuel_type)) { #  sample2 in unique(df$sample_name) / fuel2 in unique(df$fuel_type)
          if (fuel1 == fuel2) { # sample1 == sample2 / fuel1 == fuel2
            next
          }
          else {
            idx2 <- which(df$fuel_type == fuel2 & df$compound == com_grp) #  df$sample_name == sample2 / df$fuel_type == fuel2
            if (length(idx2) < 2) {
              next
            }
            else {
              if (test_choice == "ks") {
                templist[paste0(fuel1, "-", fuel2)] <- ks.test(x = df[idx1,]$Percent_Area,
                                                               y = df[idx2,]$Percent_Area,
                                                               alternative = "two.sided")$p.value
              }
              else if (test_choice == "mn") {
                templist[paste0(fuel1, "-", fuel2)] <- wilcox.test(x = df[idx1,]$Percent_Area,
                                                                   y = df[idx2,]$Percent_Area,
                                                                   alternative = "two.sided")$p.value
              }
              else {
                templist[paste0(fuel1, "-", fuel2)] <- t.test(x = df[idx1,]$Percent_Area,
                                                              y = df[idx2,]$Percent_Area)$p.value
              }
            }
          }
        }
      }
    }
    
    if (any(templist > p_val_threshold)) {
      next
    }
    else {
      # pairwise_test[[i]] <- templist
      # names(pairwise_test)[i] <- com_grp
      pairwise_test[com_grp] <- templist
    }
    i <- i + 1
  }
  return(pairwise_test)
}

# Probabilistic Quotient Normalization
pqn <- function(X, n = "median", QC = NULL) {
  X.norm <- matrix(nrow = nrow(X), ncol = ncol(X))
  colnames(X.norm) <- colnames(X)
  rownames(X.norm) <- rownames(X)
  
  if (!is.null(QC)) {
    # if QC vector exists, use this as reference spectrum
    if (length(QC) == 1) {
      # only 1 reference sample given
      mX <- as.numeric(X[QC, ])
    } else {
      if (n == "mean") {
        mX <- as.numeric(colMeans(X[QC, ]))
      }
      if (n == "median") {
        mX <- as.numeric(apply(X[QC, ], 2, median))
      }
    }
  } else {
    # otherwise use the mean or median of all samples as reference sample
    if (n == "mean") {
      mX <- as.numeric(colMeans(X))
    }
    if (n == "median") {
      mX <- as.numeric(apply(X, 2, median))
    }
  }
  
  # do the actual normalisation
  for (a in 1:nrow(X)) {
    X.norm[a, ] <- as.numeric(X[a, ] / median(as.numeric(X[a, ] / mX)))
  }
  
  return(X.norm)
}

# Relative log abundance plots
RlaPlots <- function(inputdata, type=c("ag", "wg"), cols=NULL,
                     cex.axis=0.8, las=2, ylim=c(-2, 2), oma=c(7, 4, 4, 2) + 0.1, ...) {
  type <- match.arg(type)
  groups <- factor(inputdata[, 1], levels = unique(inputdata[, 1]))
  unique.groups <- levels(groups)
  if (is.null(cols)) 
    cols <- ColList(length(unique.groups))
  box_cols <- c(rep(NA, length(rownames(inputdata))))
  for (ii in 1:length(inputdata[, 1])) 
    box_cols[ii] <- cols[which(unique.groups == inputdata[, 1][ii])]
  
  # Within groups
  if(type == "wg") {
    out_data<-data.frame()
    for (grp in unique.groups) {
      submat <- inputdata[which(inputdata[, 1] == grp), -1]
      med_vals <- apply(submat, 2, median)
      swept_mat <- sweep(submat, 2, med_vals, "-")
      out_data <- rbind(out_data, swept_mat)
    }
    # Across groups (i.e. type == "ag")
  } else  {
    med_vals <- apply(inputdata[, -1], 2, median)
    out_data <- sweep(inputdata[, -1], 2, med_vals, "-")
  }
  
  boxplot(t(out_data),
          cex.axis=cex.axis,                 # font size
          las=las,                           # label orientation
          col=box_cols,                      # colours
          ylim=ylim,                         # y-axis range
          oma=oma,                           # outer margin size
          ...
  )
  
  abline(h=0)
}

