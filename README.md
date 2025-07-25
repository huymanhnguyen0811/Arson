## Documentation

This repo is accompanying the publication: “The Use of Computational
Fingerprinting Techniques to Distinguish Sources of Accelerants Used in
Wildfire Arson”.

Users need to first install R with this
[link](https://cran.r-project.org/mirrors.html) and Rstudio with this
[link](https://posit.co/download/rstudio-desktop/).

This workflow ran on Windows 11 OS 11th Gen Intel(R) Core(TM) i7-11800H
@ 2.30GHz, 16 GB RAM;

The RStudio version used in this demo is 2023.06.0+421 “Mountain
Hydrangea” Release for Windows;

The R version used in this demo is 4.3.1

## Data processing

First, the following R packages are installed and loaded in the global
the environment along with in-house built functions to minimize
repetitiveness in the code.

Details about these functions can be found in Data processing & Normalization.R file in
this repo.

## PCA

Further details on the package and PCA functions used can be found in folder demo analysis, under the following files: demo analysis_stats_RQ1 & RQ2.R and demo_analysis_PCA_RQ1 & RQ2.R

    ## [1] "Number of significant compounds with adjusted p-value < 0.05 = 75"

    ## [1] "Number of significant compounds adjusted p-value < 0.1 = 94"
