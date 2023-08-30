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

### Gas versus Diesel

![](Embedded visualizations/figure-markdown_github/unnamed-chunk-4-1.png)

### Between groups of gas stations

![](Embedded visualizations/figure-markdown_github/unnamed-chunk-5-1.png)

## Multiple Wilcoxon tests with p-value correction for multiple testing

### Research Question 1: Gas versus Diesel

Before running multiple Wilcoxon tests, it is recommended to examine
whether criteria for univariate parametric tests, such as t-tests, are
violated. If yes, then it is safe to proceed using non-parametric
univariate test, such as Wilcoxon test.

First, equal variance and normally distributed between two populations,
here are Gas and Diesel using histogram, and Q-Q plots. Here, both Gas and
Diesel populations are NOT normally distributed.

``` r
GasData <- as.vector(t(cat_5[,c(1:21)])) # 100149 data points
DieselData <- as.vector(t(cat_5[,c(22:25)])) # 19076 data point

# Histogram
hist(GasData, col='steelblue', main='Gas')
```

![](Embedded visualizations/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
hist(DieselData, col='steelblue', main='Diesel')
```

![](Embedded visualizations/figure-markdown_github/unnamed-chunk-7-2.png)

``` r
# Q-Q plots aka. Normal Probability plots
stats::qqnorm(GasData, main='Gas')
stats::qqline(GasData)
```

![](Embedded visualizations/figure-markdown_github/unnamed-chunk-7-3.png)

``` r
stats::qqnorm(DieselData, main='Diesel')
stats::qqline(DieselData)
```

![](Embedded visualizations/figure-markdown_github/unnamed-chunk-7-4.png)

Then, equality of variance between Gas and Diesel populations can be
examined using Levene’s test and Fligner-Killeen test for non-normally
distributed data. Here, for both tests, p values are \< 0.05, and thus,
there is significant difference in variances between Gas and Diesel
populations.

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##         Df F value    Pr(>F)    
    ## group    1  34.244 5.016e-09 ***
    ##       9948                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## 
    ##  Fligner-Killeen test of homogeneity of variances
    ## 
    ## data:  data by group
    ## Fligner-Killeen:med chi-squared = 531.7, df = 1, p-value < 2.2e-16

Here, multiple Wilcoxon tests followed by p-value correction for
multiple testing were done. Different methods for p-value correction from
function *p.adjust* from package **stats** were used. After p-value
correction, the threshold for p-value can be set (for example,p.adjust
\< 0.05 or \< 0.1), to see how it affects the number of significant
compounds that can be found. Importantly, the threshold should be set so
that it can include as many ASTM reference compounds in the list of
significant compounds as possible. For example here, no significant
compounds can be found with adjusted p-value \< 0.05. But when the adjusted
p-value is \< 0.1, 127 significant compounds are found.

    ## [1] "Number of significant compounds with adjusted p-value < 0.05 = 0"

    ## [1] "Number of significant compounds adjusted p-value < 0.1 = 127"

### Research Question 2: Distinguishing samples from groups of gas stations

Here again, before running multiple Wilcoxon tests, it is recommended to
examine whether criteria for univariate parametric tests, such as t-tests,
are violated. If yes, then it is safe to proceed using non-parametric
univariate test, such as Wilcoxon test.

First, equal variance and normally distributed between two populations,
here are Gas and Diesel using histogram, and Q-Q plots. Here, both Gas
station group 1 (1, 3, 8) and Gas station group 2 (5, 7, 9) populations
are NOT normally distributed.

Equality of variance between Gas station group 1 (1, 3, 8) and Gas
station group 2 (5, 7, 9) populations can be examined using Levene’s
test and Fligner-Killeen test for non-normally distributed data. Here,
for both tests, p values are \< 0.05, and thus, there is a significant
difference in variances between Gas station group 1 (1, 3, 8) and Gas
station group 2 (5, 7, 9) populations.

![](Embedded visualizations/figure-markdown_github/unnamed-chunk-11-1.png)![](Embedded visualizations/figure-markdown_github/unnamed-chunk-11-2.png)![](Embedded visualizations/figure-markdown_github/unnamed-chunk-11-3.png)![](Embedded visualizations/figure-markdown_github/unnamed-chunk-11-4.png)

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##          Df F value Pr(>F)
    ## group     1  0.6128 0.4338
    ##       11149

    ## 
    ##  Fligner-Killeen test of homogeneity of variances
    ## 
    ## data:  data by group
    ## Fligner-Killeen:med chi-squared = 713.21, df = 1, p-value < 2.2e-16

Here, multiple Wilcoxon tests followed by p-value correction for
multiple testing were done. Different methods for p-value correction from
function *p.adjust* from package **stats** were used. After p-value
correction, the threshold for p-value can be set (for example,p.adjust
\< 0.05 or \< 0.1), to see how it affects the number of significant
compounds that can be found. Importantly, the threshold should be set so
that it can include as many ASTM reference compounds in the list of
significant compounds as possible. For example here, no significant
compounds can be found with adjusted p-value \< 0.05. But when the adjusted
p-value is \< 0.1, 127 significant compounds are found.

    ## [1] "Number of significant compounds with adjusted p-value < 0.05 = 75"

    ## [1] "Number of significant compounds adjusted p-value < 0.1 = 94"
