SpaTopic
=======

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/xiyupeng/SpaTopic/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/xiyupeng/SpaTopic/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

An R package for fast topic inference to identify tissue architecture in multiplexed images.
It implements a spatial topic model to identify immunologic topics across multiplexed images, given the cell location and cell type information as input.
Collapsed Gibbs Sampling algorithm is used for model inference.
Compared to other KNN-based methods (such as KNN-kmeans, the default in Seurat v5 R package), SpaTopic runs much faster on large-scale image dataset.


## Installation

The R package `SpaTopic` now is available in CRAN and can be installed with the following code.

``` r
install.packages("SpaTopic")
```

The development version of `SpaTopic` can be installed from the [GitHub repository](https://github.com/xiyupeng/SpaTopic).

``` r
# install.packages("devtools")
devtools::install_github("xiyupeng/SpaTopic")
```

## Dependency

SpaTopic requires dependency on the following R packages:

- [Rcpp]( https://cran.r-project.org/package=Rcpp)  for C++ codes
- [RcppProgress](https://cran.r-project.org/package=RcppProgress) for C++ codes
- [RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo) for C++ codes
- [RANN](https://cran.r-project.org/package=RANN) for fast KNN 
- [foreach](https://cran.r-project.org/package=foreach)  for parallel computing
- [sf](https://cran.r-project.org/package=sf) for spatial analysis

## Usage

The required input of SpaTopic is a data frame containing cells within on a single image or a list of data frames for multiple images. Each data frame consists of four columns: The image ID, X, Y cell coordinates, and cell type information. 

``` r
library(SpaTopic)
library(sf)
## The input can be a data frame or a list of data frames
data("lung5")
head(lung5)
#     image        X        Y           type
#1_1 image1 4215.889 158847.7      Dendritic
#2_1 image1 6092.889 158834.7     Macrophage
#3_1 image1 7214.889 158843.7 Neuroendocrine
#4_1 image1 7418.889 158813.7     Macrophage
#5_1 image1 7446.889 158845.7     Macrophage
#6_1 image1 3254.889 158838.7          CD4 T
gibbs.res<-SpaTopic_inference(lung5, ntopics = 7, sigma = 50, region_radius = 400)
```

For detailed usage of SpaTopic,
please check the [tutorial](https://xiyupeng.github.io/SpaTopic/).

## Data

The example image used in the tutorial can be downloaded from [here](https://drive.google.com/drive/folders/1_mJUjzQXWgUZlwUaLq0HKxX-aqgiQ8eD?usp=sharing).
It is stored in a Seurat v5 object. 

## Output

``` r
str(gibbs.res)
# List of 8
#  $ Perplexity   : num 11.3
#  $ Deviance     : num 485960
#  $ loglikelihood: num -242980
#  $ Beta         :'data.frame':	38 obs. of  7 variables:
#   ..$ topic1: num [1:38] 0.03587 0.02539 0.00755 0.01858 0.02585 ...
#   ..$ topic2: num [1:38] 6.51e-03 3.55e-02 2.62e-06 5.80e-04 7.75e-01 ...
#   ..$ topic3: num [1:38] 4.54e-06 4.54e-06 9.13e-04 3.45e-01 1.73e-03 ...
#   ..$ topic4: num [1:38] 0.02664 0.01743 0.00186 0.0152 0.08919 ...
#   ..$ topic5: num [1:38] 2.99e-06 2.99e-06 5.32e-03 1.91e-02 4.90e-03 ...
#   ..$ topic6: num [1:38] 6.35e-06 6.35e-06 2.04e-02 3.43e-03 6.35e-06 ...
#   ..$ topic7: num [1:38] 0.00534 0.00699 0.00604 0.01843 0.00655 ...
#  $ Theta        : num [1:971, 1:7] 0.855601 0.000232 0.999269 0.99889 0.998725 ...
#  $ Ndk          : int [1:971, 1:7] 107 0 82 54 47 72 100 0 0 0 ...
#  $ Nwk          : int [1:38, 1:7] 390 276 82 202 281 505 697 522 29 58 ...
#  $ Z.trace      :'data.frame':	100149 obs. of  7 variables:
#   ..$ topic1: num [1:100149] 0.065 0.865 0.135 0.82 0.785 0.02 0.1 0.105 0.075 0.095 ...
#   ..$ topic2: num [1:100149] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ topic3: num [1:100149] 0.275 0.005 0.21 0.005 0.005 0.77 0.02 0.015 0.085 0.075 ...
#   ..$ topic4: num [1:100149] 0.415 0 0 0.01 0.005 0.1 0.665 0.62 0.015 0.025 ...
#   ..$ topic5: num [1:100149] 0.005 0.01 0 0 0 0 0.005 0.005 0.005 0 ...
#   ..$ topic6: num [1:100149] 0 0 0.655 0.165 0.205 0.005 0 0 0 0 ...
#   ..$ topic7: num [1:100149] 0.24 0.12 0 0 0 0.105 0.21 0.255 0.82 0.805 ...
```

## Citation

Xiyu Peng, James W. Smithy, Mohammad Yosofvand, Caroline E. Kostrzewa, MaryLena Bleile, Fiona D. Ehrich, Jasme Lee, Michael A. Postow, Margaret K. Callahan, Katherine S. Panageas, Ronglai Shen. Decoding Spatial Tissue Architecture: A Scalable Bayesian Topic Model for Multiplexed Imaging Analysis.
bioRxiv. doi: https://doi.org/10.1101/2024.10.08.617293

## Contact

If you have any problems, please contact:

Xiyu Peng (pansypeng124@gmail.com, pengx1@mskcc.org)

