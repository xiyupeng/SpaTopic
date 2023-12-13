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

`SpaTopic` can be installed from the [GitHub repository](https://github.com/xiyupeng/SpaTopic) using the devtools package.

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
#List of 8
# $ Perplexity   : num 11.3
# $ Deviance     : num 485960
# $ loglikelihood: num -242980
# $ Beta         : num [1:38, 1:7] 0.03587 0.02539 0.00755 0.01858 0.02585 ...
# $ Theta        : num [1:971, 1:7] 0.855601 0.000232 0.999269 0.99889 0.998725 ...
# $ Ndk          : int [1:971, 1:7] 107 0 82 54 47 72 100 0 0 0 ...
# $ Nwk          : int [1:38, 1:7] 390 276 82 202 281 505 697 522 29 58 ...
# $ Z.trace      : int [1:100149, 1:7] 13 173 27 164 157 4 20 21 15 19 ...
```

## Citation

Coming soon......

## Contact

If you have any problems, please contact:

Xiyu Peng (pansypeng124@gmail.com, pengx1@mskcc.org)

