SpaTopic
=======

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/xiyupeng/SpaTopic/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/xiyupeng/SpaTopic/actions/workflows/R-CMD-check.yaml)
   [![](https://cranlogs.r-pkg.org/badges/grand-total/SpaTopic)](https://cran.r-project.org/package=SpaTopic)
  [![](https://cranlogs.r-pkg.org/badges/SpaTopic)](https://cran.r-project.org/package=SpaTopic)
  <!-- badges: end -->

An R package for fast topic inference to identify tissue architecture in multiplexed images.
It implements a **novel spatial topic model** to identify highly interpretable immunologic topics across **multiple** multiplexed images, simply given the cell location and cell type information as input.

In the R package, we adapt an approach originally developed for image segmentation in computer vision, incorporating spatial information into the flexible design of regions (image partitions, analogous to documents in language modeling).
We further refined the approach to address unique challenges in cellular images and provide an efficient C++ implementation of the algorithm in this R package.

Compared to other KNN-based methods (such as KNN-kmeans, the default neighborhood analysis in Seurat v5 R package), SpaTopic runs much faster on large-scale image dataset with minimal memory usage. 
For example, when working on Nanostring CosMx NSCLC dataset, Spatopic can spatial cluster **0.1 million** of cells on a single image within **1 min** on a regular Mac Air (See [tutorial](https://xiyupeng.github.io/SpaTopic/).


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
packageVersion("SpaTopic")
#> [1] '1.2.0'
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

## Example Data

The example image used in the tutorial can be downloaded from [here](https://drive.google.com/drive/folders/1_mJUjzQXWgUZlwUaLq0HKxX-aqgiQ8eD?usp=sharing).
It is stored in a Seurat v5 object. 

## Example Output

The algorithm generates two key statistics for further analysis: 

- 1) topic content, a spatially-resolved topic distribution
over cell types, and
- 2) topic assignment for each cell within images.

**Topic Spatial Distribution over images**
<div>
<img src="https://github.com/user-attachments/assets/0f116b96-6afc-473a-acbf-8137bdf54c2f" width="500" height="600"/>
</div>

**Topic Content**
<div>
<img src="https://github.com/user-attachments/assets/c80dc4b3-5388-409a-8fa1-f3c975627771" width="600" height="300"/>
</div>

## Citation

Xiyu Peng, James W. Smithy, Mohammad Yosofvand, Caroline E. Kostrzewa, MaryLena Bleile, Fiona D. Ehrich, Jasme Lee, Michael A. Postow, Margaret K. Callahan, Katherine S. Panageas, Ronglai Shen. Decoding Spatial Tissue Architecture: A Scalable Bayesian Topic Model for Multiplexed Imaging Analysis.
bioRxiv. doi: https://doi.org/10.1101/2024.10.08.617293

## Contact

If you have any problems, please contact:

Xiyu Peng (pansypeng124@gmail.com, pengx@stat.tamu.edu)

