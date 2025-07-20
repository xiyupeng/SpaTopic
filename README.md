SpaTopic (SpatialTopic)
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
For example, when working on Nanostring CosMx NSCLC dataset, SpaTopic can spatially cluster **0.1 million** of cells on a single image within **1 min** on a regular Mac Air (See [tutorial](https://xiyupeng.github.io/SpaTopic/articles/SpaTopic.html)).

**News: SpatialTopic now has been published in Nature Communications**:

Peng, X., Smithy, J.W., Yosofvand, M. et al. Scalable topic modelling decodes spatial tissue architecture for large-scale multiplexed imaging analysis. Nat Commun 16, 6619 (2025). https://doi.org/10.1038/s41467-025-61821-y

**Tutorial available to get start it**

https://xiyupeng.github.io/SpaTopic/articles/SpaTopic.html


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

```
print(gibbs.res)
#> SpaTopic Results
#> ----------------
#> Number of topics: 7 
#> Perplexity: 11.31563 
#> 
#> Topic Content(Topic distribution across cell types):
#>                                 topic1       topic2       topic3      topic4
#> Alveolar Epithelial Type 1 0.035870295 6.511503e-03 4.541367e-06 0.026643327
#> Alveolar Epithelial Type 2 0.025386476 3.553900e-02 4.541367e-06 0.017427665
#> Artery                     0.007545591 2.624548e-06 9.128148e-04 0.001856373
#> B                          0.018581190 5.800251e-04 3.446035e-01 0.015203195
#> Basal                      0.025846292 7.753466e-01 1.730261e-03 0.089193312
#>                                  topic5       topic6      topic7
#> Alveolar Epithelial Type 1 2.987411e-06 6.348481e-06 0.005341969
#> Alveolar Epithelial Type 2 2.987411e-06 6.348481e-06 0.006994451
#> Artery                     5.320579e-03 2.044846e-02 0.006041096
#> B                          1.912242e-02 3.434528e-03 0.018434717
#> Basal                      4.902342e-03 6.348481e-06 0.006549552
#> ...
#> 
#> Use $Z.trace for posterior probabilities of topic assignments for each cell
#> Use $cell_topics for final topic assignments for each cell
#> Use $parameters for accessing model parameters
```

For detailed usage of SpaTopic,
please check:

- The home page for the tutorial is available [here](https://xiyupeng.github.io/SpaTopic/).
- A step-by-step tutorial to reproduce the result on Cosmx lung cancer sample is available [here](https://xiyupeng.github.io/SpaTopic/articles/SpaTopic.html).

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

Peng, X., Smithy, J.W., Yosofvand, M. et al. Scalable topic modelling decodes spatial tissue architecture for large-scale multiplexed imaging analysis. Nat Commun 16, 6619 (2025). https://doi.org/10.1038/s41467-025-61821-y

Preprint: Xiyu Peng, James W. Smithy, Mohammad Yosofvand, Caroline E. Kostrzewa, MaryLena Bleile, Fiona D. Ehrich, Jasme Lee, Michael A. Postow, Margaret K. Callahan, Katherine S. Panageas, Ronglai Shen. Decoding Spatial Tissue Architecture: A Scalable Bayesian Topic Model for Multiplexed Imaging Analysis.
bioRxiv. doi: https://doi.org/10.1101/2024.10.08.617293

## Contact

If you have any problems, please contact:

Xiyu Peng (pansypeng124@gmail.com, pengx@stat.tamu.edu)

