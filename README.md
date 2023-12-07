SpaTopic
=======

An R package for fast topic inference to identify tissue architecture in multiplexed images.
It runs much faster on large-scale image data, compared to the default KNN-Kmeans method in the Seurat package.


## Installation

`SpaTopic` can be installed from the GitHub repository using the devtools package.

``` r
# Install SpaTopic
devtools::install_github("xiyupeng/SpaTopic")
library(SpaTopic)
```

## Usage

The required input of SpaTopic is a data frame containing cells within on a single image or a list of data frames for multiple images. Each data frame consists of four columns: The image ID, X, Y cell coordinates, and cell type. 
You can find the example input data under the data folder.

``` r
head(dataset)
#     image        X        Y           type
#1_1 image1 4215.889 158847.7      Dendritic
#2_1 image1 6092.889 158834.7     Macrophage
#3_1 image1 7214.889 158843.7 Neuroendocrine
#4_1 image1 7418.889 158813.7     Macrophage
#5_1 image1 7446.889 158845.7     Macrophage
#6_1 image1 3254.889 158838.7          CD4 T
gibbs.res<-SpaTopic_inference(dataset, ntopics = 7, sigma = 50, region_radius = 400)
```

For more detailed usage of SpaTopic,
please check the tutorial and the vignette.

## Data

The example image used in the Vignette are provided in: 

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



