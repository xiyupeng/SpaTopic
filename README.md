SpaTopic
=======

An R package for fast topic inference to identify tissue architecture in multiplexed images.
It runs much faster on large-scale image data, compared to the default KNN-Kmeans method in the Seurat package.


## Installation

`SpaTopic` can be installed from the GitHub repository using the devtools package.

``` r
# Install SpaTopic
devtools::install_github("xiyupeng/SpaTopic")
```

## Usage

The required input of SpaTopic is a data frame containing cells within on a single image or a list of data frames for multiple images. Each data frame consists of four columns: The image ID, X, Y cell coordinates, and cell type. 
You can find the example input data under the data folder.

``` r
library(SpaTopic)
gibbs.res<-SpaTopic_inference(dataset, ntopics = 7, sigma = 50, region_radius = 400)
str(gibbs.res)
```

For more detailed usage of SpaTopic,
please check the tutorial and the vignette.


## Data

The example image used in the Vignette are provided in: 


## Citation

Coming soon......

## Contact

If you have any problems, please contact:

Xiyu Peng (pansypeng124@gmail.com, pengx1@mskcc.org)



