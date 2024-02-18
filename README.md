
# SubcellularSpatialData

<!-- badges: start -->
<!-- badges: end -->

The data in this package is published in the publication Bhuva et al., _Genome Biology_, 2024. Annotations for each dataset are obtained independently of the molecular measurements. The lowest level of measurement is annotated for each dataset. For BGI STOmics, this is each DNA nanoball, and for 10x Xenium and NanoString CosMx, these are the individual transcript detections. The package provides functions to convert these unit measurements into higher level summaries such as cells, regions, square bins, or hex bins for analysis purposes. These summaries are stored in SpatialExperiment objects.

## Installation

The release version can be downloaded as follows:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SubcellularSpatialData")
```

You can install the development version of SubcellularSpatialData from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("DavisLaboratory/SubcellularSpatialData")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SubcellularSpatialData)
## basic example code
```

