#' @details
#' The data in this package is published in the publication Bhuva et al., _Genome Biology_, 2024. Annotations for each dataset are obtained independently of the molecular measurements. The lowest level of measurement is annotated for each dataset. For BGI STOmics, this is each DNA nanoball, and for 10x Xenium and NanoString CosMx, these are the individual transcript detections.
#' 
#' The package provides the \code{\link{tx2spe}} function to convert these unit measurements into higher level summaries such as cells, regions, square bins, or hex bins for analysis purposes. These summaries are stored in SpatialExperiment objects.
#' 
#' @keywords internal
"_PACKAGE"

#' Sample annotated sub-cellular localised spatial 'omics dataset
#'
#' This object stores a sample transcript table to demonstrate the format of data stored within this package. The sample dataset contains 100,000 transcript detections from 2 breast cancer samples profiled using the 10x Xenium platform.
#'
#' @format A data.frame where each row represents a transcript detection. The full description of these transcript tables is available in the vignettes or in the documentation for the `tx2spe()` function.
#' @docType data
#'
"tx_small"