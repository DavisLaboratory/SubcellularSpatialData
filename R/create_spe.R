#' Compute cell/bin-level summaries and save to a SpatialExperiment object
#'
#' This function summarises transcript detection (or DNA nanoball spot counts) into cells or hexbins and returns a SpatialExperiment object. Annotations of cells or hexbins are summarised such that the most frequent annotation is used.
#'
#' @param x a data.frame containing a sub-cellular localised dataset from the SubcellularSpatialData package.
#' @param bin a character, stating whether transcripts should be binned into cells or hexbins. The default is 'cell', however, it is advisable to use 'hexbin' when cell identification is not accurate (e.g., when cell boundaries are inferred only from nuclei and not directly measured using cytoplasmic/cell membrane stains).
#' @param nbins a numeric, stating the number of bins to create in the x and y axes
#' 
#' @details
#' 
#' ## Data description
#' This function works on all data associated with this package. For new datasets, the transcript table MUST have the following columns:
#' \itemize{
#'  \item{sample_id} {Unique identifier for the sample the transcript belongs to.}
#'  \item{cell} {Unique identifier for the cell the transcript belongs to (NA if it is not allocated to any cell).}
#'  \item{gene} {Gene name or identity of the transcript.}
#'  \item{genetype} {Type of target (e.g., gene or control probe). True genes should be annotated as "Gene". Other target types can be named as desired.}
#'  \item{x} {x-coordinate of the transcript.}
#'  \item{y} {y-coordinate of the transcript.}
#'  \item{counts} {Number of transcripts detected at this location (always 1 for NanoString CosMx and 10x Xenium as each row is an individual detection).}
#'  \item{technology} {Name of the platform the data was sequenced on. "STOmics" is treated differently as the number of spots as well as the number of nuclei/cells are counted. No specific format required for other technologies.} 
#' }
#' 
#' The "region" column is optional and only required if binning by regions. If present, it must be a character or factor. For datasets within this package, this column stored the independently annotated histological region.
#' 
#' ## Summarisation
#' When transcript counts are aggregated (using any approach), aggregation of numeric columns is performed using the `mean(..., na.rm = TRUE)` function and aggregation of character/factor columns is performed such that the most frequent class becomes the representative class. As such, for a hex bin, the highest frequency region for the transcripts allocated to the bin becomes the bin's region annotation. New coordinates for cells, bins, or regions are computed using the mean function as well, therefore represent the center of mass for the object. For hex and square bins, the average coordinate is computed by default, however, x and y indices for the bin are stored in the colData under the *bin_x* and *bin_y* columns. 
#' 
#' When aggregation is perfomed using bins or regions, an additional column, _ncells_, is computed that indicates how many unique cells are present within the bin/region. Do note that if a cell overlaps multiple bins/regions, it will be counted in each bin/region.
#' 
#' Please note that BGI only performs nuclear binning of counts therefore cellular counts obtained will represent nuclear counts only. We recommended that for fair evaluations during benchmarking, alternative binning algorithms that are fair are explored, or analysis be performed on square or hex binned counts.
#' 
#' @return a SpatialExperiment object containing cell/bin/region-level expression quantification.
#' @examples
#' data(tx_small)
#' head(tx_small)
#' tx2spe(tx_small, "region")
#' 
#' @export 
#' 
tx2spe <- function(x, bin = c('cell', 'hex', 'square', 'region'), nbins = 100) {
  required_cols = c("sample_id", "cell", "gene", "genetype", "x", "y", "counts", "technology")

  # checks
  bin = match.arg(bin)
  # missing cols
  missing_cols = setdiff(required_cols, colnames(x))
  if (length(missing_cols) > 0) {
    stop(paste0("The following columns are missing from the provided transcript table 'x': ", paste(missing_cols, collapse = ", ")))
  }
  # invalid bins
  if (nbins <= 0) stop("Value for 'nbins' should be greater than 0")
  if (bin == "region" & !"region" %in% colnames(x)) stop("Required column 'region' missing from the provided transcript table 'x'")
  
  # check types for required cols
  if (!any(x$genetype == "Gene")) stop("No genes found. Ensure genes are marked by \"Gene\" in the \'genetype\' column.")
  if (!is.numeric(x$x)) stop("\"x\" coordinates must be numeric in the transcript table 'x'.")
  if (!is.numeric(x$y)) stop("\"y\" coordinates must be numeric in the transcript table 'x'.")
  if (!is.numeric(x$counts)) stop("\"counts\" coordinates must be numeric in the transcript table 'x'.")

  if (bin == 'cell') {
    warning('Please note that cell segmentation may not be accurate in this data and it may be better to use other binning strategies. Refer to function documentation for details.')
    x = dplyr::mutate(x, cell = as.character(cell))
    x = dplyr::rename(x, bin_id = cell)
  } else if (bin %in% c('hex', 'square', 'region')) {
    allocateFun = switch(
      bin,
      hex = allocateHex,
      square = allocateSquare,
      region = allocateRegion
    )

    #compute binning and the number of cells in the bin
    x = x |> 
      split(x$sample_id) |>
      lapply(\(df) {
        df |>
          cbind(allocateFun(df$x, df$y, bins = nbins, regions = df$region))
      }) |> 
      do.call(what = 'rbind') |> 
      dplyr::group_by(bin_id, .add = TRUE) |> 
      dplyr::mutate(ncells = length(stats::na.omit(unique(cell)))) |> 
      dplyr::ungroup()
    
    if ('technology' %in% colnames(x) & x$technology[1] %in% 'STOmics') {
      x = x |> 
        dplyr::group_by(sample_id, bin_id) |> 
        dplyr::mutate(nspots = sum(!duplicated(paste(x, y, sep = '_')))) |> 
        dplyr::ungroup()
    }
  }

  # create bin annotation
  # process character annotations
  bin_annot = x |>
    dplyr::filter(!is.na(bin_id)) |>
    dplyr::select(!c(gene, genetype)) |>
    dplyr::select(sample_id, bin_id, dplyr::where(is.character)) |>
    dplyr::group_by(sample_id, bin_id) |>
    dplyr::summarise(dplyr::across(everything(),
      ~ ifelse(all(is.na(.x)), NA_character_, table(na.omit(.x)) |> sort() |> names() |> tail(1))
    )) |> 
    dplyr::ungroup()
  
  #process numeric annotations
  bin_annot = x |> 
    dplyr::filter(!is.na(bin_id)) |> 
    dplyr::select(!c(gene, genetype, counts)) |> 
    dplyr::group_by(sample_id, bin_id) |> 
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), ~ mean(.x, na.rm = TRUE))) |> 
    dplyr::full_join(bin_annot, by = dplyr::join_by(sample_id, bin_id)) |> 
    dplyr::mutate(rname = paste(sample_id, bin_id, sep = '_')) |> 
    as.data.frame()
  #add rownames
  rownames(bin_annot) = bin_annot$rname
  bin_annot = bin_annot |> dplyr::select(!rname)

  #change bin_id to cell/hex id
  colnames(bin_annot)[which(colnames(bin_annot) %in% 'bin_id')] = ifelse(bin == 'cell', 'cell_id', 'bin_id')
  
  #create gene annotation
  gene_annot = x |> 
    dplyr::filter(!is.na(bin_id)) |> 
    dplyr::select(c(gene, genetype)) |> 
    dplyr::filter(!duplicated(paste(gene)))

  #create counts
  counts = x |> 
    dplyr::filter(!is.na(bin_id)) |> 
    dplyr::mutate(
      rname = as.factor(paste(sample_id, bin_id, sep = '_')),
      gene = as.factor(gene)
    ) |> 
    dplyr::select(c(rname, gene, counts)) |> 
    dplyr::group_by(rname, gene) |> 
    dplyr::summarise(counts = sum(counts)) |> 
    with(Matrix::sparseMatrix(
      i = as.numeric(gene),
      j = as.numeric(rname),
      x = counts,
      dimnames = list(levels(gene), levels(rname))
    ))
  counts = counts[gene_annot$gene, rownames(bin_annot)]
  
  #create SPE
  spe = SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = bin_annot,
    rowData = gene_annot,
    spatialCoordsNames = c('x', 'y')
  )

  return(spe)
}

#' Allocate points to hexbins
#'
#' @param x a numeric, containing x-coordinates
#' @param y a numeric, containing y-coordinates
#' @param bins a numeric, stating the number of bins partitioning the range of x 
#' @param ... additional args if needed
#'
#' @return a numeric, containing the index of bins to which each point was allocated
#'
#' @keywords internal
allocateHex <- function(x, y, bins = 30, ...) {
  stopifnot(length(x) == length(y))
  stopifnot(bins > 0)
  
  # compute hexbins
  hx = hexbin::hexbin(x, y, xbins = bins, IDs = TRUE)
  hix = hx@cID

  # compute bin coords
  hxy = hexbin::hcell2xy(hx)
  names(hxy$x) = names(hxy$y) = as.character(hx@cell)

  # compute bin area
  # robustly compute height of triangle in hexagon
  h = diff(hx@xbnds) / hx@xbins * 0.5
  area_hex = h ^ 2 * 6 / sqrt(3)

  hix = data.frame(
    bin_id = hix,
    bin_x = hix %/% hx@dimen[2],
    bin_y = hix %% hx@dimen[2],
    coord_x = hxy$x[as.character(hix)],
    coord_y = hxy$y[as.character(hix)],
    area = area_hex
  )

  return(hix)
}

#' Allocate points to a regular square grid
#'
#' @param x a numeric, containing x-coordinates
#' @param y a numeric, containing y-coordinates
#' @param bins a numeric, stating the number of bins partitioning the range of x 
#' @param ... additional args if needed
#'
#' @return a numeric, containing the index of bins to which each point was allocated
#'
#' @keywords internal
allocateSquare <- function(x, y, bins = 30, ...) {
  stopifnot(length(x) == length(y))
  stopifnot(bins > 0)
  
  # compute bin widths
  yrng = diff(range(y))
  bw = diff(range(x)) / bins
  x = ceiling((x - min(x)) / bw)
  y = ceiling((y - min(y)) / bw)
  x[x == 0] = 1
  y[y == 0] = 1
  max_x = max(x)
  max_y = max(y)

  # compute area
  lasty_wd = yrng %% bw
  if (lasty_wd == 0) {
    area_rect = wd ^ 2
  } else {
    area_rect = ifelse(y == max_y, bw * lasty_wd, bw^2)
  }

  sqix = data.frame(bin_id = x + (y - 1) * max_x, bin_x = x, bin_y = y, area = area_rect)

  return(sqix)
}

#' Allocate points to annotated regions
#'
#' @param x a numeric, containing x-coordinates
#' @param y a numeric, containing y-coordinates
#' @param regions a character, specifying the regions
#' @param ... additional args if needed
#'
#' @return a numeric, containing the index of bins to which each point was allocated
#'
#' @keywords internal
allocateRegion <- function(x, y, regions, ...) {
  stopifnot(length(x) == length(y))
  stopifnot(length(x) == length(regions))
  regions = data.frame(bin_id = regions)
  
  return(regions)
}