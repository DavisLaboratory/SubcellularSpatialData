#' Compute cell/bin-level summaries and save to a SpatialExperiment object
#'
#' This function summarises transcript detection (or DNA nanoball spot counts) into cells or hexbins and returns a SpatialExperiment object. Annotations of cells or hexbins are summarised such that the most frequent annotation is used.
#'
#' @param x a data.frame containing a sub-cellular localised dataset from the SubcellularSpatialData package.
#' @param bin a character, stating whether transcripts should be binned into cells or hexbins. The default is 'cell', however, it is advisable to use 'hexbin' when cell identification is not accurate (e.g., when cell boundaries are inferred only from nuclei and not directly measured using cytoplasmic/cell membrane stains).
#' @param nbins a numeric, stating the number of bins to create in the x and y axes
#'
#' @return a SpatialExperiment object containing cell/hexbin-level expression quantification.
#' @examples
#' @export 
#' 
tx2spe <- function(x, bin = c('cell', 'hex', 'square', 'region'), nbins = 100) {
  #checks
  bin = match.arg(bin)
  if (nbins <= 0) stop("Value for 'nbins' should be greater than 0")
  if (bin == 'region' & !'region' %in% colnames(x)) stop("Required column 'region' not found in input data 'x'")

  if (bin == 'cell') {
    warning('Please note that cell segmentation may not be accurate in this data and it may be better to use other binning strategies.')
    x = dplyr::rename(x, bin_id = cell)
  } else if (bin %in% c('hex', 'square', 'region')) {
    allocateFun = switch(
      bin,
      hex = allocateHex,
      square = allocateSquare,
      region = allocateRegion
    )

    x = x |> 
      dplyr::group_by(sample_id) |> 
      dplyr::mutate(bin_id = allocateFun(x, y, bins = nbins, regions = region)) |> 
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

  #create bin annotation
  #process character annotations
  bin_annot = x |> 
    dplyr::filter(!is.na(bin_id)) |> 
    dplyr::select(!c(gene, genetype, counts)) |> 
    dplyr::group_by(sample_id, bin_id, dplyr::across(dplyr::where(is.character))) |> 
    dplyr::summarise(N = dplyr::n()) |> 
    dplyr::group_by(dplyr::across(c(-N))) |> 
    dplyr::arrange(N) |> 
    dplyr::ungroup() |> 
    dplyr::filter(!duplicated(paste(sample_id, bin_id))) |> 
    dplyr::select(!N)
  
  #process numeric annotations
  bin_annot = x |> 
    dplyr::filter(!is.na(bin_id)) |> 
    dplyr::select(!c(gene, genetype, counts)) |> 
    dplyr::group_by(sample_id, bin_id) |> 
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), ~ mean(.x))) |> 
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
#'
#' @return a numeric, containing the index of bins to which each point was allocated
#'
#' @examples
allocateHex <- function(x, y, bins = 30, ...) {
  stopifnot(length(x) == length(y))
  stopifnot(bins > 0)
  
  #compute hexbins
  hx = hexbin::hexbin(x, y, xbins = bins, IDs = TRUE)
  hix = hx@cID
  hxy = hexbin::hcell2xy(hx)
  #remap hex cell ids
  hixmap = 1:length(hx@cell)
  names(hixmap) = hx@cell
  hix = hixmap[as.character(hix)]
  
  return(hix)
}

#' Allocate points to a regular square grid
#'
#' @param x a numeric, containing x-coordinates
#' @param y a numeric, containing y-coordinates
#' @param bins a numeric, stating the number of bins partitioning the range of x 
#'
#' @return a numeric, containing the index of bins to which each point was allocated
#'
#' @examples
allocateSquare <- function(x, y, bins = 30, ...) {
  stopifnot(length(x) == length(y))
  stopifnot(bins > 0)
  
  #compute bin widths
  bw = diff(range(x)) / bins
  x = ceiling((x - min(x)) / bw)
  y = ceiling((y - min(y)) / bw)
  x[x == 0] = 1
  y[y == 0] = 1
  sqix = paste(x, y, sep = '_')

  return(sqix)
}

#' Allocate points to annotated regions
#'
#' @param x a numeric, containing x-coordinates
#' @param y a numeric, containing y-coordinates
#' @param regions a character, specifying the regions
#'
#' @return a numeric, containing the index of bins to which each point was allocated
#'
#' @examples
allocateRegion <- function(x, y, regions, ...) {
  stopifnot(length(x) == length(y))
  stopifnot(length(x) == length(regions))
  
  return(regions)
}