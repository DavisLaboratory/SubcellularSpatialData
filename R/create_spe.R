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
tx2spe <- function(x, bin = c('cell', 'hex'), nbins = 100) {
  bin = match.arg(bin)

  if (bin %in% 'cell') {
    x = dplyr::rename(x, bin_id = cell)
  } else {
    x = x |> 
      dplyr::group_by(sample) |> 
      dplyr::mutate(bin_id = allocateHex(x, y, nbins)) |> 
      dplyr::group_by(bin_id, .add = TRUE) |> 
      dplyr::mutate(ncell = length(stats::na.omit(unique(cell)))) |> 
      dplyr::ungroup()
  }

  #create bin annotation
  #process character annotations
  bin_annot = x |> 
    dplyr::filter(!is.na(bin_id)) |> 
    dplyr::select(!c(gene, genetype, counts)) |> 
    dplyr::group_by(sample, bin_id, dplyr::across(dplyr::where(is.character))) |> 
    dplyr::summarise(N = dplyr::n()) |> 
    dplyr::group_by(dplyr::across(c(-N))) |> 
    dplyr::arrange(N) |> 
    dplyr::ungroup() |> 
    dplyr::filter(!duplicated(paste(sample, bin_id))) |> 
    dplyr::select(!N)
  
  #process numeric annotations
  bin_annot = x |> 
    dplyr::filter(!is.na(bin_id)) |> 
    dplyr::select(!c(gene, genetype, counts)) |> 
    dplyr::group_by(sample, bin_id) |> 
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), ~ median(.x))) |> 
    dplyr::full_join(bin_annot, by = dplyr::join_by(sample, bin_id)) |> 
    dplyr::mutate(rname = paste(sample, bin_id, sep = '_')) |> 
    dplyr::rename(sample_id = sample) |> 
    as.data.frame()
  #add rownames
  rownames(bin_annot) = bin_annot$rname
  bin_annot = bin_annot |> dplyr::select(!rname)

  #change bin_id to cell/hex id
  colnames(bin_annot)[which(colnames(bin_annot) %in% 'bin_id')] = ifelse(bin == 'cell', 'cell_id', 'hex_id')
  
  #create gene annotation
  gene_annot = x |> 
    dplyr::filter(!is.na(bin_id)) |> 
    dplyr::select(c(gene, genetype)) |> 
    dplyr::filter(!duplicated(paste(gene)))

  #create counts
  counts = x |> 
    dplyr::filter(!is.na(bin_id)) |> 
    dplyr::mutate(
      rname = as.factor(paste(sample, bin_id, sep = '_')),
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
#' @param y a numeric, containing x-coordinates
#' @param bins a numeric, stating the number of bins partitioning the range of x 
#'
#' @return a numeric, containing the index of bins to which each point was allocated
#'
#' @examples
allocateHex <- function(x, y, bins = 30) {
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