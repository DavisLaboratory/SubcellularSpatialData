library(SpatialExperiment)
library(dplyr)

#download data from websites and load as table
#  - added annotation column where each detection is mapped onto annotations
#  - annotations are stored in geojson files (see paper for description on how they
#    were obtained Bhuva et al. 2023)

#xenium - brain
df = readRDS('data-local/Xenium_Mar23.rds') |> 
  dplyr::rename(cell = cell_id, gene = feature_name, x = x_location, y = y_location) |> 
  mutate(
    genetype = case_when(
      grepl('BLANK', gene) ~ 'NegBlank',
      grepl('NegControlProbe', gene) ~ 'NegPrb',
      grepl('NegControlCodeword', gene) ~ 'NegCW',
      TRUE ~ 'Gene'
    ),
    cell = ifelse(cell == -1, NA_integer_, cell),
    counts = 1,
    technology = 'Xenium'
  ) |> 
  dplyr::rename(region = RegionID, level = Level) |> 
  dplyr::rename(sample = 'sample_id') |> 
  relocate(sample_id, cell, gene, genetype, x, y, counts, region, technology, level, Level0:Level11)
saveRDS(df, 'data-local/xenium_mm_brain.rds')

#xenium - breast
df = readRDS('data-local/Xenium_breast_addon_tx_annotated.rds') |> 
  dplyr::rename(cell = cell_id, gene = feature_name, x = x_location, y = y_location) |> 
  mutate(
    genetype = case_when(
      grepl('BLANK', gene) ~ 'NegBlank',
      grepl('NegControlProbe', gene) ~ 'NegPrb',
      grepl('NegControlCodeword', gene) ~ 'NegCW',
      TRUE ~ 'Gene'
    ),
    cell = ifelse(cell == -1, NA_integer_, cell),
    counts = 1,
    technology = 'Xenium'
  ) |> 
  dplyr::rename(region = RegionID) |> 
  dplyr::rename(sample_id = 'sample') |> 
  relocate(sample_id, cell, gene, genetype, x, y, counts, region, technology)
saveRDS(df, 'data-local/xenium_hs_breast_addon.rds')

#stomics
df = readRDS('data-local/STOmics_Brain.rds') |> 
  select(!UMICount) |> 
  dplyr::rename(gene = geneID, cell = label, counts = MIDCounts) |> 
  mutate(technology = 'STOmics', genetype = 'Gene') |> 
  dplyr::rename(region = RegionID, level = Level) |> 
  dplyr::rename(sample = 'sample_id') |> 
  relocate(sample_id, cell, gene, genetype, x, y, counts, region, technology, level, Level0:Level11)
saveRDS(df, 'data-local/stomics_mm_brain.rds')

#cosmx
df = readRDS('data-local/CosMx_NSCLC_tx_annotated.rds') |> 
  select(!c(sample, x, y)) |> 
  dplyr::rename(
    sample = TMA,
    x = x_global_px,
    y = y_global_px,
    gene = target,
    cell = cell_ID
  ) |>
  dplyr::mutate(cell = ifelse(cell == 0, NA, cell), cell = paste(fov, cell, sep = '_')) |>
  dplyr::select(sample, x, y, gene, cell, fov, Level1, Level2, RegionID) |>
  mutate(
    technology = 'CosMx',
    genetype = ifelse(grepl('NegPrb', gene), 'NegPrb', 'Gene'),
    cell = ifelse(grepl('NA$', cell), NA_character_, cell),
    counts = 1
  ) |> 
  mutate(RegionID = case_when(
    RegionID %in% c('Tumour area with stroma', 'Epithelium - panCK+', 'Tumour area') ~ 'Tumour',
    RegionID %in% c('Abnormal but not tumour') ~ 'Abnormal epithelium',
    RegionID %in% 'Fibrotic area' ~ 'Fibrotic',
    RegionID %in% 'Necrotic area' ~ 'Necrotic',
    RegionID %in% 'Abnormal area' ~ 'Abnormal',
    TRUE ~ RegionID
  )) |> 
  dplyr::rename(region = RegionID) |> 
  dplyr::rename(sample = 'sample_id') |> 
  relocate(sample_id, cell, gene, genetype, x, y, counts, region, technology, Level1:Level2)
saveRDS(df, 'data-local/cosmx_hs_nsclc.rds')
