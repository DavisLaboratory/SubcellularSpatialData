##build record
meta = rbind(
  c(
    Title = 'xenium_mm_brain',
    Description = '10x Xenium dataset of 3 serial sections of fresh frozen mouse brain. Raw transcript level data is provided with region annotations for each detection obtained using image registration of the DAPI image to the Allen Brain Atlas. The Aligning Big Brain Atlases (ABBA) plugin in ImageJ was used for image registration.',
    BiocVersion = as.numeric(as.character(BiocManager::version())) + 0.01,
    Genome = NA,
    SourceType = 'CSV',
    SourceUrl = 'https://www.10xgenomics.com/resources/datasets/fresh-frozen-mouse-brain-replicates-1-standard',
    SourceVersion = '2023-03-09',
    Species = 'Mus musculus',
    TaxonomyId = '10090',
    Coordinate_1_based = TRUE,
    DataProvider = '10x',
    Maintainer = 'Dharmesh D. Bhuva <dharmesh.bhuva@adelaide.edu.au>',
    RDataClass = 'data.frame',
    DispatchClass = 'Rds',
    Location_Prefix = 'https://zenodo.org',
    RDataPath = 'record/7959787/files/xenium_mm_brain.rds'
  ),
  c(
    Title = 'stomics_mm_brain',
    Description = 'BGI STOmics dataset of the mouse brain. Raw DNA nanoball spot level data is provided with region annotations for each spot obtained using image registration of the DAPI image to the Allen Brain Atlas. The Aligning Big Brain Atlases (ABBA) plugin in ImageJ was used for image registration.',
    BiocVersion = as.numeric(as.character(BiocManager::version())) + 0.01,
    Genome = NA,
    SourceType = 'TSV',
    SourceUrl = 'https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/Bin1_matrix/Mouse_brain_Adult_GEM_bin1.tsv.gz',
    SourceVersion = '2023-03-09',
    Species = 'Mus musculus',
    TaxonomyId = '10090',
    Coordinate_1_based = TRUE,
    DataProvider = 'BGI',
    Maintainer = 'Dharmesh D. Bhuva <dharmesh.bhuva@adelaide.edu.au>',
    RDataClass = 'data.frame',
    DispatchClass = 'Rds',
    Location_Prefix = 'https://zenodo.org',
    RDataPath = 'record/7959787/files/stomics_mm_brain.rds'
  ),
  c(
    Title = 'cosmx_hs_nsclc',
    Description = 'NanoString CosMx dataset of FFPE tissues from 5 non-small cell lung cancer (NSCLC) patients. Lung 5 has 2 serial sections while Lung 9 has 2 non-serial (different location of the tumour) sections from the same patient. Raw transcript level data is provided with region annotations for each detection obtained using expert annotation from a biologist.',
    BiocVersion = as.numeric(as.character(BiocManager::version())) + 0.01,
    Genome = NA,
    SourceType = 'CSV',
    SourceUrl = 'https://nanostring.com/products/cosmx-spatial-molecular-imager/nsclc-ffpe-dataset',
    SourceVersion = '2023-03-09',
    Species = 'Homo sapiens',
    TaxonomyId = '9606',
    Coordinate_1_based = TRUE,
    DataProvider = 'NanoString',
    Maintainer = 'Dharmesh D. Bhuva <dharmesh.bhuva@adelaide.edu.au>',
    RDataClass = 'data.frame',
    DispatchClass = 'Rds',
    Location_Prefix = 'https://zenodo.org',
    RDataPath = 'record/7959787/files/cosmx_hs_nsclc.rds'
  ),
  c(
    Title = 'xenium_hs_breast_addon',
    Description = '10x Xenium dataset of an intraductal carcinoma (IDC) and intralobular carcinoma (ILC) breast cancers. Raw transcript level data is provided with region annotations for each detection obtained using by annotating the accompanying histology image (H&E).',
    BiocVersion = as.numeric(as.character(BiocManager::version())) + 0.01,
    Genome = NA,
    SourceType = 'CSV',
    SourceUrl = 'https://www.10xgenomics.com/datasets/ffpe-human-breast-with-custom-add-on-panel-1-standard',
    SourceVersion = '2024-01-31',
    Species = 'Homo sapiens',
    TaxonomyId = '9606',
    Coordinate_1_based = TRUE,
    DataProvider = '10x',
    Maintainer = 'Dharmesh D. Bhuva <dharmesh.bhuva@adelaide.edu.au>',
    RDataClass = 'data.frame',
    DispatchClass = 'Rds',
    Location_Prefix = 'https://zenodo.org',
    RDataPath = "record/10516814/files/xenium_hs_breast_addon.rds"
  )
)

#write and test metadata file
write.csv(meta, file = 'inst/extdata/metadata.csv', row.names = FALSE)
ExperimentHubData::makeExperimentHubMetadata('../SubcellularSpatialData')
