library(SeuratDisk)
library(Seurat)
# Convert Seurat object to Loom
convertSeuratRDSToLoom <- function(input_path, output_path=NULL, isSeurat=TRUE, raw=TRUE){
  seurat_obj <- readRDS(input_path)
  if (!isSeurat) {
    seurat_obj <- Seurat::CreateSeuratObject(counts = seurat_obj$counts, assay = "Spatial")
  }
  if (raw) {DefaultAssay(seurat_obj) <- "RNA"}
  file_name <- tools::file_path_sans_ext(input_path)
  if (is.null(output_path)) {output_path <- dirname(input_path)}
  file_name <- stringr::str_split(basename(input_path), "\\.")[[1]][1]
  
  SeuratDisk::as.loom(seurat_obj, filename = paste0(output_path, "/", file_name, ".loom"))
}

convertSeuratRDSToh5ad <- function(input_path, output_path=NULL, isSeurat=TRUE,
                                   raw=TRUE, PP=FALSE){
  seurat_obj = readRDS(input_path)
  
  if (!isSeurat){
    seurat_obj <- createSeuratFromCounts(seurat_obj$counts, PP=PP)
  }
  if (raw) { DefaultAssay(seurat_obj) <- "RNA"}
  file_name <- tools::file_path_sans_ext(input_path)
  if (is.null(output_path)) {output_path <- dirname(input_path)}
  SaveH5Seurat(seurat_obj, filename = paste0(output_path, "/", file_name, ".h5Seurat"))
  Convert(paste0(output_path, "/", file_name, ".h5Seurat"), dest = "h5ad")
  file.remove(paste0(output_path, "/", file_name, ".h5Seurat"))
}

path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"
repl <- "rep1/"
datasets <- c("brain_cortex_generation.rds", "cerebellum_cell_generation.rds",
              "cerebellum_nucleus_generation.rds", "hippocampus_generation.rds",
              "kidney_generation.rds", "pbmc_generation.rds", "scc_p5_generation.rds")

dataset_types = c("artificial_uniform_distinct", "artificial_diverse_distinct", 
                  "artificial_uniform_overlap", "artificial_diverse_overlap",
                  "artificial_dominant_celltype_diverse",
                  "artificial_partially_dominant_celltype_diverse",
                  "artificial_dominant_rare_celltype_diverse",
                  "artificial_regional_rare_celltype_diverse")

for (dataset in datasets){
  dataset = stringr::str_split(dataset, "\\.")[[1]][1]
  for (dataset_type in dataset_types){
    print(paste(dataset, dataset_type))
    convertSeuratRDSToLoom(paste0(path, dataset, "/", repl, dataset, "_",
                                  dataset_type, "_synthvisium.rds"),
                           isSeurat=FALSE, raw=FALSE)
  }
}

# reference files
test_set <- c('brain_cortex_test.rds', 'cerebellum_cell_test.rds', 'cerebellum_nucleus_test.rds',
              'hippocampus_test.rds', 'kidney_test.rds',
              'pbmc_test.rds', 'scc_p5_test.rds')

for (dataset in test_set[2:3]){
  convertSeuratRDSToLoom(paste0(path, "test_set/", dataset))
}

