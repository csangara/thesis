library(SeuratDisk)
library(Seurat)
# Convert Seurat object to Loom
convertSeuratRDSToLoom <- function(input_path, output_path=NULL, rep="", isSeurat=TRUE, raw=TRUE){
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

convertSeuratRDSToh5ad <- function(input_path, output_path=NULL, isSeurat=TRUE, raw=TRUE, update=FALSE){
  seurat_obj <- readRDS(input_path)
  if (!isSeurat){
    seurat_obj <- CreateSeuratObject(counts = seurat_obj$counts, assay = "Spatial")
  }
  if (raw) { DefaultAssay(seurat_obj) <- "RNA" }
  if (update){
    seurat_obj <- UpdateSeuratObject(seurat_obj)
    seurat_obj <- CreateSeuratObject(counts = GetAssayData(seurat_obj), assay = DefaultAssay(seurat_obj),
                                     meta.data=seurat_obj@meta.data)
  }
  
  file_name <- tools::file_path_sans_ext(input_path)
  if (is.null(output_path)) {output_path <- dirname(input_path)}
  file_name <- stringr::str_split(basename(input_path), "\\.")[[1]][1]
  if (file.exists(paste0(output_path, "/", file_name, ".h5ad"))){return ("h5ad file exists") }
  SaveH5Seurat(seurat_obj, filename = paste0(output_path, "/", file_name, ".h5Seurat"))
  Convert(paste0(output_path, "/", file_name, ".h5Seurat"), dest = "h5ad")
  file.remove(paste0(output_path, "/", file_name, ".h5Seurat"))
}

path <- "/home/chananchidas/data/"
repl <- "rep6/"
datasets <- c("brain_cortex_generation", "cerebellum_cell_generation",
              "cerebellum_nucleus_generation", "hippocampus_generation",
              "kidney_generation", "pbmc_generation", "scc_p5_generation")

dataset_types = c("artificial_uniform_distinct", "artificial_diverse_distinct", 
                  "artificial_uniform_overlap", "artificial_diverse_overlap",
                  "artificial_dominant_celltype_diverse",
                  "artificial_partially_dominant_celltype_diverse",
                  "artificial_dominant_rare_celltype_diverse",
                  "artificial_regional_rare_celltype_diverse")
# for (i in 7:10){
#   repl <- paste0("rep", i, "/")
#   for (dataset in datasets){
#     for (dataset_type in dataset_types){
#       print(paste(dataset, dataset_type))
#       convertSeuratRDSToLoom(paste0(path, dataset, "/", repl, dataset, "_",
#                                     dataset_type, "_synthvisium.rds"),
#                              isSeurat=FALSE, raw=FALSE)
#       convertSeuratRDSToh5ad(paste0(path, dataset, "/", repl, dataset, "_",
#                                     dataset_type, "_synthvisium.rds"),
#                              isSeurat=FALSE, raw=FALSE)
#     }
#   }
# }

# generation files
# path <- "/group/irc/shared/synthetic_visium/generation/"
# for (dataset in datasets){
#     convertSeuratRDSToh5ad(paste0(path, dataset, ".rds"), output_path="data/generation_h5ad/", update=TRUE)
# }

# reference files
# test_set <- c('brain_cortex_test.rds', 'cerebellum_cell_test.rds', 'cerebellum_nucleus_test.rds',
#               'hippocampus_test.rds', 'kidney_test.rds',
#               'pbmc_test.rds', 'scc_p5_test.rds')

# for (dataset in test_set[2:3]){
#   convertSeuratRDSToLoom(paste0(path, "test_set/", dataset))
# }

# convertSeuratRDSToh5ad("/group/irc/shared/synthetic_visium/raw_data/scRNAseq/seurat_obj_scrnaseq_kidney_filtered.rds",
#   output_path="data/", update=TRUE)

convertSeuratRDSToh5ad("/group/irc/shared/synthetic_visium/raw_data/brain_cortex/scRNAseq/seurat_obj_scrnaseq_cortex_filtered.rds",
   output_path="thesis/downsampling/brain_cortex/", update=TRUE)


