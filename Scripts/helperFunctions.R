#### HELPER FUNCTIONS #####

library(Seurat)
library(synthvisium)
library(dplyr)
library(Biobase)
library(stringr)
library(loomR)
library(SeuratDisk)
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

# Creates synthetic visium data given a reference scRNA-seq data and dataset type
# Saves it to subfolder synthvisium_spatial
createSynthvisiumRDS <- function(inputscRNA_rds, dataset_type, output_folder="synthvisium_spatial/"){
  seurat_obj_scRNA =  readRDS(inputscRNA_rds)

  # Create synthetic visium data from scRNA data
  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj_scRNA, dataset_type = dataset_type, 
                                                    clust_var = "subclass", region_var = "brain_subregion" , n_regions = NULL,
                                                    n_spots_min = 50, n_spots_max = 200, visium_mean = 20000, visium_sd = 5000)
  
  directory = dirname(inputscRNA_rds)
  inputscRNA_name = str_split(basename(inputscRNA_rds), "\\.")[[1]][1]
  
  saveRDS(synthetic_visium_data, paste0(directory, "/", output_folder, inputscRNA_name, "_", dataset_type, "_synthvisium.rds"))
  print(paste0("Dataset saved at ", directory, "/", output_folder, inputscRNA_name, "_", dataset_type, "_synthvisium.rds"))
}

# Create seurat object from count data (meant to be used with synthetic visium data) and perform normalization,
# dimensionality reduction, and clustering
createAndPPSeuratFromVisium <- function(counts_data){
  seurat_obj_visium = CreateSeuratObject(counts = counts_data, min.cells = 2, min.features = 200, assay = "Spatial")
  seurat_obj_visium = SCTransform(seurat_obj_visium, assay = "Spatial", verbose = FALSE)
  seurat_obj_visium = RunPCA(seurat_obj_visium, assay = "SCT", verbose = FALSE)
  seurat_obj_visium = RunTSNE(seurat_obj_visium, reduction = "pca", dims = 1:30)
  seurat_obj_visium = RunUMAP(seurat_obj_visium, reduction = "pca", dims = 1:30)
  seurat_obj_visium = FindNeighbors(seurat_obj_visium, reduction = "pca", dims = 1:30)
  seurat_obj_visium = FindClusters(seurat_obj_visium, verbose = FALSE, resolution = 0.5)
  
  return(seurat_obj_visium)
}

# Create ExpressionSet object from seurat object, assuming the Ident of the object is the cell type
SeuratToExprSet <- function(seurat_object){
  sc.pheno <- data.frame(check.names=F, check.rows=F,
                         stringsAsFactors=F,
                         row.names=names(Idents(seurat_object)),
                         samples=names(Idents(seurat_object)))
  if (!is.null(seurat_object@meta.data$celltype)){
    sc.pheno$celltype = seurat_object@meta.data$celltype
  }
  
  sc.pdata <- new("AnnotatedDataFrame",
                  data=sc.pheno)
  sc.data <- as.matrix(GetAssayData(seurat_object)[,names(Idents(seurat_object)),drop=F])
  sc.eset <- ExpressionSet(assayData=sc.data, phenoData=sc.pdata)
  return(sc.eset)
}

# Write file in a CIBERSORT-compatible format given a Seurat object, by default use SCT column
# Other arguments are "raw" and "CPM"
WriteFileCS <- function(seurat_obj, output_file, assay="SCT"){
  if (assay=="SCT"){
    count_matrix <-  as.matrix(GetAssayData(seurat_obj)) # SCTransformed
  } else if (assay=="CPM"){
    count_matrix <- as.matrix(GetAssayData(seurat_obj[["RNA"]]))
    count_matrix <- count_matrix/sum(count_matrix)*10^6
  } else if (assay=="raw"){
    count_matrix <- as.matrix(GetAssayData(seurat_obj[["RNA"]]))
  } else {
    return(NULL)
  }
  
  colnames(count_matrix) <- seurat_obj@meta.data$subclass
  
  # Write to CIBERSORT-compatible format
  write.table("Gene", file = output_file, sep = "\t",eol="\t", append = FALSE, quote = FALSE, row.names=FALSE, col.names = FALSE)
  write.table(SCT_matrix, file = output_file, sep = "\t", append = TRUE, quote = FALSE, row.names = TRUE, col.names = TRUE)
  print(paste0("Wrote file successfully at ", output_file))
}

reduceSpotsCS <- function(input_path, no_spots){
  input_file = read.table(input_path, sep="\t", row.names=1, header=TRUE)
  # Randomly sample no_spots
  reduced_cols <- sample(seq(1, ncol(input_file)), no_spots)
  input_file <- input_file[,reduced_cols]
  
  file_name <- tools::file_path_sans_ext(input_path); file_ext <- tools::file_ext(input_path)
  output_file <- paste0(file_name, "_", no_spots, ".", file_ext)
  
  write.table("GeneSymbol", file = output_file, sep = "\t",eol="\t", append = FALSE, quote = FALSE, row.names=FALSE, col.names = FALSE)
  write.table(input_file, file = output_file, sep = "\t", append = TRUE, quote = FALSE, row.names = TRUE, col.names = TRUE)
  print(paste0("Wrote file successfully at ", output_file))
}

convertSeuratRDSToLoom <- function(input_path, createSeuratFromRDS=FALSE){
  seurat_obj =  readRDS(input_path)
  if (createSeuratFromRDS){ seurat_obj <- createAndPPSeuratFromVisium(seurat_obj$counts) }
  file_name <- tools::file_path_sans_ext(input_path)
  as.loom(seurat_obj, filename = paste0(file_name, ".loom"))
}

convertSeuratRDSToh5ad <- function(input_path,createSeuratFromRDS=FALSE){
  seurat_obj =  readRDS(input_path)
  file_name <- tools::file_path_sans_ext(input_path)
  if (createSeuratFromRDS){ seurat_obj <- createAndPPSeuratFromVisium(seurat_obj$counts) }
  SaveH5Seurat(seurat_obj, filename = paste0(file_name, ".h5Seurat"))
  Convert(paste0(file_name, ".h5Seurat"), dest = "h5ad")
}
