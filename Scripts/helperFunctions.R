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

# Perform normalization, dimensionality reduction, and clustering on a seurat object
preprocessSeurat <- function(seurat_obj, assay="RNA"){
  seurat_obj = SCTransform(seurat_obj, assay = assay, verbose = FALSE)
  seurat_obj = RunPCA(seurat_obj, assay = "SCT", verbose = FALSE)
  seurat_obj = RunTSNE(seurat_obj, reduction = "pca", dims = 1:30, check_duplicates = FALSE)
  seurat_obj = RunUMAP(seurat_obj, reduction = "pca", dims = 1:30)
  seurat_obj = FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
  seurat_obj = FindClusters(seurat_obj, verbose = FALSE, resolution = 0.5)  
  return(seurat_obj)
}

# Create seurat object from count data (meant to be used with synthetic visium data)
createSeuratFromCounts <- function(counts_data, PP=TRUE){
  seurat_obj_visium = CreateSeuratObject(counts = counts_data, min.cells = 2, min.features = 200, assay = "Spatial")
  if (PP){ seurat_obj_visium <- preprocessSeurat(seurat_obj_visium, "Spatial") }
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

# Convert Seurat object to Loom
convertSeuratRDSToLoom <- function(input_path, raw=TRUE, isSeurat=TRUE, PP=FALSE){
  seurat_obj =  readRDS(input_path)
  if (!isSeurat){ seurat_obj <- createAndPPSeuratFromCounts(seurat_obj$counts, PP=PP) }
  if (raw) { DefaultAssay(seurat_obj) <- "RNA"}
  file_name <- tools::file_path_sans_ext(input_path)
  as.loom(seurat_obj, filename = paste0(file_name, ".loom"))
}

# Convert Seurat object to h5ad (pp = preprocess and sctransform, if FALSE, raw counts will be saved)
convertSeuratRDSToh5ad <- function(input_path, raw=TRUE, isSeurat=TRUE, PP=FALSE){
  seurat_obj = readRDS(input_path)
  file_name <- tools::file_path_sans_ext(input_path)
  if (!isSeurat){ seurat_obj <- createAndPPSeuratFromCounts(seurat_obj$counts, PP=PP) }
  if (raw) { DefaultAssay(seurat_obj) <- "RNA"}
  SaveH5Seurat(seurat_obj, filename = paste0(file_name, ".h5Seurat"))
  Convert(paste0(file_name, ".h5Seurat"), dest = "h5ad")
  file.remove(paste0(file_name, ".h5Seurat"))
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

# Downsample mixture matrix to no_spots (new file written in CIBERSORT-compatible format)
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

# Create a list with each element containing the proportion matrix returned from each method
createDeconvResultList <- function(methods, celltypes, result_path, dataset){
  results <- list()
  
  # Divide methods depending on their output file
  rds_methods <- c("spotlight", "music", "RCTD")
  csv_methods <- c("cell2location")
  tsv_methods <- c("stereoscope", "cibersort")
  
  for (method in methods){
    file_name <- paste0(result_path, method, "/", dataset, "_", dataset_type, "_", method)
    
    if (method %in% rds_methods){
      temp_deconv <- readRDS(paste0(file_name, ".rds"))
      
      if (method == "spotlight"){temp_deconv <- temp_deconv[[2]]}
      colnames(temp_deconv) <- str_replace_all(colnames(temp_deconv), "[/ ]", ".")
      temp_deconv <- temp_deconv[, match(celltypes, colnames(temp_deconv))]
      
    } else if (method %in% csv_methods){
      temp_deconv <- read.csv(paste0(file_name, ".csv"), row.names=1)
      
      if (method == "cell2location"){
        temp_deconv <- temp_deconv/rowSums(temp_deconv) # So everything sums to one
        colnames(temp_deconv) <- str_replace(colnames(temp_deconv), "q05_spot_factors", "")
        temp_deconv <- temp_deconv[, match(celltypes, colnames(temp_deconv))]
      } else {
        colnames(temp_deconv)[1:length(celltypes)] <- celltypes
      }
      
    } else if (method %in% tsv_methods){
      temp_deconv = read.table(paste0(paste0(file_name, ".tsv")), sep="\t", row.names=1, header=TRUE)
      
      if (method == "stereoscope"){
        temp_deconv <- temp_deconv[, match(celltypes, colnames(temp_deconv))]
      } else {
        colnames(temp_deconv)[1:length(celltypes)] <- celltypes
      }
    }
    results[[method]] <- temp_deconv
  }
  return(results)
}

# Get the confusion matrix given ground truth and deconvoluted result
getConfusionMatrix <- function(known_props, test_props){
  test_props <- round(test_props, 2)
  tp <- 0; tn <- 0; fp <- 0; fn <- 0
  missing_rows <- which(rowSums(is.na(known_props)) > 0)
  
  if (length(missing_rows) > 0){
    test_props <- test_props[-missing_rows,]
    known_props <- known_props[-missing_rows,]
  }
  for (i in 1:nrow(known_props)){
    for (j in 1:ncol(known_props)){
      if (known_props[i, j] > 0 & test_props[i, j] > 0){
        tp <- tp + 1
      } else if (known_props[i, j] == 0 & test_props[i, j] == 0){
        tn <- tn + 1
      } else if (known_props[i, j] > 0 & test_props[i, j] == 0){
        fn <- fn + 1
      } else if (known_props[i, j] == 0 & test_props[i, j] > 0){
        fp <- fp + 1
      }
    }
  }
  return(list(tp=tp, tn=tn, fn=fn, fp=fp))
}

