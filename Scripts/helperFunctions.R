#### HELPER FUNCTIONS #####

library(Seurat)
library(dplyr)
library(stringr)
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

# Creates synthetic visium data given a reference scRNA-seq data and dataset type
# Saves it to subfolder synthvisium_spatial
createSynthvisiumRDS <- function(inputscRNA_rds, dataset_type, output_folder="synthvisium_spatial/"){
  seurat_obj_scRNA =  readRDS(inputscRNA_rds)

  # Create synthetic visium data from scRNA data
  synthetic_visium_data = synthvisium::generate_synthetic_visium(seurat_obj = seurat_obj_scRNA, dataset_type = dataset_type, 
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
createSeuratFromCounts <- function(counts_data, PP=TRUE, metadata=NULL){
  seurat_obj_visium = CreateSeuratObject(counts = counts_data, assay = "Spatial", meta.data=metadata)
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
  
  sc.pdata <- new("AnnotatedDataFrame", data=sc.pheno)
  sc.data <- as.matrix(GetAssayData(seurat_object)[,names(Idents(seurat_object)),drop=F])
  sc.eset <- Biobase::ExpressionSet(assayData=sc.data, phenoData=sc.pdata)
  return(sc.eset)
}

# Convert Seurat object to Loom
convertSeuratRDSToLoom <- function(input_path, output_path=NULL,isSeurat=TRUE, raw=TRUE){
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

# Convert Seurat object to h5ad
convertSeuratRDSToh5ad <- function(input_path, output_path=NULL, isSeurat=TRUE, raw=TRUE, update=FALSE){
  seurat_obj <- readRDS(input_path)
  if (!isSeurat){
    seurat_obj <- CreateSeuratObject(counts = seurat_obj$counts, assay = "Spatial")
  }
  if (raw) { DefaultAssay(seurat_obj) <- "RNA" }
  if (update){
    seurat_obj <- UpdateSeuratObject(seurat_obj)
    seurat_obj <- CreateSeuratObject(counts = GetAssayData(seurat_obj), assay = DefaultAssay(seurat_obj))
  }
  
  file_name <- tools::file_path_sans_ext(input_path)
  if (is.null(output_path)) {output_path <- dirname(input_path)}
  file_name <- stringr::str_split(basename(input_path), "\\.")[[1]][1]
  if (file.exists(paste0(output_path, "/", file_name, ".h5ad"))){return ("h5ad file exists") }
  SeuratDisk::SaveH5Seurat(seurat_obj, filename = paste0(output_path, "/", file_name, ".h5Seurat"))
  SeuratDisk::Convert(paste0(output_path, "/", file_name, ".h5Seurat"), dest = "h5ad")
  file.remove(paste0(output_path, "/", file_name, ".h5Seurat"))
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
  rds_methods <- c("spotlight", "music", "RCTD", "spotlight_optim")
  csv_methods <- c("cell2location")
  tsv_methods <- c("stereoscope", "cibersort", "stereoscope_optim")
  
  for (method in methods){
    file_name <- paste0(result_path, method, "/", dataset, "_", dataset_type, "_", str_split(method, "_")[[1]][1])
    
    if (method %in% rds_methods){
      temp_deconv <- readRDS(paste0(file_name, ".rds"))
      
      if (grepl("spotlight", method)){temp_deconv <- temp_deconv[[2]]}
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
      
      if (grepl("stereoscope", method)){
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

# Get region composition from synthvisium dataset
getregionComp <- function(synthetic_visium_data, verbose=TRUE,
                           plotUMAP=FALSE, labelUMAP=TRUE){
  # Create seurat obj and UMAP
  spot_comp <- synthetic_visium_data$relative_spot_composition
  
  if (plotUMAP){
    seurat_obj_visium <- createSeuratFromCounts(synthetic_visium_data$counts)
    p <- DimPlot(seurat_obj_visium, reduction = "umap", label = FALSE, group.by = "orig.ident")
  }
  
  region_comp = c()
  for (region in paste0("priorregion", 1:5)){
    temp_spot_comp <- spot_comp[spot_comp$region==region,1:(ncol(spot_comp)-2)]  
    mean_comp <- apply(temp_spot_comp, 2, mean)
    mean_comp <- mean_comp[mean_comp != 0]
    comp_text <- paste(names(mean_comp), round(mean_comp, 2), sep = ":", collapse = "; ")
    region_comp[region] <- comp_text
    
    if (verbose){
      print(region)
      print(comp_text)
    }
    
    # sd_comp <- apply(temp_spot_comp, 2, sd)
    # sd_comp <- sd_comp[sd_comp != 0]
    # print(paste(names(sd_comp), round(sd_comp, 2), sep = ":", collapse = "; "))
    
    if (plotUMAP & labelUMAP){
      umap_coords <-seurat_obj_visium@reductions$umap@cell.embeddings[seurat_obj_visium$orig.ident==region,]
      median_umap <- apply(umap_coords, 2, median)
      p <- p + annotate("text", label=comp_text, x=median_umap[1], y=median_umap[2])
    }
    
  }
  if (plotUMAP) {print(p + ggtitle(dataset_type))}
  return (region_comp)
}
