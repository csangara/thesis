#### HELPER FUNCTIONS #####

library(Seurat)
library(synthvisium)
library(dplyr)
library(Biobase)
library(stringr)

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


######## PLOT PROPS #########
# library(ggplot2)
# 
# n = nrow(decon_mtrx)
# knownP = synthetic_visium_data$relative_spot_composition[1:n,1:23]
# predP = decon_mtrx[1:n,1:23]
# new_pred = data.frame("pred" = c(t(predP)))
# new_pred$known = as.numeric(c(t(knownP)))
# new_pred$celltype = rep(colnames(knownP), n)
# 
# ggplot(new_pred, aes(x=known, y=pred, shape=celltype)) +
#   geom_abline(slope=1, intercept=0, linetype=2, colour="gray20") +
#   scale_shape_manual(values = 1:23) + geom_point(size=1) +
#   labs(x="Known Proportions", y="Predicted Proportions")
