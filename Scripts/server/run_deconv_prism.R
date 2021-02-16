setwd("~/thesis/")
path <- "~/thesis/"
library(Seurat)
library(synthvisium)
library(dplyr)
source("Scripts/helperFunctions.R")

recommended_dataset_types = c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap",
                              "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
                              "artificial_partially_dominant_celltype_diverse", "artificial_dominant_rare_celltype_diverse",
                              "artificial_regional_rare_celltype_diverse")

datasets <- c('allen_cortex_dwn', 'brain_cortex_generation', 'cerebellum_cell_generation', 'cerebellum_nucleus_generation',
              'hippocampus_generation', 'kidney_generation', 'pbmc_generation', 'scc_p5_generation')

######## GENERATING NEEDED DATA ############
# Create synthetic data from scRNA data and save as RDS
for (dataset_type in possible_dataset_types){
  print(dataset_type)
  createSynthvisiumRDS("rds/allen_cortex_dwn_original.rds", dataset_type,
                       output_folder="Data/synthetic_datasets/allen_cortex_dwn/")
}

dataset <- datasets[2]
scrna_dir <- "/group/irc/shared/synthetic_visium/generation/"
scrna_path <- paste0(scrna_dir, dataset, ".rds")
dir.create(paste0("results/", dataset))

######### SPOTLIGHT ######### 

library(SPOTlight)
dir.create(paste0("results/", dataset, "/spotlight/"))
# Extract the top marker genes from each cluster
seurat_obj_scRNA <- readRDS(scrna_path)
seurat_obj_scRNA <- preprocessSeurat(seurat_obj_scRNA)
Idents(object = seurat_obj_scRNA) <- seurat_obj_scRNA@meta.data$subclass
cluster_markers_all <- FindAllMarkers(object = seurat_obj_scRNA, assay = "SCT", slot = "data", verbose = TRUE, 
                                      only.pos = TRUE, logfc.threshold = 1, min.pct = 0.9)

for (dataset_type in possible_dataset_types){
  set.seed(123)
  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", dataset, "_", dataset_type, "_synthvisium.rds"))
  seurat_obj_visium <- createSeuratFromCounts(synthetic_visium_data$counts)
  start_time <- Sys.time()
  spotlight_deconv <- spotlight_deconvolution(se_sc = seurat_obj_scRNA, counts_spatial = seurat_obj_visium@assays$Spatial@counts,
                                              clust_vr = "subclass", cluster_markers = cluster_markers_all, cl_n = 50,
                                              hvg = 3000, ntop = NULL, transf = "uv", method = "nsNMF", min_cont = 0.09)
  end_time <- Sys.time()
  decon_mtrx <- spotlight_deconv[[2]]
  
  res = cor(t(synthetic_visium_data$relative_spot_composition[,1:23]), t(decon_mtrx[,1:23]))
  print(mean(diag(res), na.rm=TRUE))
  
  #res2 = cor(synthetic_visium_data$relative_spot_composition[,1:23], decon_mtrx[,1:23], use="complete.obs")
  #mean(diag(res2), na.rm=TRUE)
  
  saveRDS(spotlight_deconv, paste0("results/", dataset, "/spotlight/", dataset, "_", dataset_type, "_spotlight.rds"))
}

######### RCTD ##########

library(RCTD)
library(Matrix)
dir.create(paste0("results/", dataset, "/RCTD/"))
seurat_obj_scRNA <- readRDS(scrna_path)
seurat_obj_scRNA@meta.data$liger_ident_coarse <- factor(seurat_obj_scRNA@meta.data$subclass)
seurat_obj_scRNA@meta.data$nUMI <- colSums(seurat_obj_scRNA@assays$RNA@counts)

for (dataset_type in possible_dataset_types[2:length(possible_dataset_types)]){
  
  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", dataset, "_", dataset_type, "_synthvisium.rds"))
  spatialRNA_obj_visium <- RCTD:::SpatialRNA(counts=as(as(synthetic_visium_data$counts,"matrix"),"dgCMatrix"))
  
  start_time <- Sys.time()
  RCTD_deconv <- create.RCTD(spatialRNA_obj_visium, seurat_obj_scRNA, max_cores = 4, CELL_MIN_INSTANCE = 5)
  RCTD_deconv <- run.RCTD(RCTD_deconv, doublet_mode = FALSE)
  end_time <- Sys.time()
  res = as.matrix(sweep(RCTD_deconv@results$weights, 1, rowSums(RCTD_deconv@results$weights), '/'))
  
  saveRDS(res, paste0("results/", dataset, "/RCTD/", dataset, "_", dataset_type, "_RCTD.rds"))
}